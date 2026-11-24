import os
import sys
import argparse
import numpy as np
from scipy import stats
from concurrent.futures import ProcessPoolExecutor, wait, FIRST_COMPLETED
from collections import defaultdict

import numpy as np
import pandas as pd
from tqdm import tqdm

from pastri import SMaHTid
from pastri.plume import perform_clustering, coverage_filter
from pastri.qdpi_io import stream_qdpi


# Don't like separating these from the perform_clustering outputs
GERM_COLS = ["chrom", "start", "end", "donor", "hap", "nclust", "ngerm", "nassign", 
             "coverage", "bandwidth", "mean", "mid", "lower", "upper", "norm_entropy", "state"]
UNST_COLS = ["chrom", "start", "end", "donor", "protocol", "hap", "delta_min",
             "delta_mean", "delta_mid", "delta_max", "spread", "alt_reads", 
             "coverage", "vaf"]

def calc_germ(reads, grp_cols, length_col='length'):
    """
    Calculate all the GERM_COLS for a set of reads
    The reads must already be subset. So really this just works no a set of lengths.
    coverage mean mid lower upper norm_entropy
    And you can subset/fillna for the others.
    """
    return (
        reads.groupby(grp_cols, as_index=False).agg(
            coverage=(length_col, 'count'),
            mean=(length_col, 'mean'),
            mid=(length_col, 'median'),
            lower=(length_col, 'min'),
            upper=(length_col, 'max'),
            norm_entropy=(length_col, normalized_entropy)
        )
        .round({'mean': 1, 'norm_entropy':4})
        .astype({'mid': int})
    )
        
    
def bounded_map(executor, fn, arg_iter, max_pending):
    """
    Lazily submit fn(*args) for each args tuple in arg_iter, keeping at most
    `max_pending` futures in flight. Yields results as they complete
    (unordered). arg_iter is only advanced when there's room for more work.
    """
    pending = set()
    for args in arg_iter:
        pending.add(executor.submit(fn, *args))
        if len(pending) >= max_pending:
            done, pending = wait(pending, return_when=FIRST_COMPLETED)
            for fut in done:
                yield fut.result()
    # Source exhausted; drain what's left
    while pending:
        done, pending = wait(pending, return_when=FIRST_COMPLETED)
        for fut in done:
            yield fut.result()

def parse_args(args):
    parser = argparse.ArgumentParser(prog="scan", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_qdpi", type=str,
                        help="File with paths of input qdpi files")
    parser.add_argument("-s", "--subset", type=str,
                        help="Bed file of loci to analyze")
    parser.add_argument("-o", "--output", default="output",
                        help="Output prefix (%(default)s)")
    parser.add_argument("-b", "--bandwidth", type=int, default=10,
                        help="Bandwidth for counting out_range read (%(default)s)")
    #parser.add_argument("--min-cov", type=int, default=20,
    #                    help="Minimum coverage over the locus (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (%(default)s)")

    args = parser.parse_args(args)
    return args

def normalized_entropy(lengths: np.ndarray, bins: int = 50) -> float:
    counts, _ = np.histogram(lengths, bins=bins)
    probs = counts / counts.sum()
    probs = probs[probs > 0]
    ent = -(probs * np.log(probs)).sum()
    return ent / np.log(bins)  # 1.0 = perfectly flat across bins

def run_analysis(locus, h1_parts, h2_parts, smhtids, args):
    chrom, start, end = locus
    span = int(end) - int(start)
    
    # Turn to DataFrame compatible with perform_clustering
    h1_counts = [x.size for x in h1_parts]
    h2_counts = [x.size for x in h2_parts]

    # 2. Repeat sample IDs based on the counts of each sub-array
    h1_smhtids = np.repeat(smhtids, h1_counts)
    h2_smhtids = np.repeat(smhtids, h2_counts)

    # 3. Concatenate arrays
    all_smhtids = np.concatenate([h1_smhtids, h2_smhtids])
    all_haps = np.repeat([1, 2], [sum(h1_counts), sum(h2_counts)])

    # Handle empty lists gracefully before concatenation
    h1_flat = np.concatenate(h1_parts) if sum(h1_counts) > 0 else np.array([], dtype=int)
    h2_flat = np.concatenate(h2_parts) if sum(h2_counts) > 0 else np.array([], dtype=int)
    all_lengths = np.concatenate([h1_flat, h2_flat]) + span

    # 4. Construct DataFrame instantly
    reads = pd.DataFrame({'smhtid': all_smhtids,
                         'hap': all_haps,
                         'length': all_lengths
    })
    reads['donor'] = reads['smhtid'].apply(lambda x: x.donor)
    reads['protocol'] = reads['smhtid'].apply(lambda x: x.protocol)

    germ = calc_germ(reads, ['donor', 'hap'])
    germ.insert(0, 'end', locus[2])
    germ.insert(0, 'start', locus[1])
    germ.insert(0, 'chrom', locus[0])
    
    tissue = calc_germ(reads, ['donor', 'protocol', 'hap'])
    tissue.insert(0, 'end', locus[2])
    tissue.insert(0, 'start', locus[1])
    tissue.insert(0, 'chrom', locus[0])
    
    # Does any donor-tissue have at least min_reads above min_vaf outside germ median±bandwidth.
    # 1. Bring the donor/hap mean onto each read
    merged = reads[['donor', 'protocol', 'hap', 'length']].merge(germ[['donor', 'hap', 'mid']], 
                                                                 on=['donor', 'hap'], 
                                                                 how='left')
    
    # 2. Flag reads
    merged['out_range'] = ~merged['length'].between(merged['mid'] - args.bandwidth, merged['mid'] + args.bandwidth)
    
    # 3. Group and sum
    summary = merged.groupby(['donor', 'protocol', 'hap'], as_index=False).agg(
        total_reads=('length', 'count'),
        out_range=('out_range', 'sum')
    )
    summary['vaf'] = (summary['out_range'] / summary['total_reads']).round(4)
    
    tissue = tissue.merge(
        summary[['donor', 'protocol', 'hap', 'out_range', 'vaf']],
        on=['donor', 'protocol', 'hap'],
        how='left'
    )
    return germ, tissue


def scan_main(args):
    args = parse_args(args)

    to_analyze = None
    if args.subset:
        to_analyze = set()
        with open(args.subset, 'r') as fh:
            for line in fh:
                to_analyze.add(tuple(line.strip().split('\t')[:3]))

    inputs = open(args.in_qdpi, 'r').read().strip().split('\n')
    smhtids = [SMaHTid.from_re(os.path.basename(_)) for _ in inputs]

    germ_out = open(f'{args.output}.germline.tsv', 'w')
    f_germ = True # For header writing
    unst_out = open(f'{args.output}.tissues.tsv', 'w')
    f_filt = True # For header writing

    # Debug
    #for locus, h1_parts, h2_parts in stream_qdpi(inputs, to_analyze):
        #data, germ, filtered = run_analysis(locus, h1_parts, h2_parts, smhtids, args)
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        tasks = (
            (locus, h1_parts, h2_parts, smhtids, args)
            for locus, h1_parts, h2_parts in stream_qdpi(inputs, to_analyze)
        )
        results = bounded_map(executor, run_analysis, tasks, max_pending=args.threads * 4)
        
        for germ, filtered in tqdm(results):
            if germ is None:
                continue
            
            germ.to_csv(germ_out, sep='\t', index=False, header=f_germ)
            f_germ = False

            if filtered is not None:
                filtered.to_csv(unst_out, sep='\t', index=False, header=f_filt)
                f_filt = False

if __name__ == '__main__':
    scan_main(sys.argv[1:])
