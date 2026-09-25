import os
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import numpy as np
from scipy import stats

import numpy as np
import pandas as pd
from tqdm import tqdm

from pastri import SMaHTid
from pastri.plume import perform_clustering, coverage_filter
from pastri.qdpi_io import stream_qdpi

def parse_args(args):
    parser = argparse.ArgumentParser(prog="scan", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_qdpi", type=str,
                        help="File with paths of input qdpi files")
    parser.add_argument("-s", "--subset", type=str,
                        help="Bed file of loci to analyze")
    parser.add_argument("-o", "--output", default="output",
                        help="Output prefix (%(default)s)")
    parser.add_argument("-r", "--min-reads-donor", type=int, default=3,
                        help="Minimum number of non-germ reads required in donor (%(default)s)")
    parser.add_argument("-R", "--min-reads-tissue", type=int, default=3,
                        help="Minimum number of non-germ reads required in tissue (%(default)s)")
    parser.add_argument("-b", "--min-bandwidth", type=int, default=10,
                        help="Minimum MeanShift clustering bandwidth (%(default)s)")
    parser.add_argument("-g", "--germ-vaf", type=float, default=0.8,
                        help="Minimum fraction of reads to collect germline cluster (%(default)s)")
    parser.add_argument("-q", "--germ-q", type=float, default=0.05,
                        help="Germline length interval [0-1] for masking haplotagging errors (%(default)s)")
    parser.add_argument("-p", "--plump-threshold", type=int, default=50,
                        help="Plump type annotation placed on alleles with germline spread ≥ (%(default)s)")
    parser.add_argument("--no-mask", action='store_true',
                        help="Don't mask potential haplotagging errors (%(default)s)")
    parser.add_argument("--abs-delta", action="store_true",
                        help="Calculate abs(∆) (%(default)s)")
    parser.add_argument("--min-cov", type=int, default=30,
                        help="Minimum coverage over the locus (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (%(default)s)")

    args = parser.parse_args(args)
    return args



def flatness_metrics(lengths: np.ndarray) -> dict:
    """
    How flat is the distribution
    """
    lengths = np.asarray(lengths, dtype=float)
    lengths = lengths[~np.isnan(lengths)]

    rng = lengths.max() - lengths.min()
    p10, p50, p90 = np.percentile(lengths, [10, 50, 90])
    q1, q3 = np.percentile(lengths, [25, 75])

    return {
        "range": rng,
        "excess_kurtosis": stats.kurtosis(lengths, fisher=True, bias=False),
        "skewness": stats.skew(lengths, bias=False),
        #"skewtest": stats.skewtest(lengths).pvalue if len(lengths) >= 30 else np.nan,
        # core-spread / total-range: near 1 = flat, near 0 = peaked w/ outlier tails
        "p10_90_over_range": (p90 - p10) / rng if rng > 0 else np.nan,
        "iqr_over_range": (q3 - q1) / rng if rng > 0 else np.nan,
        # normalized Shannon entropy of the histogram: 1 = perfectly uniform
        "norm_entropy": _normalized_entropy(lengths),
        #"cv": lengths.std(ddof=1) / lengths.mean() if lengths.mean() != 0 else np.nan,
    }

def _normalized_entropy(lengths: np.ndarray, bins: int = 50) -> float:
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

    if len(reads) < args.min_cov:
        return None, None, None

    reads['donor'] = reads['smhtid'].apply(lambda x: x.donor)
    reads['protocol'] = reads['smhtid'].apply(lambda x: x.protocol)
    try:
        reads, germ = perform_clustering(reads,
                                        min_bandwidth=args.min_bandwidth,
                                        germ_vaf=args.germ_vaf,
                                        germ_q=args.germ_q,
                                        absolute=args.abs_delta,
                                        fix_haps=not args.no_mask,
                                        logging=False)
    except Exception as e:
        print(f"Exception on {chrom}:{start}-{end} {e}", file=sys.stderr)
        return None, None, None

    reads.insert(0, 'end', locus[2])
    reads.insert(0, 'start', locus[1])
    reads.insert(0, 'chrom', locus[0])

    germ.insert(0, 'end', locus[2])
    germ.insert(0, 'start', locus[1])
    germ.insert(0, 'chrom', locus[0])
    donor = all_smhtids[0].donor

    # for each hap, calculate plump metrics
    #for i in germ[(germ['upper'] - germ['lower']) >= args.plump_threshold].index.levels[1].unique():
    for i in germ.index.levels[1].unique():
        # only check the germline reads
        sub = reads[(reads['hap'] == i) & reads['is_germ']]
        metrics = flatness_metrics(sub['length'].values)
        germ.loc[(donor, i), list(metrics.keys())] = list(metrics.values())
    
    # state of 1 if there's a plume
    germ['state'] = ((germ['upper'] - germ['lower']) >= args.plump_threshold).astype(int)

    filtered = None
    if (~reads['is_germ']).any():
        try:
            filtered = coverage_filter(reads, args.min_reads_donor, args.min_reads_tissue)
        except Exception as e:
            print(f"Exception on {chrom}:{start}-{end} {e}", file=sys.stderr)
            return None, None, None

        if not filtered.empty:
            filtered.insert(0, 'end', locus[2])
            filtered.insert(0, 'start', locus[1])
            filtered.insert(0, 'chrom', locus[0])
            # Found a plume as well
            germ.loc[germ.index.levels[1].isin(filtered['hap'].unique()), 'state'] += 2

    return reads, germ, filtered

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
    reads_out = open(f'{args.output}.anno_reads.tsv', 'w')
    unst_out = open(f'{args.output}.unstable.tsv', 'w')
    f_read = True # For header writing, reads and unst are done together
    # Debug
    #for locus, h1_parts, h2_parts in stream_qdpi(inputs, to_analyze):
        #data, germ, filtered = run_analysis(locus, h1_parts, h2_parts, smhtids, args)
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(run_analysis, locus, h1_parts, h2_parts, smhtids, args)
            for locus, h1_parts, h2_parts in stream_qdpi(inputs, to_analyze)
        ]

        for future in tqdm(as_completed(futures), total=len(futures)):
            data, germ, filtered = future.result()

            if data is None:
                continue
            
            germ.to_csv(germ_out, sep='\t', header=f_germ)
            f_germ = False

            if filtered is not None:
                data.drop(columns=['donor', 'protocol'], inplace=True)
                data.to_csv(reads_out, sep='\t', index=False, header=f_read, float_format="%.1f")
                filtered.to_csv(unst_out, sep='\t', index=False, header=f_read)
                f_read = False

if __name__ == '__main__':
    scan_main(sys.argv[1:])
