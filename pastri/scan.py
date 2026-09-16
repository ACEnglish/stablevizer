import os
import sys
import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed

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
    parser.add_argument("--no-mask", action='store_true',
                        help="Don't mask potential haplotagging errors (%(default)s)")
    parser.add_argument("--min-cov", type=int, default=30,
                        help="Minimum coverage over the locus (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of threads (%(default)s)")

    args = parser.parse_args(args)
    return args


def run_analysis(locus, h1_parts, h2_parts, smhtids, min_cov):
    _, start, end = locus
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
    data = pd.DataFrame({'smhtid': all_smhtids,
                         'hap': all_haps,
                         'length': all_lengths
    })

    if len(data) < min_cov:
        return locus, None, None

    data['donor'] = data['smhtid'].apply(lambda x: x.donor)
    data['protocol'] = data['smhtid'].apply(lambda x: x.protocol)
    data, germ = perform_clustering(data,
                                    min_bandwidth=10,
                                    germ_vaf=0.80,
                                    germ_q=1,
                                    absolute=False,
                                    fix_haps=True,
                                    logging=False)
    data['chrom'] = locus[0]
    data['start'] = locus[1]
    data['end'] = locus[2]

    germ['chrom'] = locus[0]
    germ['start'] = locus[1]
    germ['end'] = locus[2]

    return locus, data, germ

def scan_main(args):
    args = parse_args(args)

    to_analyze = None
    if args.subset:
        to_analyze = set()
        with open(args.subset, 'r') as fh:
            for line in fh:
                to_analyze.add(tuple(line.strip().split('\t')))

    inputs = open(args.in_qdpi, 'r').read().strip().split('\n')
    smhtids = [SMaHTid.from_re(os.path.basename(_)) for _ in inputs]

    germ_out = open(f'{args.output}.germline.tsv', 'w')
    f_germ = True # For header writing
    reads_out = open(f'{args.output}.anno_reads.tsv', 'w')
    unst_out = open(f'{args.output}.unstable.tsv', 'w')
    f_read = True # For header writing, reads and unst are done together
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(run_analysis, locus, h1_parts, h2_parts, smhtids, args.min_cov)
            for locus, h1_parts, h2_parts in stream_qdpi(inputs, to_analyze)
        ]

        for future in tqdm(as_completed(futures), total=len(futures)):
            locus, data, germ = future.result()
            
            # Write
            if data is None:
                continue

            germ.to_csv(germ_out, sep='\t', header=f_germ)
            f_germ = False

            if (~data['is_germ']).any():
                filtered = coverage_filter(data, 3, 3) 
                if not filtered.empty:
                    data.drop(columns=['donor', 'protocol'], inplace=True)
                    data.to_csv(reads_out, sep='\t', index=False, header=f_read, float_format="%.1f")
                    filtered['chrom'] = locus[0]
                    filtered['start' ] = locus[1]
                    filtered['end'] = locus[2]
                    filtered.to_csv(unst_out, sep='\t', index=False, header=f_read)
                    f_read = False

if __name__ == '__main__':
    scan_main(sys.argv[1:])
