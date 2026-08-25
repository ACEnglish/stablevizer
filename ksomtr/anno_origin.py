"""
Analyze the anno_reads.tsv from stablevizer using smahtkit origin
"""
import sys
import argparse
import pandas as pd

from smahtkit import ReadCounts, SMaHTid
from smahtkit.origin import annotate_origin

parser = argparse.ArgumentParser(prog="anno_origin", description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("in_tsv", type=str,
                    help="Input reads.tsv")
parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                    help="Output origin.tsv annotations (%(default)s)")
parser.add_argument("-s", "--square-reads", default=None,
                    help="Save the per-sample ref/alt read counts to file (off)")
parser.add_argument("-m", "--min-sup", type=int, default=3,
                    help="Min number of reads in a sample to allow presence (%(default)s)")
args = parser.parse_args()
reads = pd.read_csv(args.in_tsv, sep='\t')

read_parts = []
parts = []
for (donor, hap), sub_reads in reads.groupby(["donor", "hap"]):
    donor_key = f'{donor}.{hap}'
    all_cols = sub_reads['full_name'].unique()

    ref = (sub_reads[sub_reads['is_germ']]
           .groupby(['full_name']).size()
           .reindex(all_cols, fill_value=0))
    ref.name = donor_key
    ref = ref.to_frame().T

    alt = (sub_reads[~sub_reads['is_germ']]
           .groupby(['full_name']).size()
           .reindex(all_cols, fill_value=0))
    # MASKING - fewer than min_sup isn't considered present
    alt[alt < args.min_sup] = 0
    alt.name = donor_key
    alt = alt.to_frame().T
    read_counts = ReadCounts(ref.values, alt.values,
                             index=ref.index,
                             columns=[SMaHTid(_) for _ in ref.columns])
    
    if alt.iloc[0].sum() != 0:
        origin = annotate_origin(read_counts)
        parts.append(origin)
    
    ref = ref.T[donor_key]
    ref.name = 'ref'

    alt = alt.T[donor_key]
    alt.name = 'alt'
    read_parts.append(pd.concat([ref, alt], axis=1).reset_index())

    
output = pd.concat(parts)
output.to_csv(args.output, index=False, sep='\t')

if args.square_reads is not None:
    output = pd.concat(read_parts)
    output.to_csv(args.square_reads, sep='\t', index=False)
