"""
For working on a single locus' stablevizer output
"""
import sys
import pandas as pd

from smahtkit import ReadCounts, SMaHTid
from smahtkit.origin import annotate_origin

MINSUP = 3
reads = pd.read_csv(sys.argv[1], sep='\t')

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
    # MASKING - Risky!
    alt[alt < MINSUP] = 0
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
output.to_csv("origin.tsv", index=False, sep='\t')

output = pd.concat(read_parts)
output.to_csv("square_reads.tsv", sep='\t', index=False)
