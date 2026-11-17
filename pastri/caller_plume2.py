import os
import re
import sys
import pandas as pd

results = []
for bdir in sys.argv[1:]:
    motif_fn = os.path.join(bdir, 'motif_annotations.tsv')
    unstable_fn = os.path.join(bdir, 'instability_table.tsv')
    if not os.path.exists(motif_fn) or not os.path.exists(unstable_fn):
        continue
    motifs = pd.read_csv(motif_fn, sep='\t')
    unstable = pd.read_csv(unstable_fn, sep='\t')

    motifs.set_index(['locus', 'region'], inplace=True)
    keep_unstable = unstable[unstable['instability'] >= 0.02].set_index(['locus', 'region'])
    sub_motifs = motifs.loc[keep_unstable.index]

    def summarize(state, values):
        # Quick length summary
        def subdo(arr):
            return [arr.mean().round(1), arr.median().round(1), arr.max().round(1), arr.min().round(1)]
        return subdo(values[state]) + subdo(values[~state]) + [(~state).sum()]

    rows = []
    for (l,r), sub in sub_motifs.groupby(['locus', 'region']):
        m_pur = sub.groupby(['is_somatic'])['purity'].mean().reindex([True, False]).round(4)
        rows.append([l, r]
                    + [",".join(sub['motif'].unique()), len(sub), m_pur.loc[False], m_pur.loc[True]]
                    + summarize(~sub['is_somatic'], sub['delta_units'])
                   )
    motif_result = pd.DataFrame(rows, columns=['locus', 'region', 'motif', 'coverage', 'germ_purity', 'soma_purity',
                                         'germ_mean', 'germ_median', 'germ_max', 'germ_min',
                                         'soma_mean', 'soma_median', 'soma_max', 'soma_min',
                                        'soma_reads',]).set_index(['locus', 'region'])
    motif_result['instability'] = keep_unstable['instability'].round(4)
    motif_result.reset_index(inplace=True)
    motif_result[['donor', 'protocol', 'hap']] = motif_result['locus'].str.split('_', expand=True)
    motif_result['vaf'] = (motif_result['soma_reads'] / motif_result['coverage']).round(4)


    # I could just do this at the end, but whatever
    reorder_cols = ['locus', 'donor', 'protocol', 'hap', 'region', 'soma_reads', 
                    'coverage', 'vaf', 'instability', 'germ_purity', 'soma_purity']
    for col_name in reorder_cols[::-1]:
        col = motif_result.pop(col_name)
        motif_result.insert(0, col_name, col)

    motif_result = motif_result.drop(columns=['locus'])
    chrom, start, end = re.split(':|-', bdir.split('/')[-2])
    motif_result.insert(0, 'chrom', chrom)
    motif_result.insert(1, 'start', start)
    motif_result.insert(2, 'end', end)
    results.append(motif_result)

results = pd.concat(results)
results.to_csv("delme.tsv", sep='\t', index=False)


