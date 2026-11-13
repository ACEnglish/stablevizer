"""
"""
import os
import sys
import json
import numpy as np
import pandas as pd

def main():
    bdir = sys.argv[1]
    meta = json.load(open(os.path.join(bdir, "meta.json")))
    read_fn = os.path.join(bdir, "reads.tsv")
    unst_fn = os.path.join(bdir, "instability_table.tsv")
    if not os.path.exists(read_fn) or not os.path.exists(unst_fn):
        sys.exit(0)

    reads = pd.read_csv(read_fn, sep='\t')
    unstable = pd.read_csv(unst_fn, sep='\t')

    # Subtract low-quality donors
    remove = meta['covanom_donors'] + meta['poorphase_donors'] + meta['poor_allelebalance_donors']
    sub_reads = reads[~reads['donor'].isin(remove)]

    # Subtract low-instability donors
    keep_unstable = unstable[unstable['instability'] >= 0.02]
    keys = sub_reads['donor'] + '_' + sub_reads['protocol'] + '_' + sub_reads['hap'].astype(str) 
    sub_reads = sub_reads[keys.isin(keep_unstable['locus'])]

    def summarize(state, values):
        # Quick length summary
        def subdo(arr):
            return [arr.mean().round(1), arr.median().round(1), arr.max().round(1), arr.min().round(1)]
        return subdo(values[state]) + subdo(values[~state]) + [(~state).sum()]

    def reformat(val):
        if val == 'NaN':
            return np.nan
        if ',' not in val:
            return float(val)
        vals = np.fromiter(val.split(','), float)
        return vals[np.argmax(np.abs(vals))]
        
    if sub_reads['delta_units'].dtype == 'str':
        sub_reads['delta_units'] = sub_reads['delta_units'].apply(reformat)

    out_rows = []
    for (d, t, h), df in sub_reads.groupby(['donor', 'protocol', 'hap']):
        out_rows.append([d, t, h] 
                        + summarize(df['is_germ'], df['delta'])
                        + [','.join(df['motif'].unique())]
                        + summarize(~df['is_somatic'], df['delta_units'])
                        + [len(df)]
                       )

    out = pd.DataFrame(out_rows, columns=['donor', 'protocol', 'haplotype',
                                          'is_germ_mean', 'is_germ_median', 'is_germ_max', 'is_germ_min',
                                          'not_germ_mean', 'not_germ_median', 'not_germ_max', 'not_germ_min',
                                          'not_germ_count', 
                                          'motifs',
                                          'not_soma_mean', 'not_soma_median', 'not_soma_max', 'not_soma_min',
                                          'is_soma_mean', 'is_soma_median', 'is_soma_max', 'is_soma_min',
                                          'is_soma_count',
                                          'coverage'
                                         ])

    unstable['donor'] = unstable['locus'].apply(lambda x: x.split('_')[0])
    m_score = unstable.groupby(['donor', 'protocol', 'haplotype'])['instability'].max().round(5)
    out.set_index(['donor', 'protocol', 'haplotype'], inplace=True)
    out['instability'] = m_score
    out.to_csv(os.path.join(bdir, "plume_caller_results.tsv"), sep='\t')

if __name__ == '__main__':
    try:
        main()
    except Exception as e:
        print(e)
        print(sys.argv[1])
        pass
