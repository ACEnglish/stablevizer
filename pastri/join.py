"""
Glue script to put a TR pattern annotation (plump/plump/both) as well as cnv annotation
"""
import pandas as pd
import glob
from tqdm import tqdm
import numpy as np

for germ_fn in tqdm(glob.glob("*.germline.tsv")):
    donor = germ_fn.split('.')[0]
    tqdm.write(donor)
    cnv_fn = f"~/fritz/english/SMaHT_donors/TR_Work/tr_cnv/{donor}.regions.bed.gz"

    germ = pd.read_csv(germ_fn, sep='\t').set_index(['chrom', 'start', 'end'])
    cnv = pd.read_csv(cnv_fn, sep='\t', names=['chrom', 'start', 'end', 'cnv_p', 'cnv_n']).set_index(['chrom', 'start', 'end'])
    germ[['cnv_p', 'cnv_n']] = cnv

    # Also add in plume/plump annotations
    germ = germ.reset_index().set_index(['chrom', 'start', 'end', 'hap'])
    unst = pd.read_csv(f"{donor}.unstable.tsv", sep='\t').set_index(['chrom', 'start', 'end', 'hap'])
    is_plume = germ.index.isin(unst.index)
    is_plump = (germ['upper'] - germ['lower']) >= 50
    cond = [
        is_plume & is_plump,
        is_plume,
        is_plump
    ]
    choices = ['both', 'plume', 'plump']
    germ['state'] = np.select(cond, choices, default='neither')
    germ.to_csv(f"{donor}.germline_anno.tsv", sep='\t')
