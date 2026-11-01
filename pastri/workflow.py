import os
import sys
import glob
import json
from collections import defaultdict

import pysam
import numpy as np
import pandas as pd
from tqdm.contrib.concurrent import process_map

# Todo: Reorganize repository around these core functions
from pastri import SMaHTid
from pastri.scan import flatness_metrics
from pastri.motif_caller import make_msa
from pastri.plume import perform_clustering
from pastri.annotate_tr_seqs import annotate_locus
from pastri.experimental_caller import instability_metrics
from pastri.readcounts import ReadCounts
from pastri.origin import annotate_origin

# ToDo parameterize
# Maximum percent of reads from a donor that can be unphased before we are no longer interested in analyzing it
MAX_UNPHASED_PCT = 0.50
# Maximum allele balance allowed to process a donor
MAX_ALLELE_BALANCE = 0.60
# Minimum span of germline read lengths to assume the site has some overdispersion in germline signal
MIN_OVERDISPERSION_SPAN = 50
# Minimum normalized entropy of germline reads with high overdispersion span to consider the site as 'plump'
MIN_OVERDISPERSION_ENTROPY = 0.65
# Minimum number of outlier reads in a single tissue to consider a plume
MIN_OUTLIER_SUPPORT = 3
# Minimum number of basepairs away from germline to be considered somatic
MIN_LENGTH_CHANGE = 10
# Alpha for origin annotations
ORIGIN_ALPHA = 0.05
# Multiprocessing
THREADS = 8

def check_reference_filters(in_dir):
    """
    - Start with the ard and segdup and just don't go into those loci

    Returns a list of Flags this overlaps from ARD, SegDup, chrX
    """
    is_ard = os.path.exists(os.path.join(in_dir, "ard.tsv"))
    is_segdup = os.path.exists(os.path.join(in_dir, "segdup.tsv"))
    is_chrX = os.path.basename(in_dir).startswith('chrX')
    
    ret = []
    if is_ard:
        ret.append('ARD')
    if is_segdup:
        ret.append('SegDup')
    if is_chrX:
        ret.append('chrX')
    return ret


def make_read_table(in_dir):
    """
    - Then go straight to the the fasta and get their lengths so that you can make a table

    Returns a table with columns:
    - read: Full fasta entry name
    - length: length of the read
    - sample: SMaHTid sample name
    - locus: chrom:start-end locus identifier
    - hap: assigned haplotype
    - read_name: Sequencing read name
    """
    fasta_file = os.path.join(in_dir, "reads.fasta")
    fasta = pysam.FastaFile(fasta_file)
    table = pd.DataFrame(zip(fasta.references, fasta.lengths), columns=['read', 'length'])
    table[['sample', 'locus', 'hap', 'read_name']] = table['read'].str.split('_', n=3, expand=True)
    table['sample'] = table['sample'].apply(SMaHTid)
    table['donor'] = table['sample'].apply(lambda x: x.donor)
    table['protocol'] = table['sample'].apply(lambda x: x.protocol)
    table['hap'] = table['hap'].astype(int)
    return table, fasta

def read_depth_anomaly(in_dir):
    """
    - Check the depth anomaly and mask any of those donors completely from the table

    Updates the table in-place with CovAnom boolean where True means we don't like the coverage here
    Returns the table as well as a pd.Series of the violating donors
    """
    cov_anom_fn = os.path.join(in_dir, "depth_anomaly.txt")
    df = pd.read_csv(cov_anom_fn, sep='\t', names=['chrom', 'start', 'end', 'pval', 'null', 'donor'])
    return df

def measure_read_counts(table, max_up_pct=MAX_UNPHASED_PCT, max_ab=MAX_ALLELE_BALANCE):
    """
    Counts number of reads by haplotype per-donor
    Sets the allele balance and the unphased read percent
    """
    counts = (
        table.groupby(['donor', 'hap'])
             .size()
             .unstack(fill_value=0)
             .reindex(columns=[0, 1, 2], fill_value=0)
    )
    counts['total_coverage'] = counts.sum(axis=1)
    counts['unphased_pct'] = counts[0] / counts['total_coverage']
    counts['allele_balance'] = counts[[1, 2]].max(axis=1) / counts[[1, 2]].sum(axis=1)
    counts['highly_unphased'] = counts['unphased_pct'] > max_up_pct
    counts['poor_allele_balance'] = counts['allele_balance'] > max_ab
    return counts.reset_index()

def check_dispersion(reads):
    """
    Calculates the flatness metrics of germline read distributions and returns as a table
    """
    rows = []
    for (donor, hap), sub in reads[reads['is_germ']].groupby(['donor', 'hap']):
        # If you commit to dropping the other measurements, you can probably just call norm_entropy directly
        metrics = flatness_metrics(sub['length'].values)
        metrics['donor'] = donor
        metrics['hap'] = hap
        rows.append(metrics)
    plump_table = pd.DataFrame(rows)
    plump_table['plump_signal'] = (plump_table['range'] >= MIN_OVERDISPERSION_SPAN) & (plump_table['norm_entropy'] >= MIN_OVERDISPERSION_ENTROPY)
    return plump_table.set_index(['donor', 'hap'])

def check_outliers(reads, min_reads=MIN_OUTLIER_SUPPORT):
    """
    Ensure that there are at least min_reads in a single protocol
    Need up update a Series with index of donor hap and boolean plume_signal
    And a table of the donor-tissue-haplotypes that pass the filter.
    """
    counts = reads.groupby(['donor', 'protocol', 'hap', 'is_germ']).size().unstack()
    return counts[counts[False] >= min_reads]

def calculate_origin(reads):
    """
    Hooks into origin annotation.
    returns a dict with origin_table and global_origin keys
    """
    ret = {}
    origin_parts = []
    read_parts = []
    for (donor, hap), sub_reads in reads.groupby(["donor", "hap"]):
        donor_key = f'{donor}.{hap}'
        all_cols = sub_reads['sample'].unique()
        ref = (sub_reads[~sub_reads['is_somatic']]
                          .groupby(['sample']).size()
                          .reindex(all_cols, fill_value=0))
        ref.name = donor_key
        ref = ref.to_frame().T

        alt = (sub_reads[sub_reads['is_somatic']]
                          .groupby(['sample']).size()
                          .reindex(all_cols, fill_value=0))

        alt.name = donor_key
        alt = alt.to_frame().T
        read_counts = ReadCounts(ref.values, alt.values,
                                 index=ref.index,
                                 columns=ref.columns)

        # If there's nothing somatic, there's nothing to look at
        if alt.iloc[0].sum() != 0:
            try:
                origin = annotate_origin(read_counts, ORIGIN_ALPHA, MIN_OUTLIER_SUPPORT)
                origin_parts.append(origin)
            except Exception:
                # I got edge cases around when e.g. only Blood/Fibro have read support
                # Which is weird, but need to just skip them
                pass
        read_parts.append(read_counts)

    if origin_parts:
        ret['origin_table'] = pd.concat(origin_parts)

    # Safety
    if not read_parts:
        return ret

    # Experimental Global origin annotation across consolidated reads
    # Just to see if we can find an easy classification of the region
    all_protocol = reads['sample'].apply(lambda x: x.protocol).unique()
    all_plat = reads['sample'].apply(lambda x: x.platform).unique()
    new_donor = []
    for i in all_protocol:
        for j in all_plat:
            new_donor.append(SMaHTid(f"SMHTALL-{i}-XX-XX-{j}001-all"))

    ref = np.zeros(len(new_donor))
    alt = np.zeros(len(new_donor))

    dest_map = {f"{name.protocol}-{name.platform}": i for i, name in enumerate(new_donor)}

    for sample in read_parts:
        idx = np.array([dest_map[f"{name.protocol}-{name.platform}"] for name in sample.col])
        np.add.at(ref, idx, sample.ref[0])
        np.add.at(alt, idx, sample.alt[0])

    read_counts = ReadCounts([ref], [alt], index=['SMHTALL'], columns=new_donor)
    ret['global_origin'] = annotate_origin(read_counts)
    return ret

def main(in_dir):
    locus = os.path.basename(in_dir)
    o_dir = os.path.join(in_dir, "caller_results")
    if not os.path.exists(o_dir):
        os.mkdir(o_dir)
    meta = {'locus': locus}
    meta['ref_filters'] = check_reference_filters(in_dir)
    # We don't process reference filtered loci anymore
    if meta['ref_filters']:
        with open(os.path.join(o_dir, 'meta.json'), 'w') as fout:
                json.dump(meta, fout)
        return meta

    reads, fasta = make_read_table(in_dir)

    # Drop coverage anomaly donors
    cov_anom = read_depth_anomaly(in_dir)
    meta['covanom_donors'] = list(cov_anom['donor'].unique())
    reads.drop(reads.index[reads['donor'].isin(cov_anom['donor'])], inplace=True)

    # Drop donors that are highly unphased 
    read_counts = measure_read_counts(reads)
    bad_donors = read_counts[read_counts['highly_unphased']]['donor'].unique()
    
    meta['preclustering_read_counts'] = json.loads(read_counts.set_index('donor').to_json())
    meta['poorphase_donors'] = list(bad_donors)
    reads.drop(reads.index[reads['donor'].isin(bad_donors)], inplace=True)
        

    # Need to drop unphased here before clustering
    reads.drop(reads.index[reads['hap'] == 0], inplace=True)

    # Perform clustering
    # - Group by donor tissue haplotype and run the plume annotation
    # This should produce an updated set of haplotypes and is_somatic information per-read
    reads, germ = perform_clustering(reads, germ_q=1, logging=False)

    # - Check the overall allele balance (post re-haplotyping) and drop any donors which don't seem normal
    read_counts = measure_read_counts(reads)
    bad_donors = read_counts[read_counts['poor_allele_balance']]['donor'].unique()
    meta['poor_allelebalance_donors'] = list(bad_donors)
    
    reads.drop(reads.index[reads['donor'].isin(bad_donors)], inplace=True)
    # Don't carry around the filtered donors
    germ.drop(index=bad_donors, level='donor', inplace=True)
    
    plump_table = check_dispersion(reads)
    germ['plump_signal'] = plump_table["plump_signal"]
    
    plume_table = check_outliers(reads)
    germ['plume_signal'] = germ.index.isin(plume_table.reset_index(level=1).index)

    # Expand the read annotations with their motif
    read_sequences = defaultdict(list)
    reads.set_index(['donor', 'protocol', 'hap'], inplace=True)
    
    for key, read_row in reads.loc[plume_table.index].iterrows():
        key = "%s_%s_%d" % key
        read_name = read_row['read']
        read_sequences[key].append((read_name, fasta.fetch(read_name)))
    
    motif_annotations = []
    for donor_tiss, sequences in read_sequences.items():
        names, seqs = make_msa(sequences)
        motif_annotations += annotate_locus(names, seqs, donor_tiss, min_allele_frac=0.01)
        
    motif_annotations = pd.DataFrame(motif_annotations)

    ret = {"reads": reads,
            "germ_table": germ, 
            "plume_table": plume_table, 
            "plump_table": plump_table, 
            "motif_annotations": motif_annotations,
            "meta": meta}

    # We need this guard, but I don't like what this is doing to the code structure
    # But this section is all about pullin over motif annotated and is somatic information
    # Plus some summary information
    if not motif_annotations.empty:
        motif_annotations.drop(columns=['status'], inplace=True)
        motif_annotations['is_somatic'] = (motif_annotations['tr_len_bp'] - motif_annotations['germline_alleles_bp']).abs() >= MIN_LENGTH_CHANGE
        
        # Mark the reads that we got motif annotations on
        agg = (
            motif_annotations
                .groupby('read')[['region', 'motif', 'delta_units']]
                .agg(lambda s: ','.join(map(str, s)))
        )
        reads = reads.reset_index().join(agg, on='read')
        
        # And if they are somatic
        reads = reads.join(
            motif_annotations.groupby('read')['is_somatic'].agg('any'),
            on='read')
        reads['is_somatic'] = reads['is_somatic'].fillna(False).astype(bool)

        # Record the old is_germ just so we have it
        motif_annotations['is_germ'] = motif_annotations['read'].map(reads.set_index(['read'])['is_germ'])

        motif_table = reads[~reads['delta_units'].isna()].groupby(['donor', 'protocol', 'hap', 'is_somatic']).size().unstack()
        motif_table.columns = motif_table.columns.map({False: 'germ', True: 'soma'})
        if len(motif_table.columns) == 2:
            motif_table['vaf'] = motif_table['soma'] / motif_table[['germ', 'soma']].sum(axis=1)

        instability_table = instability_metrics(motif_annotations, key_on='is_somatic')
        instability_table[['donor', 'protocol', 'haplotype']] = instability_table['locus'].str.split('_', expand=True)
        instability_table.drop(columns=['donor'], inplace=True)
        
        ret['reads'] = reads
        ret['meta']['motif_failures'] = list(set(plume_table.index) - set(motif_table.index))
        ret['meta']['unstable_protocols'] = json.loads(motif_table.reset_index()['protocol'].value_counts().to_json())

        ret["instability_table"] = instability_table
        ret["motif_table"] = motif_table
        
        # Origin
        ret.update(calculate_origin(reads))
        if 'global_origin' in ret:
            go = ret['global_origin']
            go = go[(go['direction'] == '+') | go['direction'].isna()]
            ret['meta']['origin'] = {
                    'protocols':list(go['col'].unique()),
                    'annotation':go.iloc[0]['annotation']
            }

    # I guess just save this in a subdir named caller_results/
    
    for k, d in ret.items():
        if k != 'meta':
            d.to_csv(os.path.join(o_dir, k + '.tsv'), sep='\t')
        else:
            with open(os.path.join(o_dir, 'meta.json'), 'w') as fout:
                json.dump(d, fout)
    return ret

if __name__ == '__main__':
    #b_dir = '/Users/English/code/SMaHT_TR_explore/analysis/debug_oct22_catalog'
    b_dir = sys.argv[1]
    results = process_map(main, glob.glob(os.path.join(b_dir, '*')), max_workers=THREADS, chunksize=1)
    #result = main("/Users/English/code/SMaHT_TR_explore/analysis/debug_oct22_catalog/chr2/chr2:69096232-69096559/")

