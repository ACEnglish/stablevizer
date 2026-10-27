"""
load reads.fasta
load unstable.tsv
load germline.tsv
for each donor/tissue/haplotype in unstable.tsv, 
    Fetch its reads,
         assign to haplotype based on its read-length overlap to the germline (both means keep unaltered)
             And then msa / annotate.
                 if the annotation returns... well, anything, we can keep the unstable.
                     if the annotation returns nothing, then we can drop the unstable.
"""
import os
import sys
import argparse 
from collections import defaultdict

import pysam
import pyabpoa
import pandas as pd
from pastri import SMaHTid
from pastri.annotate_tr_seqs import annotate_locus
from tqdm.contrib.concurrent import process_map

def parse_args(args):
    """
    Argument parser
    """
    parser = argparse.ArgumentParser(prog="motif_caller", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--fasta", required=True,
                        help="Reads.fasta")
    parser.add_argument("-u", "--unstable", required=True,
                        help="Unstable.tsv")
    parser.add_argument("-g", "--germline", required=True,
                        help="Germline.tsv")
    parser.add_argument("-o", "--output", default="read_motifs.tsv",
                        help="Read output (%(default)s)")
    parser.add_argument("-c", "--calls", default="motif_calls.tsv",
                        help="Donor-tissue-haplotype calls (%(default)s)")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Threads (%(default)s)")
    args = parser.parse_args(args)
    if not os.path.exists(args.fasta) or not os.path.exists(args.unstable) or not os.path.exists(args.germline):
        print("Missing File(s)")
        sys.exit(1)
    return args

def make_msa(reads):
    """
    Returns tuple of all the names and all their sequences
    """
    all_seqs = sorted(reads, reverse=True, key=lambda x: len(x[1]))
    all_names, all_seqs = map(list, zip(*all_seqs))
    aligner = pyabpoa.msa_aligner()
    aln_result = aligner.msa(all_seqs, False, True)
    return all_names, aln_result.msa_seq

def core_process(item):
    """
    """
    donor_tiss, reads = item
    names, seqs = make_msa(reads)
    # trf, ploidy, min_allele_frac
    return annotate_locus(names, seqs, donor_tiss)

def motif_caller_main(args):
    args = parse_args(args)
    
    fasta = pysam.FastaFile(args.fasta)
    read_table = pd.Series(fasta.references, name='read_name').to_frame()
    read_table['smahtid'] = read_table['read_name'].apply(SMaHTid.from_re)
    read_table['donor'] = read_table['smahtid'].apply(lambda x: x.donor)
    read_table['protocol'] = read_table['smahtid'].apply(lambda x: x.protocol)
    read_table['orig_hap'] = read_table['read_name'].apply(lambda x: int(x.split('_')[2]))
    
    read_table.set_index(['donor', 'protocol'], inplace=True)
    unstable = pd.read_csv(args.unstable, sep='\t').set_index(['donor', 'protocol'])

    read_table = read_table.loc[unstable.index].reset_index(level=1)
    germline = pd.read_csv(args.germline, sep='\t').set_index(['donor', 'hap'])
    
    m_reads = defaultdict(list)

    for _, row in read_table.reset_index().iterrows():
        read = fasta.fetch(row['read_name'])

        # If this is part of an unstable donor tissue
        unstable_donor_tissue = row['orig_hap'] in unstable.loc[(row['donor'], row['protocol'])][['hap']].values
        if not unstable_donor_tissue:
            continue

        # and it isn't a mask candidate
        read_len = len(read)
        m_germ = germline.loc[(row['donor'], row['orig_hap'])]
        is_m_germ = m_germ['lower'] <= read_len <= m_germ['upper']

        other_hap_key = (row['donor'], row['orig_hap'] % 2 + 1)
        if other_hap_key in germline.index:
            o_germ = germline.loc[other_hap_key]
            is_o_germ = o_germ['lower'] <= read_len <= o_germ['upper']
        else:
            is_o_germ = False
        # somatic but was flipped
        if not is_m_germ and is_o_germ:
            continue

        # use the read
        m_key = "%s_%s_%d" % (row['donor'], row['protocol'], row['orig_hap'])
        m_reads[m_key].append((row['read_name'], read))

    rows = []
    for result in process_map(core_process, m_reads.items(), chunksize=min(1, max_workers=args.threads):
        rows += result
    
    if len(rows) == 0:
        sys.exit(0)

    result = pd.DataFrame(rows)
    result.to_csv(args.output, sep='\t', index=False)
    # Make calls
    calls = result[result['status'] != 'germline']
    if calls.empty:
        with open(args.calls, 'w') as fout:
            fout.write("donor\tprotocol\thaplotype\tcount\n")
        sys.exit(0)
    # This can be better.
    # I want donor tissue haplotype Motif(N±M) Support Coverage VAF
    # Where N±M is the Mean and STD of the copy number? and they're comma separated? Meh, hard to parse
    # Just be able to pull the reads from `--output`
    calls = calls[['locus', 'read']].sort_values(by=['locus', 'read']).drop_duplicates()['locus'].value_counts()
    calls = calls.to_frame().reset_index()
    calls[['donor', 'protocol', 'haplotype']] = calls['locus'].str.split('_', expand=True)
    calls = calls.drop(columns=['locus'])[['donor', 'protocol', 'haplotype', 'count']]
    calls.to_csv(args.calls, sep='\t', index=False)

if __name__ == '__main__':
    motif_caller_main(sys.argv[1:])
