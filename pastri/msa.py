"""
Utility for making MSA from a locus' reads.fasta with helpers to subset by donor/tissue_abv/haplotype
Outputs to f'{args.output_prefix}_{donor}_{tissue_abv}_{haplotype}.msa'

Allow `--cross-tissue` and/or `--cross-haplotype`

And allow `--consensus` mode to write output.consensus.fa

And encourage alv for visualizing (I might want to hook over that so that the less -SR or whatever is possible
"""
import os
import sys
import argparse

import pysam
import pyabpoa
import pandas as pd

from tqdm import tqdm
from pastri import SMaHTid

def parse_args(args):
    parser = argparse.ArgumentParser(prog="msa", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-f", "--fasta", default=None,
                        help="Fasta file of reads")
    parser.add_argument("--save-fasta", action="store_true",
                        help="Save a plain fasta file")
    parser.add_argument("-o", "--output", required=True,
                        help="Output prefix")
    filtsg = parser.add_argument_group("Fasta Subsetting Arguments")
    filtsg.add_argument("-d", "--donor", default=None, type=str,
                        help="Subset to donor(s). comma-separated")
    filtsg.add_argument("-p", "--protocol", default=None, type=str,
                        help="Subset to protocol(s). comma-separated")
    filtsg.add_argument("-t", "--tissue", default=None, type=str,
                        help="Subset to tissue_abv(s). comma-separated")
    filtsg.add_argument("-H", "--haplotype", default=None, type=str,
                        help="Subset to haplotype(s). comma-separated")
    filtsg.add_argument("-P", "--platform", default=None, type=str,
                        help="Subset to platform_name(s). comma-separated")
    filtsg.add_argument("--allow-unphased", action="store_true",
                        help="Allow unphased reads to be MSA'd")

    consog = parser.add_argument_group("Fasta Consolidation Arguments")
    consog.add_argument("-C", "--cross-tissue", action="store_true",
                        help="Don't separate by tissue")
    consog.add_argument("-c", "--cross-haplotype", action="store_true",
                        help="Don't separate by haplotype")

    args = parser.parse_args(args)
    
    # Split comma-separated fields
    args.donor = args.donor if not args.donor else args.donor.split(',')
    args.protocol = args.protocol if not args.protocol else args.protocol.split(',')
    args.tissue = args.tissue if not args.tissue else args.tissue.split(',')
    args.haplotype = args.haplotype if not args.haplotype else args.haplotype.split(',')
    args.platform = args.platform if not args.platform else args.platform.split(',')

    return args

def subset_exclude(mid, hap, args):
    """
    For each of the subsetting flags, exclude if it is not None and value is in the list
    """
    if args.donor and not mid.donor in args.donor:
        return True
    if args.protocol and not mid.protocol in args.protocol:
        return True
    if args.tissue and not mid.tissue_abv in args.tissue:
        return True
    if args.haplotype and not hap in args.haplotype:
        return True
    if args.platform and not mid.platform_name in args.platform:
        return True
    return False

def msa_main(args):
    args = parse_args(args)
    
    fasta = pysam.FastaFile(args.fasta)
    table = []
    for i in fasta.references:
        sample, region, hap, rname = i.split('_', maxsplit=3)
        sample = SMaHTid.from_re(i)
        if subset_exclude(sample, hap, args):
            continue
        if hap == '0' and not args.allow_unphased:
            continue
        table.append([sample, region, int(hap), rname, i])

    table = pd.DataFrame(table, columns=['sample', 'region', 'hap', 'rname', 'entry'])
    table['donor'] = table['sample'].apply(lambda x: x.donor)
    table['protocol'] = table['sample'].apply(lambda x: x.protocol)

    # Setup the groupby
    m_group = ['donor']
    if not args.cross_tissue:
        m_group.append('protocol')
    if not args.cross_haplotype:
        m_group.append('hap')

    aligner = pyabpoa.msa_aligner()
    for values, sub in tqdm(table.groupby(m_group)):
        values = list(values)
        if args.cross_tissue:
            values.insert(1, 'all')
        if args.cross_haplotype:
            values.insert(2, 'all')

        all_names = []
        all_seqs = []
        for e_name in sub['entry']:
            all_names.append(e_name)
            all_seqs.append(fasta.fetch(e_name))

        all_seqs = sorted(all_seqs, reverse=True, key=lambda x: len(x))
        aln_result = aligner.msa(all_seqs, False, True)
        
        if args.save_fasta:
            fa_out = open(f'{args.output}_{values[0]}_{values[1]}_{values[2]}.fasta', 'w')

        with open(f'{args.output}_{values[0]}_{values[1]}_{values[2]}.msa', 'w') as fout:
            for name, msa in zip(all_names, aln_result.msa_seq):
                print(f'>{name}\n{msa}', file=fout)
                if args.save_fasta:
                    print(f">{name}\n{msa.replace('-','')}", file=fa_out)
        if args.save_fasta:
            fa_out.close()



