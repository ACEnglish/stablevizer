"""
Utility to setup a locus for analysis by extractng deltas/sequences from qdpi output
pastri setup chr8:36230992-36231078 -q qdpi_paths.txt -f fasta_paths.txt -o loci/
Do I want the CNV Filter available here or???
And do we put the CNV into pastri?
"""
import os
import sys
import argparse

import pysam
import tabix
from tqdm import tqdm
from pastri.smahtid import SMaHTid

def parse_args(args):
    parser = argparse.ArgumentParser(prog="tr", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--region", required=True,
                        help="Region to extract chrom:start-end")
    parser.add_argument("-q", "--qdpi", required=True,
                        help="File of list of qdpi files")
    parser.add_argument("-f", "--fasta", default=None,
                        help="File of list of fasta files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output prefix (makes a directory)")

    args = parser.parse_args(args)
    return args

def locus_setup_main(args):
    args = parse_args(args)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    out_reads = os.path.join(args.output, 'reads.tsv')
    # Build reads.tsv - this is just extractor.py
    print("Extracing qdpi deltas", file=sys.stderr)
    with open(args.qdpi, 'r') as files, open(out_reads, 'w') as fout:
        print("donor\tprotocol\tfull_name\tlength\thap", file=fout)
        for fn in tqdm(files.read().strip().split('\n')):
            mid = SMaHTid.from_re(os.path.basename(fn))
            fh = tabix.open(fn)
            try:
                result = fh.querys(args.region)
            except Exception:
                tqdm.write(f"ERROR: Unable to parse reads for {mid.full_name} on {fn}", file=sys.stderr)
                continue

            for (chrom, start, end, cov, h0, h1, h2) in result:
                span = int(end) - int(start)
                for hap_idx, h in enumerate((h0, h1, h2)):
                    if h == '.':
                        continue
                    for i in h.split(','):
                        print(mid.donor, mid.protocol, mid.full_name, span + int(i), hap_idx, sep='\t', file=fout)
    
    if not args.fasta:
        return

    # Build seqs.fa
    out_fasta = os.path.join(args.output, 'reads.fasta')
    args.region += '_' # Ensure we can find it
    print("Extracing fasta sequences", file=sys.stderr)
    with open(args.fasta, 'r') as files, open(out_fasta, 'w') as fout:
        for fn in tqdm(files.read().strip().split('\n')):
            mid = SMaHTid.from_re(os.path.basename(fn))
            fasta = pysam.FastaFile(fn)

            for i in fasta.references:
                if i.startswith(args.region):
                    fout.write(f">{mid.full_name}_{i}\n{fasta.fetch(i)}\n")
