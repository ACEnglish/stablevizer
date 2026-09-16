"""
Utility to setup a locus for analysis by extractng deltas/sequences from qdpi output
pastri setup chr8:36230992-36231078 -q qdpi_paths.txt -f fasta_paths.txt -o loci/
Do I want the CNV Filter available here or???
And do we put the CNV into pastri?
"""
import os
import re
import sys
import argparse

import pysam
import tabix
from tqdm import tqdm
from pastri.smahtid import SMaHTid

def parse_args(args):
    parser = argparse.ArgumentParser(prog="tr", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-r", "--region", default=None,
                        help="Region to extract chrom:start-end")
    parser.add_argument("-R", "--regions", default=None,
                        help="Regions to annotate (bed file)")
    parser.add_argument("-q", "--qdpi", required=True,
                        help="File of list of qdpi files")
    parser.add_argument("-f", "--fasta", default=None,
                        help="File of list of fasta files")
    parser.add_argument("-o", "--output", required=True,
                        help="Output prefix (makes a directory with sub-dirs per-chrom)")

    args = parser.parse_args(args)
    if not (args.region is not None) ^ (args.regions is not None):
        print("Only set one of --region or --Regions", file=sys.stderr)
        sys.exit(1)
    return args

def locus_setup_main(args):
    args = parse_args(args)

    if not os.path.exists(args.output):
        os.mkdir(args.output)

    regions = []
    if args.region:
        chrom, start, end = re.split(':|-', args.region)
        start = int(start)
        end = int(end)
        regions = [[chrom, start, end]]
    else:
        with open(args.regions, 'r') as fh:
            for line in fh:
                d = line.strip().split('\t')
                d[1] = int(d[1])
                d[2] = int(d[2])
                regions.append(d[:3])

    # for each region, setup the sub-directory
    print(f"Making {len(regions)} directories", file=sys.stderr)
    delta_fns = []
    fasta_fns = []
    for chrom, start, end in regions:
        bdir = os.path.join(args.output, chrom)
        if not os.path.exists(bdir):
            os.mkdir(bdir)
        dest = os.path.join(bdir, f'{chrom}:{start}-{end}')
        #os.mkdir(dest)

        delta_fns.append(os.path.join(dest, 'reads.tsv'))
        fasta_fns.append(os.path.join(dest, 'reads.fasta'))
    
    print("Extracing qdpi deltas", file=sys.stderr)
    """
    # For each region
    for (reg_chrom, reg_start, reg_end), de_fn in tqdm(
            zip(regions, delta_fns), 
            total=len(regions), 
            desc="loci", 
    ):
        fout = open(de_fn, 'w')
        print("donor\tprotocol\tfull_name\tlength\thap", file=fout)
        with open(args.qdpi, 'r') as files:
            for qdpi_fn in tqdm(files.read().strip().split('\n'), desc="qdpi", leave=False):
                tbx_fh = tabix.open(qdpi_fn)
                mid = SMaHTid.from_re(os.path.basename(qdpi_fn))

                try:
                    result = tbx_fh.query(reg_chrom, reg_start, reg_end)
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
        fout.close()
    """
        

    # File opening / closing, but on smaller files.
    with open(args.qdpi, 'r') as files:
        # Open qdpi files once
        first = True
        for qdpi_fn in tqdm(files.read().strip().split('\n'), desc="qdpi"):
            tbx_fh = tabix.open(qdpi_fn)
            mid = SMaHTid.from_re(os.path.basename(qdpi_fn))

            # Extract all the regions from it
            for (chrom, start, end), de_fn in tqdm(
                    zip(regions, delta_fns), 
                    total=len(regions), 
                    desc="loci", 
                    leave=False
            ):
                fout = open(de_fn, 'a')
                if first:
                    print("donor\tprotocol\tfull_name\tlength\thap", file=fout)
                try:
                    result = tbx_fh.query(chrom, start, end)
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
                fout.close()
            first = False
    
    if not args.fasta:
        return

    print("Extracing fasta sequences", file=sys.stderr)
    # Build seqs.fa
    with open(args.fasta, 'r') as files:
        # Open fasta files once
        for fasta_fn in tqdm(files.read().strip().split('\n'), desc="seqs"):
            fasta_fh = pysam.FastaFile(fasta_fn)
            mid = SMaHTid.from_re(os.path.basename(fasta_fn))

            # Extract all the regions from it
            for (chrom, start, end), fa_fn in tqdm(
                    zip(regions, fasta_fns), 
                    total=len(regions), 
                    desc="loci", 
                    leave=False
            ):
                fout = open(fa_fn, 'a')
                m_region = f'{chrom}:{start}-{end}_'
                for i in fasta_fh.keys():
                    if i.startswith(m_region):
                        fout.write(f">{mid.full_name}_{i}\n{fasta_fh.fetch(i)}\n")
                fout.close()
