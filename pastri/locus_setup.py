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
import resource
from collections import defaultdict
from bisect import bisect_left, bisect_right

import tabix
import pysam
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
        if not os.path.exists(dest):
            os.mkdir(dest)

        delta_fns.append(os.path.join(dest, 'reads.tsv'))
        fasta_fns.append(os.path.join(dest, 'reads.fasta'))
    
    # ~20k output handles need a higher fd limit than the usual default of 1024
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    target = 65536 if hard == resource.RLIM_INFINITY else min(65536, hard)
    if soft < target:
        resource.setrlimit(resource.RLIMIT_NOFILE, (target, hard))

    print("Extracing qdpi deltas", file=sys.stderr)
    HEADER = "donor\tprotocol\tfull_name\tlength\thap\n"
    # Per-chromosome sorted locus arrays, so each record can find its loci by bisection
    by_chrom = defaultdict(list)
    for idx, (chrom, start, end) in enumerate(regions):
        by_chrom[chrom].append((int(start), int(end), idx))

    starts, ends, idxs = {}, {}, {}
    for chrom, loci in by_chrom.items():
        loci.sort()
        starts[chrom] = [l[0] for l in loci]
        ends[chrom] = [l[1] for l in loci]
        idxs[chrom] = [l[2] for l in loci]
        # bisect trick needs ends sorted too, i.e. no locus nested inside another
        assert all(a <= b for a, b in zip(ends[chrom], ends[chrom][1:])), \
            f"nested loci on {chrom}; use an interval tree instead"

    # Open every locus file once, header included
    outs = [open(fn, 'w') for fn in delta_fns]
    for fh in outs:
        fh.write(HEADER)

    try:
        with open(args.qdpi) as f:
            qdpi_fns = f.read().strip().split('\n')

        for qdpi_fn in tqdm(qdpi_fns, desc="qdpi"):
            mid = SMaHTid.from_re(os.path.basename(qdpi_fn))
            prefix = f"{mid.donor}\t{mid.protocol}\t{mid.full_name}\t"
            buf = defaultdict(list)  # locus index -> text for this sample

            try:
                with pysam.TabixFile(qdpi_fn) as tbx:
                    for chrom in tbx.contigs:
                        if chrom not in starts:
                            continue
                        c_starts, c_ends, c_idxs = starts[chrom], ends[chrom], idxs[chrom]

                        # One sequential pass over the contig
                        for line in tbx.fetch(chrom):
                            _, rs, rend, _, h0, h1, h2 = line.split('\t')
                            rs, rend = int(rs), int(rend)

                            # Loci overlapping [rs, rend): same semantics as tabix query
                            lo = bisect_right(c_ends, rs)
                            hi = bisect_left(c_starts, rend)
                            if lo >= hi:
                                continue

                            span = rend - rs
                            text = ''.join(
                                f"{prefix}{span + int(d)}\t{hap}\n"
                                for hap, h in enumerate((h0, h1, h2)) if h != '.'
                                for d in h.split(',')
                            )
                            if text:
                                for k in range(lo, hi):
                                    buf[c_idxs[k]].append(text)
            except Exception as e:
                tqdm.write(f"ERROR: Unable to parse {qdpi_fn} ({mid.full_name}): {e}", file=sys.stderr)
                continue

            # Flush once per sample
            for k, chunks in buf.items():
                outs[k].writelines(chunks)
    finally:
        for fh in outs:
            fh.close()
      
    if not args.fasta:
        return

    print("Extracing fasta sequences", file=sys.stderr)
    # Matches the "chrom:start-end" prefix of names like "chr1:100-200_readname"
    # (non-greedy so contigs like chrUn_KI270742v1 still work)
    key_re = re.compile(r'^(.+?:\d+-\d+)_')
    
    # Open every locus file once
    out = {
        f'{chrom}:{start}-{end}': open(fa_fn, 'w')
        for (chrom, start, end), fa_fn in zip(regions, fasta_fns)
    }
    try:
        with open(args.fasta) as f:
            sample_fns = f.read().strip().split('\n')

        for fasta_fn in tqdm(sample_fns, desc="seqs"):
            mid = SMaHTid.from_re(os.path.basename(fasta_fn))
            prefix = mid.full_name

            # Stream the file sequentially; no .fai lookups or random access
            with pysam.FastxFile(fasta_fn) as fx:
                for rec in fx:
                    m = key_re.match(rec.name)
                    if m is None:
                        continue
                    fh = out.get(m.group(1))
                    if fh is not None:
                        fh.write(f">{prefix}_{rec.name}\n{rec.sequence}\n")
    finally:
        for fh in out.values():
            fh.close()
