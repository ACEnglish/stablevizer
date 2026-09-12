"""
Select loci from QDPI files that show signs of instability
Assumes inputs have same loci in same order
"""
import os
import sys
import gzip
import argparse
from collections import defaultdict

import numpy as np
from pastri import SMaHTid

EMPTY = np.array([], dtype=int)

def to_int(data):
    """
    String to numpy
    """
    if data == '.':
        return EMPTY
    return np.fromstring(data, dtype=int, sep=",")


def tissue_dev(deltas, smhtids, global_median, spread_min, min_alt, min_vaf):
    """
    Consolidate read support by tissue and check if there are at least min_alt reads in the tissue deviating from the
    median and they are above some minimum VAF
    """
    by_tiss = defaultdict(list)
    for m_id, p in zip(smhtids, deltas):
        if p.size:
            by_tiss[m_id.tissue_abv].append(p)

    for ds in by_tiss.values():
        ds = np.concatenate(ds) if len(ds) > 1 else ds[0]
        # How many reads outside bounds
        d_count = int(((ds > global_median + spread_min) |
                       (ds < global_median - spread_min)).sum())
        if d_count >= min_alt and (d_count / len(ds) >= min_vaf):
            return True  # Only need one to pass
    return False


def parse_args(args):
    """
    UI
    """
    parser = argparse.ArgumentParser(prog="locus_select", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("inputs", type=str, nargs="+",
                        help="Read lengths tsv")
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="Output loci (stdout)")
    parser.add_argument("-c", "--cov", dtype=int, default=30,
                        help="Minimum haplotype coverage (%(default)s)")
    parser.add_argument("-a", "--alt", dtype=int, default=3,
                        help="Minimum deviating coverage (%(default)s)")
    parser.add_argument("-v", "--vaf", dtype=float, default=0.01,
                        help="Minimum deviating VAF (%(default)s)")
    parser.add_argument("-s", "--spread", dtype=int, default=10,
                        help="Minimum spread of deltas (%(default)s)")
    return parser.parse_args(args)


def selector_main(args):
    """
    Main
    """
    args = parse_args(args)

    files = [gzip.open(_) for _ in args.inputs]
    smhtids = [SMaHTid.from_re(os.path.basename(_)) for _ in args.inputs]
    args.output = open(args.output, 'w')
    while True:
        # Join each qdpi bed file locus
        try:
            lines = [next(_).decode() for _ in files]
        except StopIteration:
            break

        locus = None
        h1_parts = []
        h2_parts = []
        for line in lines:
            data = line.strip().split('\t')
            locus = data[:3]
            h1_parts.append(to_int(data[5]))
            h2_parts.append(to_int(data[6]))

        for parts in [h1_parts, h2_parts]:
            h = np.concatenate(parts) if len(parts) > 1 else parts[0]
            # Sufficiently covered overall
            if h.size < args.cov:
                continue

            # Global spread
            if h.max() - h.min() < args.spread:
                continue

            # Minimum deviating coverage and minimum VAF of that deviating coverage
            if not tissue_dev(parts, smhtids, np.median(h), args.spread, args.alt, args.vaf):
                continue

            # Only need one to pass
            print(*locus, sep='\t', file=args.output)
            break

    args.output.close()
