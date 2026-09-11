import os
import sys
import gzip
from collections import defaultdict

import numpy as np
from pastri import SMaHTid

# Should parameterize
MIN_TOT_COV = 30 # Only analyze haplotypes with at least this much coverage
SPREAD_MIN = 10 # Need the deltas to span at least this far
MIN_COV = 3 # We need at least this many reads IN A TISSUE that are SPREAD_MIN from global median size
MIN_VAF = 0.01 # And IN THAT TISSUE it is at least this VAF
EMPTY = np.array([], dtype=int)

def to_int(data):
    """
    String to numpy
    """
    if data == '.':
        return EMPTY
    return np.fromstring(data, dtype=int, sep=",")

def parse_line(line):
    """
    Get the locus fields and the two haps deltas to int
    """
    fields = line.strip().split('\t')
    return fields[:3], to_int(fields[5]), to_int(fields[6])

def spread(d):
    """
    total spread of deltas
    """
    return (d.max() - d.min()) if d.size else 0

def dev_count(d, median=None):
    """
    Are there at least MIN_COV reads at least SPREAD_MIN away from the mean
    """
    if not d.size:
        return 0
    if not median:
        median = np.median(d)
    return int(((d > median + SPREAD_MIN) | (d < median - SPREAD_MIN)).sum())

def dev_count_mad(d, mad_k=3):
    """
    Ignore, test that explodes number of selected loci to like everything
    Alternative to dev_count: MAD-based robust threshold instead of a fixed
    SPREAD_MIN offset from the median. Not currently called - swap in for
    dev_count() if needed. See mad_k note below.
    """
    if not d.size:
        return 0
    m = np.median(d)
    mad = np.median(np.abs(d - m))
    if mad == 0:
        return int((d != m).sum())
    scaled_mad = mad * 1.4826
    return int((np.abs(d - m) > mad_k * scaled_mad).sum())

def by_tissue_vaf_check(parts, smhtids, global_median):
    """
    Consolidate read support by tissue and check if there are at least MIN_COV reads in the tissue deviating from the
    median and they are above some minimum VAF
    """
    separated = defaultdict(list)
    for m_id, p in zip(smhtids, parts):
        if not p.size:
            continue
        separated[m_id.tissue_abv].append(p)

    for k, deltas in separated.items():
        deltas = np.concatenate(deltas) if len(deltas) > 1 else deltas[0]
        d_count = dev_count(deltas, global_median)
        if d_count >= MIN_COV and (d_count / len(deltas) >= MIN_VAF):
            return True
    return False

if __name__ == '__main__':
    files = [gzip.open(_) for _ in sys.argv[1:]]
    smhtids = [SMaHTid.from_re(os.path.basename(_)) for _ in sys.argv[1:]]
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
            fields, d1, d2 = parse_line(line)
            if locus is None:
                locus = fields
            h1_parts.append(d1)
            h2_parts.append(d2)

        h1 = np.concatenate(h1_parts) if len(h1_parts) > 1 else h1_parts[0]
        h2 = np.concatenate(h2_parts) if len(h2_parts) > 1 else h2_parts[0]
        passing = False
        for h, parts in zip([h1, h2], [h1_parts, h2_parts]):
            # Sufficiently covered overall
            if h.size < MIN_TOT_COV:
                continue

            # Global spread
            if spread(h) < SPREAD_MIN:
                continue
            
            # Minimum deviating coverage and minimum VAF of that deviating coverage
            if not by_tissue_vaf_check(parts, smhtids, np.median(h)):
                continue

            passing = True
            break # Only need one to pass
        
        if not passing:
            continue

        print(*locus, sep='\t')
