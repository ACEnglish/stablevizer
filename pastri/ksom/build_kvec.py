"""
Build kmer vectors from reference regions
"""
import argparse

import pysam
import numpy as np
from tqdm import tqdm

## LLM rewrite of kanpig code
# Lookup table: maps ASCII byte value -> 2-bit code (A=0, G=1, C=2, T=3, default=0)
_NUC_LUT = np.zeros(256, dtype=np.uint64)
_NUC_LUT[ord('A')] = 0
_NUC_LUT[ord('a')] = 0
_NUC_LUT[ord('G')] = 1
_NUC_LUT[ord('g')] = 1
_NUC_LUT[ord('C')] = 2
_NUC_LUT[ord('c')] = 2
_NUC_LUT[ord('T')] = 3
_NUC_LUT[ord('t')] = 3

COMPLEMENT = str.maketrans("ATCG", "TAGC")

def parse_bed_regions(fn):
    """
    Simple bed parser
    """
    ret = []
    with open(fn, 'r') as fh:
        for line in fh:
            chrom, start, end = line.strip().split('\t')[:3]
            start = int(start)
            end = int(end)
            ret.append((chrom, start, end))
    return ret

def seq_to_kvec(seq, k_size=4, kmers=None, count=False):
    """
    Build the forward/rc kmer vector
    Provide `kmers` (a single 4**k_size) to place into an existing array
    Default reports the frequency, count=True reports the count
    """
    if kmers is None:
        kmers = np.zeros(4 ** k_size, dtype=np.float16)

    # Forward
    pos, cnt = seq_to_kmer(seq.encode(), k_size)
    kmers[pos] = cnt
    # Reverse compliment
    pos, cnt = seq_to_kmer(rev_comp(seq).encode(), k_size)
    kmers[pos] = cnt
    
    if not count:
        total = (len(seq) - k_size + 1) * 2
        kmers /= total
    
    return kmers

def seq_to_kmer(sequence: bytes, kmer: int = 4):
    """
    Vectorized k-mer encoder.

    sequence: bytes or bytearray of nucleotides
    kmer: k-mer length
    negative: if True, counts are -1.0, else +1.0

    Returns: (keys, values) as numpy arrays:
        keys   -> uint64 encoded k-mers, sorted ascending
        values -> float32 aggregated counts (zeros removed)
    """
    seq = np.frombuffer(sequence, dtype=np.uint8)
    n = seq.shape[0]

    if n < kmer:
        return np.empty(0, dtype=np.uint64), np.empty(0, dtype=np.float32)

    # Encode every base to its 2-bit code in one vectorized lookup
    codes = _NUC_LUT[seq]  # uint64 array, shape (n,)

    num_kmers = n - kmer + 1

    # Build a rolling 2-bit-packed representation using a prefix trick:
    # Compute cumulative "rolling" value via strided sliding window sum
    # of codes[i] << shift, using the standard rolling-hash approach:
    #
    #   value[0] = sum_{j=0}^{k-1} codes[j] << (2*(k-1-j))
    #   value[i] = ((value[i-1] & mask) << 2) + codes[i+k-1]   for i>=1
    #
    # Vectorize by building it with a strided-sum via as_strided (fast, no python loop).

    shifts = (2 * np.arange(kmer - 1, -1, -1)).astype(np.uint64)  # shape (k,)

    # sliding_window_view gives an (num_kmers, k) view without copying data
    windows = np.lib.stride_tricks.sliding_window_view(codes, kmer)  # (num_kmers, k) uint64

    values = np.left_shift(windows, shifts).sum(axis=1, dtype=np.uint64)

    counts = np.full(num_kmers, 1.0, dtype=np.float32)

    # Sort by k-mer key
    order = np.argsort(values, kind='stable')
    sorted_keys = values[order]
    sorted_counts = counts[order]

    # Aggregate (sum) counts for identical keys, analogous to dedup_by + sum
    unique_keys, start_idx = np.unique(sorted_keys, return_index=True)
    summed_counts = np.add.reduceat(sorted_counts, start_idx)

    # Retain nonzero counts (mirrors retain(|&(_, v)| v != 0.0))
    nonzero_mask = summed_counts != 0.0
    final_keys = unique_keys[nonzero_mask]
    final_values = summed_counts[nonzero_mask]

    return final_keys, final_values

def rev_comp(seq):
    return seq.translate(COMPLEMENT)[::-1]

def build_kvec(args):
    """
    Main
    """
    parser = argparse.ArgumentParser(prog="kvec", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--bed-fn", required=True,
                        help="Bed file of TR regions")
    parser.add_argument("-r", "--ref-fn", required=True,
                        help="Reference fasta file")
    parser.add_argument("-o", "--out-fn", required=True,
                        help="Output kmervectors.npz")
    parser.add_argument("-k", "--kmer", type=int, default=4,
                        help="Kmer size (%(default)s)")
    parser.add_argument("-c", "--count", action="store_true",
                        help="Report kmer count instead of frequency")
    args = parser.parse_args(args)

    regions = parse_bed_regions(args.bed_fn)
    ref = pysam.FastaFile(args.ref_fn)

    kmers = np.zeros((len(regions), 4 ** args.kmer), dtype=np.float16)
    idx = 0
    for reg in tqdm(regions):
        seq_to_kvec(ref.fetch(*reg),
                    args.kmer,
                    kmers[idx],
                    args.count)
        idx += 1

    np.savez(args.out_fn,
             kmers=kmers,
             bed_fn=args.bed_fn,
             ref_fn=args.ref_fn,
             k=args.kmer,
             count=args.count)
