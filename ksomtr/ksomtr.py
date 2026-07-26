import sys
import pickle
import argparse

import pysam
import minisom
from truvari.annotations.lcr import sequence_entropy
import numpy as np
from tqdm import tqdm

KMER = 3
COMPLEMENT = str.maketrans("ATCG", "TAGC")

import numpy as np

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


def seq_to_kmer(sequence: bytes, kmer: int = KMER):
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

def rev_comp(seq):
    return seq.translate(COMPLEMENT)[::-1]

def seq_to_kvec(seq, k_size=KMER, kmers=None, normalize=False):
    """
    Build the forward/rc normalized kmer vector
    Provide `kmers` (a single 4**k_size) to place into an existing array
    """
    if kmers is None:
        kmers = np.zeros(4 ** k_size, dtype=np.float16)


    pos, cnt = seq_to_kmer(seq.encode(), k_size)
    kmers[pos] = cnt
    pos, cnt = seq_to_kmer(rev_comp(seq).encode(), k_size)
    
    if normalize:
        total = (len(seq) - k_size + 1) * 2
        kmers /= total
    
    return kmers

def build_kvecs(args):
    """
    Building the initial kmer vectors
    """
    parser = argparse.ArgumentParser(prog="kvec")
    parser.add_argument("-b", "--bed-fn", required=True)
    parser.add_argument("-r", "--ref-fn", required=True)
    parser.add_argument("-o", "--out-fn", required=True)
    parser.add_argument("-k", "--kmer", type=int, default=KMER,
                        help="Kmer size (%(default)s)")
    parser.add_argument("-n", "--normalize", action="store_true",
                        help="Normalize")
    args = parser.parse_args(args)

    regions = parse_bed_regions(args.bed_fn)
    ref = pysam.FastaFile(args.ref_fn)

    kmers = np.zeros((len(regions), 4 ** args.kmer), dtype=np.float16)
    idx = 0
    for reg in tqdm(regions):
        seq = ref.fetch(*reg)

        # This now needs to be checked
        seq_to_kvec(seq, args.kmer, kmers[idx], args.normalize)
        idx += 1

    np.savez(args.out_fn,
             kmers=kmers,
             bed_fn=args.bed_fn,
             ref_fn=args.ref_fn,
             k=args.kmer,
             normalize=args.normalize)

def build_som(args):
    """
    Train the SOM on the kmers
    """
    parser = argparse.ArgumentParser(prog="som")
    parser.add_argument("kmers_fn")
    parser.add_argument("-o", "--out-fn", required=True)
    parser.add_argument("-s", "--sigma", type=float, default=1.5,
                        help="%(default)s")
    parser.add_argument("-l", "--learning-rate", type=float, default=1,
                        help="%(default)s")
    parser.add_argument("-i", "--iters", type=int, default=1_000_000,
                        help="%(default)s")
    parser.add_argument("-p", "--pca-init", action="store_true")
    parser.add_argument("--seed", default=None)
    args = parser.parse_args(args)

    data = np.load(args.kmers_fn)
    som = minisom.MiniSom(25, 25, data['kmers'].shape[1],
                          sigma=args.sigma,
                          learning_rate=args.learning_rate,
                          topology='hexagonal',
                          neighborhood_function='gaussian',
                          activation_distance='euclidean',
                          random_seed=args.seed
                         )

    if args.pca_init:
        som.pca_weights_init(data['kmers'])

    som.train_batch(data['kmers'], args.iters, verbose=True)
    
    output = {'som': som,
              'k': int(data['k']),
              'normalize': bool(data['normalize']),
              }

    with open(args.out_fn, 'wb') as fout:
        pickle.dump(output, fout)

def bed_map(args):
    """
    Map a bed file onto a som
    """
    parser = argparse.ArgumentParser(prog="bed-map")
    parser.add_argument("-b", "--bed-fn", required=True)
    parser.add_argument("-r", "--ref-fn", required=True)
    parser.add_argument("-s", "--som-fn", required=True)
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="%(default)s")

    args = parser.parse_args(args)

    regions = parse_bed_regions(args.bed_fn)
    ref = pysam.FastaFile(args.ref_fn)
    som = pickle.load(open(args.som_fn, 'rb'))
    fout = open(args.output, 'w')

    for region in tqdm(regions):
        seq = ref.fetch(*region)
        vec = seq_to_kvec(seq, som['k'], normalize=som['normalize'])
        w = som['som'].winner(vec)
        print(*region, *w, sep='\t', file=fout)
    fout.close()

def homopolymer_percent(seq, min_run=3, ignore_n=True):
    """
    Calculate the percentage of a nucleotide sequence that lies within
    homopolymer runs (consecutive repeats of the same base).

    Parameters
    ----------
    seq : str
        Nucleotide sequence (A, C, G, T, N, case-insensitive).
    min_run : int, default 2
        Minimum run length to count as a "homopolymer"
        (2 means any repeated pair counts; 3+ is a stricter definition
        commonly used for flagging problematic runs).
    ignore_n : bool, default True
        If True, 'N' bases are excluded from both the homopolymer count
        and the total length (since N is ambiguous, not a real repeat).
        If False, N's are treated like any other base and runs of N
        count as homopolymers too.

    Returns
    -------
    float
        Percentage (0-100) of the sequence that is part of a homopolymer run.
    """
    seq = seq.upper()

    if ignore_n:
        effective_seq = seq.replace('N', '')
    else:
        effective_seq = seq

    total_len = len(effective_seq)
    if total_len == 0:
        return 0.0

    homopolymer_bases = 0
    i = 0
    n = len(effective_seq)

    while i < n:
        j = i
        while j < n and effective_seq[j] == effective_seq[i]:
            j += 1
        run_length = j - i
        if run_length >= min_run:
            homopolymer_bases += run_length
        i = j

    return homopolymer_bases / total_len

def bed_seqstat(args):
    """
    Calculate GC% of bed file regions
    """
    parser = argparse.ArgumentParser(prog="bed-seqstat")
    parser.add_argument("-b", "--bed-fn", required=True)
    parser.add_argument("-r", "--ref-fn", required=True)
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="%(default)s")
    args = parser.parse_args(args)

    regions = parse_bed_regions(args.bed_fn)
    ref = pysam.FastaFile(args.ref_fn)
    fout = open(args.output, 'w')

    for region in tqdm(regions):
        seq = ref.fetch(*region).upper()
        gc = seq.count('G') + seq.count('C')
        gc_pct = gc / len(seq)
        hom_pct = homopolymer_percent(seq)
        entropy = sequence_entropy(seq)
        print(*region, "%.5f" % (gc / len(seq)),
                       "%.5f" % (hom_pct),
                       "%.5f" % (entropy),
              sep='\t', file=fout)
    fout.close()

def sample_farthest_point(df, cols, n_samples, random_state=None):
    rng = np.random.default_rng(random_state)
    points = df[cols].to_numpy()
    n = len(points)

    mask = np.ones(n, dtype=bool)
    min_dists = np.full(n, np.inf)

    start = rng.integers(n)
    selected_idx = [start]
    mask[start] = False
    tot = n_samples - 1
    for _ in tqdm(range(tot), total=tot):
        last_point = points[selected_idx[-1]]
        dists = np.linalg.norm(points - last_point, axis=1)
        min_dists = np.minimum(min_dists, dists)

        masked_dists = np.where(mask, min_dists, -np.inf)
        next_idx = np.argmax(masked_dists)

        selected_idx.append(next_idx)
        mask[next_idx] = False

    return df.iloc[selected_idx]

def downsample(args):
    parser = argparse.ArgumentParser(prog="downsample")
    parser.add_argument("-b", "--bed-fn", required=True)
    parser.add_argument("-n", "--num", required=True)
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="%(default)s")
    args = parser.parse_args(args)

    # Should be asserting that this is the shape..
    df = pd.read_csv(in_fn, sep='\t', names=['chrom', 'start', 'end', 'a', 'b', 'c'])
    fps_sample = sample_farthest_point(df, ['a', 'b', 'c'], n_samples=NSAMP, random_state=8811)
    fps_sample.to_csv(args.output, sep='\t', index=False, header=False)

if __name__ == '__main__':

    TOOLS = {
            'kvec': build_kvecs,
            'som': build_som,
            'bed-map': bed_map,
            'bed-seqstat': bed_seqstat,
            'downsample': downsample,
            }
    parser = argparse.ArgumentParser(prog="ksomtr", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str, default=None,
                        help="Command to execute: "  + "\n".join(TOOLS.keys()))
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")
    args = parser.parse_args()
    TOOLS[args.cmd](args.options)
