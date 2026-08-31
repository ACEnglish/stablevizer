"""
Calculate (G+C)%, homopolymer %, and entropy of bed file region sequences
"""
import math
import argparse

import pysam
from tqdm import tqdm
from pastri.ksom.build_kvec import parse_bed_regions, seq_to_kvec

def sequence_entropy(sequence, N=4):
    """
    Computes the Shannon Entropy of a given sequence of a
    biopolymer with `N` possible residues. See (Wooton, 1993)
    for more.

    :param sequence: the nucleotide or protein sequence whose Shannon Entropy is to calculated.
    :param N: the total number of possible residues in the biopolymer `sequence` belongs to.
    """
    encountered_residues = set()
    repvec = []

    for residue in sequence:
        if residue not in encountered_residues:
            residue_count = sequence.count(residue)

            repvec.append(residue_count)

            encountered_residues.add(residue)

        if len(encountered_residues) == N:
            break

    while len(repvec) < N:
        repvec.append(0)
    
    repvec = sorted(repvec, reverse=True)
    L = len(sequence)
    entropy = sum((-1*(n/L)*math.log((n/L), N) for n in repvec if n != 0))

    return entropy

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


def seqstat_bed(args):
    parser = argparse.ArgumentParser(prog="seqstat-bed", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--bed-fn", required=True,
                        help="Input bed file")
    parser.add_argument("-r", "--ref-fn", required=True,
                        help="Input reference fasta")
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="Output `seqstat.bed` %(default)s")
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

