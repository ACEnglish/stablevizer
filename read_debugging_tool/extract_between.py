# extract_between.py  (pip install pyfaidx)
import argparse, sys, collections
from pyfaidx import Fasta

COMP = str.maketrans("ACGTacgtNn", "TGCAtgcaNn")
rc = lambda s: s.translate(COMP)[::-1]

def best_hits(paf, min_frac, min_ident):
    best = {}
    for line in open(paf):
        f = line.rstrip("\n").split("\t")
        read, strand, locus = f[0], f[4], f[5]
        qs, qe = int(f[2]), int(f[3])
        tlen, ts, te = int(f[6]), int(f[7]), int(f[8])
        nmatch, blen = int(f[9]), int(f[10])
        if (te - ts) / tlen < min_frac or nmatch / blen < min_ident:
            continue
        key = (read, locus)
        if key not in best or (te - ts) > best[key][3]:
            best[key] = (qs, qe, strand, te - ts)
    return best

ap = argparse.ArgumentParser()
ap.add_argument("reads")
ap.add_argument("left_paf")
ap.add_argument("right_paf")
ap.add_argument("-o", "--inserts", default="inserts.fa")
ap.add_argument("-c", "--clipped", default="clipped.fa")
ap.add_argument("--min-frac", type=float, default=0.8)    # fraction of flank aligned
ap.add_argument("--min-ident", type=float, default=0.75)  # matches / alignment block
ap.add_argument("--min-clip", type=int, default=1)        # min length to report in clipped.fa
ap.add_argument("--keep-flanks", action="store_true")     # include the flank hit in the output
ap.add_argument("--one-based", action="store_true")       # report 1-based inclusive coordinates
args = ap.parse_args()

def coords(a, b):
    """Read coordinates for the header: 0-based half-open by default,
    1-based inclusive with --one-based."""
    return f"{a + 1}-{b}" if args.one_based else f"{a}-{b}"

L = best_hits(args.left_paf, args.min_frac, args.min_ident)
R = best_hits(args.right_paf, args.min_frac, args.min_ident)
reads = Fasta(args.reads, as_raw=True)
n = collections.Counter()

with open(args.inserts, "w") as ins, open(args.clipped, "w") as clip:
    for key in sorted(L.keys() | R.keys()):
        read, locus = key
        seq = str(reads[read])
        rlen = len(seq)

        if key in L and key in R:                         # both flanks: extract insert
            lqs, lqe, ls, _ = L[key]
            rqs, rqe, rs, _ = R[key]
            if ls != rs:
                n["both_orientation_conflict"] += 1
                continue
            if ls == "+":
                if rqs < lqe:
                    n["both_order_conflict"] += 1
                    continue
                a, b = (lqs, rqe) if args.keep_flanks else (lqe, rqs)
                insert = seq[a:b]
            else:
                if lqs < rqe:
                    n["both_order_conflict"] += 1
                    continue
                a, b = (rqs, lqe) if args.keep_flanks else (rqe, lqs)
                insert = rc(seq[a:b])
            ins.write(f">{locus}|{read}|{ls}|{coords(a, b)}|{len(insert)}\n{insert}\n")
            n["insert"] += 1
            continue

        # only one flank hit: output everything from the flank edge to the read end,
        # in reference orientation
        if key in L:
            qs, qe, strand, _ = L[key]
            side = "left"
            if strand == "+":
                a, b = (qs if args.keep_flanks else qe), rlen
                rest = seq[a:b]
            else:
                a, b = 0, (qe if args.keep_flanks else qs)
                rest = rc(seq[a:b])
        else:
            qs, qe, strand, _ = R[key]
            side = "right"
            if strand == "+":
                a, b = 0, (qe if args.keep_flanks else qs)
                rest = seq[a:b]
            else:
                a, b = (qs if args.keep_flanks else qe), rlen
                rest = rc(seq[a:b])

        if len(rest) < args.min_clip:
            n["clipped_too_short"] += 1
            continue
        clip.write(f">{locus}|{read}|{strand}|{side}_only|{coords(a, b)}|{len(rest)}\n{rest}\n")
        n[f"clipped_{side}"] += 1

print(dict(n), file=sys.stderr)
