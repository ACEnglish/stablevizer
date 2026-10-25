#!/usr/bin/env python3
"""
Harmonized tandem-repeat annotation across all reads in a locus MSA, then
per-read copy number + germline/somatic classification.

Strategy (why it's built this way)
----------------------------------
Running a TR detector independently on every read gives inconsistent motif
calls/boundaries, and low-VAF expansion reads are exactly the ones that get
mis-annotated. So instead:

  1. Build a majority-rule consensus from the MSA (columns that are mostly gap,
     i.e. the low-VAF insertions, are dropped -> consensus ~ germline).
  2. Detect TRs ONCE on that consensus (TRF by default, TRAL optionally).
     -> one harmonized set of regions + motifs per locus.
  3. Map each consensus TR back to MSA columns, so every read gets the
     residues falling in those columns (including insertion columns) as its
     "TR segment". Same boundaries, same motif, for every read.
  4. For each read/segment: copy number of the harmonized motif (edlib infix
     alignment against motif^N; falls back to length/period), purity, length.
  5. Call germline alleles (high-support length clusters, up to ploidy) and flag
     reads outside them as somatic expansion/contraction, with delta in units.

Usage
-----
  pastri anno_msa msas/ -o tr_calls.tsv 

Deps: pandas; edlib (pip install edlib); TRF binary on PATH for trf;
"""
import argparse
import math
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path

import pandas as pd
from tqdm import tqdm

MSA_EXTS = {".fa", ".fasta", ".fna", ".aln", ".afa", ".msa"}


# ----------------------------------------------------------------------------
# IO / MSA helpers
# ----------------------------------------------------------------------------
def read_fasta(path):
    names, seqs, cur = [], [], []
    for line in Path(path).read_text().splitlines():
        if line.startswith(">"):
            if names:
                seqs.append("".join(cur).upper())
            names.append(line[1:].split()[0])
            cur = []
        elif line.strip():
            cur.append(line.strip())
    if names:
        seqs.append("".join(cur).upper())
    assert len({len(s) for s in seqs}) == 1, f"{path}: sequences are not aligned (unequal lengths)"
    return names, seqs


def msa_consensus(seqs, gap="-"):
    """Majority consensus. Columns where gaps are >50% are dropped.
    Returns (consensus_str, col_of_cons) where col_of_cons[i] = MSA column of consensus base i."""
    n = len(seqs)
    cons, cols = [], []
    for c in range(len(seqs[0])):
        counts = Counter(s[c] for s in seqs)
        if counts.get(gap, 0) > n / 2:
            continue
        counts.pop(gap, None)
        cons.append(counts.most_common(1)[0][0])
        cols.append(c)
    return "".join(cons), cols


def ungapped_prefix_counts(aligned, gap="-"):
    """cum[c] = number of non-gap residues in columns < c (length ncol+1)."""
    cum = [0]
    for ch in aligned:
        cum.append(cum[-1] + (ch != gap))
    return cum


# ----------------------------------------------------------------------------
# Motif utilities
# ----------------------------------------------------------------------------
def primitive_root(motif):
    """ACAC -> AC"""
    n = len(motif)
    for p in range(1, n + 1):
        if n % p == 0 and motif[:p] * (n // p) == motif:
            return motif[:p]
    return motif


def revcomp(s):
    return s[::-1].translate(str.maketrans("ACGTN", "TGCAN"))


def canonical_motif(motif):
    """Name a motif by lexicographically smallest rotation (strand-specific; reads are
    assumed to be in the same orientation). Reverse-complement is reported separately."""
    motif = primitive_root(motif)
    rots = [motif[i:] + motif[:i] for i in range(len(motif))]
    return min(rots)


# ----------------------------------------------------------------------------
# TR detection on the consensus (one harmonized annotation per locus)
# ----------------------------------------------------------------------------
def detect_trs_trf(seq, name, trf="trf", params=(2, 7, 7, 80, 10, 50, 500)):
    """Tandem Repeats Finder. Returns list of dicts (0-based half-open on `seq`)."""
    with tempfile.TemporaryDirectory() as d:
        fa = Path(d) / "cons.fa"
        fa.write_text(f">{name}\n{seq}\n")
        # TRF returns a non-zero exit status by design; don't check
        subprocess.run([trf, fa.name, *map(str, params), "-d", "-h"], cwd=d,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, check=False)
        dat = Path(d) / (fa.name + "." + ".".join(map(str, params)) + ".dat")
        hits = []
        if dat.exists():
            for line in dat.read_text().splitlines():
                p = line.split()
                # start end period copies consize matches indels score A C G T entropy motif repeat ...
                if len(p) >= 15 and p[0].isdigit() and p[1].isdigit():
                    hits.append(dict(start=int(p[0]) - 1, end=int(p[1]), period=int(p[2]),
                                     score=int(p[7]), motif=p[13], source="TRF"))
        return hits


def harmonize_regions(hits, max_overlap=0.5):
    """Keep best-scoring non-overlapping TR calls; collapse each to primitive motif."""
    kept = []
    for h in sorted(hits, key=lambda h: -h["score"]):
        ok = True
        for k in kept:
            ov = min(h["end"], k["end"]) - max(h["start"], k["start"])
            if ov > max_overlap * min(h["end"] - h["start"], k["end"] - k["start"]):
                ok = False
                break
        if ok:
            h = dict(h)
            h["motif"] = canonical_motif(h["motif"])
            h["period"] = len(h["motif"])
            kept.append(h)
    return sorted(kept, key=lambda h: h["start"])


# ----------------------------------------------------------------------------
# Per-read copy number
# ----------------------------------------------------------------------------
def count_copies(seg, motif):
    """Return (copy_number, edit_distance|None). Infix-aligns the read segment to motif^N."""
    p = len(motif)
    if not seg:
        return 0.0, 0
    try:
        import edlib
    except ImportError:
        return len(seg) / p, None
    target = motif * (len(seg) // p + 3)
    r = edlib.align(seg, target, mode="HW", task="locations")
    s, e = r["locations"][0]
    return (e - s + 1) / p, r["editDistance"]


# ----------------------------------------------------------------------------
# Germline vs somatic calling for one region
# ----------------------------------------------------------------------------
def call_alleles(lengths, period, ploidy=2, min_allele_frac=0.15, tol_bp=None):
    """Greedy length clustering -> germline allele lengths (<= ploidy)."""
    tol = tol_bp if tol_bp is not None else max(2, period - 1)
    remaining = list(lengths)
    total = len(lengths)
    alleles = []
    while remaining and len(alleles) < ploidy:
        mode_len, _ = Counter(remaining).most_common(1)[0]
        members = [x for x in remaining if abs(x - mode_len) <= tol]
        if len(members) / total < min_allele_frac:
            break
        alleles.append(int(round(sum(members) / len(members))))
        remaining = [x for x in remaining if abs(x - mode_len) > tol]
    return alleles, tol


def classify(seg_len, alleles, period, tol):
    if not alleles:
        return "NA", None, float("nan")
    nearest = min(alleles, key=lambda a: abs(seg_len - a))
    delta_bp = seg_len - nearest
    delta_units = delta_bp / period
    if abs(delta_bp) <= tol:
        return "germline", nearest, delta_units
    return ("somatic_expansion" if delta_bp > 0 else "somatic_contraction"), nearest, delta_units


# ----------------------------------------------------------------------------
# Locus driver
# ----------------------------------------------------------------------------
def annotate_locus(names, seqs, locus, trf="trf", ploidy=1, min_allele_frac=0.15):
    cons, col_of_cons = msa_consensus(seqs)
    if not cons:
        return []

    detect = lambda s, n: detect_trs_trf(s, n, trf=trf)
    regions = harmonize_regions(detect(cons, locus + "_consensus"))
    if not regions:
        print(f"[{locus}] no TR detected in consensus", file=sys.stderr)
        return []

    cums = [ungapped_prefix_counts(s) for s in seqs]
    ungapped = ["".join(ch for ch in s if ch != "-") for s in seqs]

    rows = []
    for ri, reg in enumerate(regions):
        c0, c1 = col_of_cons[reg["start"]], col_of_cons[reg["end"] - 1]  # inclusive MSA columns
        per_read = []
        for name, u, cum in zip(names, ungapped, cums):
            seg = u[cum[c0]:cum[c1 + 1]]  # includes insertion columns inside the TR span
            cn, ed = count_copies(seg, reg["motif"])
            purity = (1 - ed / len(seg)) if (ed is not None and seg) else float("nan")
            per_read.append((name, seg, cn, purity))

        alleles, tol = call_alleles([len(s) for _, s, _, _ in per_read], reg["period"],
                                    ploidy=ploidy, min_allele_frac=min_allele_frac)
        for name, seg, cn, purity in per_read:
            status, allele, d_units = classify(len(seg), alleles, reg["period"], tol)
            rows.append(dict(
                locus=locus, region=ri,
                cons_start=reg["start"], cons_end=reg["end"],
                motif=reg["motif"], motif_rc=revcomp(reg["motif"]), period=reg["period"],
                read=name, tr_len_bp=len(seg), copy_number=round(cn, 2), purity=round(purity, 3),
                germline_alleles_bp=",".join(map(str, alleles)),
                nearest_allele_bp=allele, delta_units=round(d_units, 2) if d_units == d_units else None,
                status=status,
            ))
    return rows


def anno_msa_main(args):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("msa", help="MSA FASTA file, or a directory of them")
    ap.add_argument("-o", "--out", default="tr_calls.tsv", help="Output tsv (%(default)s))")
    ap.add_argument("--trf", default="trf", help="path to TRF binary (%(default)s)")
    ap.add_argument("--ploidy", type=int, default=1, help="Germline ploidy (%(default)s)")
    ap.add_argument("--min-allele-frac", type=float, default=0.15,
                    help="a length cluster needs this read fraction to count as a germline allele (%(default)s)")
    args = ap.parse_args(args)

    p = Path(args.msa)
    files = sorted(f for f in p.iterdir() if f.suffix in MSA_EXTS) if p.is_dir() else [p]
    rows = []
    for f in tqdm(files):
        names, seqs = read_fasta(f)
        locus = Path(f).stem
        rows += annotate_locus(names, seqs, locus, args)

    df = pd.DataFrame(rows)
    df.to_csv(args.out, sep="\t", index=False)
    print(f"wrote {len(df)} read-level rows -> {args.out}")

    if df.empty:
        return
    # Example queries ---------------------------------------------------------
    som = df[df.status.str.startswith("somatic")]
    print("\nSomatic reads per locus/region:")
    summ = (df.groupby(["locus", "region", "motif"])
              .agg(n_reads=("read", "size"),
                   n_somatic=("status", lambda s: s.str.startswith("somatic").sum()))
              .assign(vaf=lambda x: (x.n_somatic / x.n_reads).round(3)))
    print(summ[summ.n_somatic > 0].sort_values(by='n_somatic', ascending=False).head(30).to_string())
    #print("\nSomatic expansions (motif, copy number, delta units):")
    #print(som.groupby(['locus', 'region', 'status','motif'])['delta_units'].describe())

if __name__ == "__main__":
    anno_msa_main(sys.argv[1:])
