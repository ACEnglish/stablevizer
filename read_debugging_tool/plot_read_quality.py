#!/usr/bin/env python3
"""
plot_read_quality.py   (pip install numpy matplotlib)

Plot per-base Phred quality across the regions written by extract_between.py
(inserts.fa and clipped.fa), using the original reads in FASTQ format.

Coordinates are taken from the fasta headers:
    inserts.fa : locus|read|strand|start-end|length
    clipped.fa : locus|read|strand|left_only (or right_only)|start-end|length
The coordinate convention (0-based half-open vs 1-based inclusive) is detected
per record from the length field, so no flag is needed.

Usage:
    python plot_read_quality.py reads.fastq.gz -i inserts.fa -c clipped.fa
    python plot_read_quality.py reads.fastq.gz --context 200 --locus locus1

Outputs (prefix defaults to read_quality):
    <prefix>.inserts.pdf   one panel per insert record (multi-page)
    <prefix>.clipped.pdf   one panel per clipped record (multi-page)
    <prefix>.overlay.png   all regions overlaid, smoothed quality vs position
    <prefix>.summary.tsv   per-record quality summary (all records, not capped)
"""
import argparse
import collections
import gzip
import os
import random
import sys

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

ROWS, COLS = 4, 2
COLOURS = {"insert": "tab:blue", "clipped": "tab:orange"}
FILEKEY = {"insert": "inserts", "clipped": "clipped"}


def open_text(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_headers(path, kind):
    """Parse records written by extract_between.py. Fields are read from the
    right-hand end so a '|' inside a locus name does not break parsing."""
    recs = []
    with open_text(path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            p = line[1:].rstrip("\n").split("|")
            try:
                if kind == "insert":
                    length, coords, strand, read = int(p[-1]), p[-2], p[-3], p[-4]
                    side, locus = "both", "|".join(p[:-4])
                else:
                    length, coords = int(p[-1]), p[-2]
                    side = p[-3].replace("_only", "")
                    strand, read, locus = p[-4], p[-5], "|".join(p[:-5])
                s, e = (int(x) for x in coords.split("-"))
            except (ValueError, IndexError):
                print(f"WARNING: cannot parse header {line.strip()!r}", file=sys.stderr)
                continue
            if e - s == length:            # 0-based, half-open
                start, end = s, e
            elif e - s + 1 == length:      # 1-based, inclusive
                start, end = s - 1, e
            else:
                print(f"WARNING: coordinates {coords} do not match length {length}: "
                      f"{line.strip()!r}", file=sys.stderr)
                continue
            recs.append(dict(kind=kind, locus=locus, read=read, strand=strand,
                             side=side, start=start, end=end, length=length))
    return recs


def load_quals(fastq, needed):
    """Stream the FASTQ once, keeping quality strings only for reads in `needed`.
    Assumes 4-line FASTQ records. Stops early once every read has been found."""
    quals = {}
    with open_text(fastq) as fh:
        while len(quals) < len(needed):
            h = fh.readline()
            if not h:
                break
            fh.readline()
            fh.readline()
            q = fh.readline().rstrip("\n")
            parts = h[1:].split()
            if parts and parts[0] in needed:
                quals[parts[0]] = q
    return quals


def phred(qstr):
    return np.frombuffer(qstr.encode("ascii"), dtype=np.uint8).astype(np.int16) - 33


def smooth(x, w):
    """Moving average with edge normalisation (no drop-off at the ends)."""
    x = np.asarray(x, dtype=float)
    if w <= 1 or x.size == 0:
        return x
    w = min(w, x.size)
    k = np.ones(w)
    return np.convolve(x, k, "same") / np.convolve(np.ones(x.size), k, "same")


def summarise(qr, lowq):
    if qr.size == 0:
        return dict(mean_Q=np.nan, median_Q=np.nan, mean_err_Q=np.nan,
                    min_Q=np.nan, frac_low=np.nan)
    qf = qr.astype(float)
    return dict(mean_Q=qf.mean(),
                median_Q=float(np.median(qf)),
                mean_err_Q=-10 * np.log10(np.mean(10 ** (-qf / 10))),
                min_Q=qf.min(),
                frac_low=float((qf < lowq).mean()))


def flipped(rec, args):
    return rec["strand"] == "-" and not args.read_orientation


def region_quals(q, rec, args):
    qr = q[rec["start"]:rec["end"]].astype(float)
    return qr[::-1] if flipped(rec, args) else qr


def plot_panel(ax, q, rec, args):
    ctx = args.context
    ws, we = max(0, rec["start"] - ctx), min(q.size, rec["end"] + ctx)
    pos = np.arange(ws, we)
    qq = q[ws:we].astype(float)
    if flipped(rec, args):
        off = (rec["end"] - 1) - pos           # offset 0 = first base of the record
        off, qq = off[::-1], qq[::-1]
    else:
        off = pos - rec["start"]
    ax.plot(off, qq, lw=0.4, color="0.75")
    ax.plot(off, smooth(qq, args.smooth), lw=1.1, color=COLOURS[rec["kind"]])
    ax.axhline(args.lowq, ls=":", lw=0.7, color="red")
    if ctx:
        ax.axvspan(0, rec["length"], color="0.5", alpha=0.10, lw=0)
    ax.set_ylim(0, max(40, qq.max() + 2))
    ax.set_title(f"{rec['locus']}  {rec['read']}\n"
                 f"{rec['strand']} {rec['side']}  read:{rec['start']}-{rec['end']}  "
                 f"mean Q {rec['mean_Q']:.1f}", fontsize=7, loc="left")
    ax.set_xlabel("bp relative to region start", fontsize=7)
    ax.set_ylabel("Phred Q", fontsize=7)
    ax.tick_params(labelsize=7)


def make_pdf(path, sel, quals, args):
    per_page = ROWS * COLS
    with PdfPages(path) as pdf:
        for i in range(0, len(sel), per_page):
            batch = sel[i:i + per_page]
            fig, axes = plt.subplots(ROWS, COLS, figsize=(8.5, 11), squeeze=False)
            for ax, rec in zip(axes.flat, batch):
                plot_panel(ax, quals[rec["read"]], rec, args)
            for ax in list(axes.flat)[len(batch):]:
                ax.axis("off")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


def make_overlay(path, by_kind, quals, args):
    kinds = [k for k in ("insert", "clipped") if by_kind.get(k)]
    if not kinds:
        return
    rng = random.Random(0)
    fig, axes = plt.subplots(len(kinds), 1, figsize=(10, 3.6 * len(kinds)), squeeze=False)
    for ax, kind in zip(axes[:, 0], kinds):
        recs = by_kind[kind]
        shown = rng.sample(recs, args.max_overlay) if len(recs) > args.max_overlay else recs
        alpha = min(0.8, max(0.15, 6 / len(shown)))   # fewer traces -> more opaque
        for rec in shown:
            qr = region_quals(quals[rec["read"]], rec, args)
            step = max(1, qr.size // 2000)
            ax.plot(np.arange(qr.size)[::step], smooth(qr, args.smooth)[::step],
                    lw=0.8, alpha=alpha, color=COLOURS[kind])
        ax.axhline(args.lowq, ls=":", lw=0.7, color="red")
        ax.set_ylim(bottom=0)
        ax.set_xlabel("position in extracted region (bp)")
        ax.set_ylabel(f"Phred Q ({args.smooth} bp mean)")
        ax.set_title(f"{kind} regions: {len(shown)} of {len(recs)} shown")
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("fastq", help="original reads, FASTQ (optionally .gz)")
    ap.add_argument("-i", "--inserts", default="inserts.fa")
    ap.add_argument("-c", "--clipped", default="clipped.fa")
    ap.add_argument("-p", "--prefix", default="read_quality")
    ap.add_argument("--context", type=int, default=0,
                    help="bp of surrounding read to show either side of the region "
                         "(region is shaded when > 0)")
    ap.add_argument("--smooth", type=int, default=25, help="moving-average window (bp)")
    ap.add_argument("--lowq", type=float, default=10, help="low-quality threshold (Phred)")
    ap.add_argument("--max-plots", type=int, default=200,
                    help="max panels per PDF (summary TSV always has every record)")
    ap.add_argument("--max-overlay", type=int, default=300,
                    help="max traces drawn per overlay panel")
    ap.add_argument("--locus", help="only plot records for this locus")
    ap.add_argument("--read-orientation", action="store_true",
                    help="plot in the read's own orientation instead of the orientation "
                         "of the fasta record (minus-strand records are flipped by default)")
    args = ap.parse_args()

    recs = []
    for kind, path in (("insert", args.inserts), ("clipped", args.clipped)):
        if os.path.exists(path):
            recs += parse_headers(path, kind)
        else:
            print(f"note: {path} not found, skipping {kind} records", file=sys.stderr)
    if args.locus:
        recs = [r for r in recs if r["locus"] == args.locus]
    if not recs:
        sys.exit("No records to plot.")

    raw = load_quals(args.fastq, {r["read"] for r in recs})
    quals = {name: phred(s) for name, s in raw.items()}

    kept, skipped = [], collections.Counter()
    for r in recs:
        q = quals.get(r["read"])
        if q is None:
            skipped["read_not_in_fastq"] += 1
            continue
        if r["start"] < 0 or r["end"] > q.size:
            skipped["coords_outside_read"] += 1
            continue
        r.update(summarise(q[r["start"]:r["end"]], args.lowq))
        kept.append(r)

    # summary table (every record)
    na = lambda v, f: "NA" if v != v else format(v, f)
    with open(f"{args.prefix}.summary.tsv", "w") as out:
        out.write("\t".join(["type", "locus", "read", "strand", "side", "read_start",
                             "read_end", "length", "mean_Q", "median_Q", "mean_err_Q",
                             "min_Q", "frac_below_lowq"]) + "\n")
        for r in kept:
            out.write("\t".join([r["kind"], r["locus"], r["read"], r["strand"], r["side"],
                                 str(r["start"]), str(r["end"]), str(r["length"]),
                                 na(r["mean_Q"], ".2f"), na(r["median_Q"], ".1f"),
                                 na(r["mean_err_Q"], ".2f"), na(r["min_Q"], ".0f"),
                                 na(r["frac_low"], ".3f")]) + "\n")

    by_kind = {}
    for kind in ("insert", "clipped"):
        sel = [r for r in kept if r["kind"] == kind and r["length"] > 0]
        skipped[f"empty_{kind}"] += sum(1 for r in kept if r["kind"] == kind and r["length"] == 0)
        if not sel:
            continue
        by_kind[kind] = sel
        make_pdf(f"{args.prefix}.{FILEKEY[kind]}.pdf", sel[:args.max_plots], quals, args)
        if len(sel) > args.max_plots:
            print(f"note: {kind}: plotted first {args.max_plots} of {len(sel)} records "
                  f"(raise --max-plots for more)", file=sys.stderr)
    make_overlay(f"{args.prefix}.overlay.png", by_kind, quals, args)

    print(f"{len(kept)} records summarised; skipped: {dict(skipped)}", file=sys.stderr)


if __name__ == "__main__":
    main()
