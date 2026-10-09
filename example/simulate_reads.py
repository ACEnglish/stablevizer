#!/usr/bin/env python3
"""Simulate per-read tandem-repeat allele lengths for donors x tissues.

Output TSV columns: donor, protocol, length, hap
  - each row is one read; protocol == tissue
  - ~COVERAGE reads per donor/tissue (Poisson-distributed)
  - hap is 1 or 2 with p = 0.5
  - each donor has one germline allele length ~ Uniform{80..120}, shared across tissues;
    read lengths = germline + Uniform{-10..10}
  - in one tissue, a random subset of donors carries a low-VAF expansion: a few reads
    (min 3) have lengths Uniform{120..200}
"""
import argparse
import sys

import numpy as np
import pandas as pd
from stablevizer.protocols import PROTOCOLS
import random

random.seed(123)

# ---- defaults (edit here or override on the command line) ----
N_DONORS = 25
TISSUES = random.sample(list(PROTOCOLS.keys()), 10)
N_TISSUES = None          # None -> every donor gets all tissues in TISSUES
COVERAGE = 30             # mean reads per donor/tissue

GERMLINE_RANGE = (80, 120)    # inclusive, uniform per donor
WIGGLE = 10                   # read length = germline +/- WIGGLE (uniform, integer)
EXPANSION_RANGE = (120, 200)  # inclusive, length of expanded reads
EXPANSION_MIN_READS = 5
EXPANSION_DONOR_FRAC = (0.10, 0.50)  # random fraction of donors carrying the expansion
EXPANSION_READ_FRAC = (0.05, 0.15)   # random VAF per expanded donor (floor of MIN_READS reads applies)


def simulate(n_donors, tissues, n_tissues, coverage, expansion_tissue, seed):
    rng = np.random.default_rng(seed)
    n_tissues = len(tissues) if n_tissues is None else n_tissues
    if not 1 <= n_tissues <= len(tissues):
        raise ValueError(f"N_TISSUES must be between 1 and {len(tissues)}")

    width = max(2, len(str(n_donors)))
    donors = [f"donor_{i + 1:0{width}d}" for i in range(n_donors)]
    germline = dict(zip(donors, rng.integers(GERMLINE_RANGE[0], GERMLINE_RANGE[1] + 1, n_donors)))

    # which tissues each donor has
    donor_tissues = {
        d: list(rng.choice(tissues, size=n_tissues, replace=False)) if n_tissues < len(tissues) else list(tissues)
        for d in donors
    }

    # expansion tissue + carrier donors (only donors that actually have that tissue)
    if expansion_tissue is None:
        expansion_tissue = str(rng.choice(tissues))
    eligible = [d for d in donors if expansion_tissue in donor_tissues[d]]
    carriers = set()
    if eligible:
        frac = rng.uniform(*EXPANSION_DONOR_FRAC)
        k = min(len(eligible), max(1, round(frac * len(eligible))))
        carriers = set(rng.choice(eligible, size=k, replace=False))

    rows, truth = [], []
    for d in donors:
        for t in donor_tissues[d]:
            n = max(1, rng.poisson(coverage))
            lengths = germline[d] + np.clip(np.rint(rng.normal(0, WIGGLE / 3, n)), -WIGGLE, WIGGLE).astype(int)
            n_exp = 0
            if t == expansion_tissue and d in carriers:
                vaf = rng.uniform(*EXPANSION_READ_FRAC)
                n_exp = min(n, max(EXPANSION_MIN_READS, round(vaf * n)))
                idx = rng.choice(n, size=n_exp, replace=False)
                lengths[idx] = rng.integers(EXPANSION_RANGE[0], EXPANSION_RANGE[1] + 1, n_exp)
            hap = rng.integers(1, 3, n)
            rows.append(pd.DataFrame({"donor": d, "protocol": t, "length": lengths, "hap": hap}))
            truth.append((d, t, int(germline[d]), n, n_exp))

    reads = pd.concat(rows, ignore_index=True)
    truth = pd.DataFrame(truth, columns=["donor", "protocol", "germline_length", "n_reads", "n_expanded_reads"])
    return reads, truth, expansion_tissue


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--n-donors", type=int, default=N_DONORS)
    ap.add_argument("--tissues", default=",".join(TISSUES), help="comma-separated tissue list")
    ap.add_argument("--n-tissues", type=int, default=N_TISSUES, help="tissues per donor (default: all)")
    ap.add_argument("--coverage", type=float, default=COVERAGE)
    ap.add_argument("--expansion-tissue", default=None, help="tissue carrying the expansion (default: random)")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--out", default="simulated_reads.tsv")
    ap.add_argument("--truth", default=None, help="optional per-donor/tissue ground-truth TSV")
    a = ap.parse_args()

    tissues = [t.strip() for t in a.tissues.split(",") if t.strip()]
    reads, truth, exp_t = simulate(a.n_donors, tissues, a.n_tissues, a.coverage, a.expansion_tissue, a.seed)

    reads.to_csv(a.out, sep="\t", index=False)
    if a.truth:
        truth.to_csv(a.truth, sep="\t", index=False)

    n_carriers = truth.loc[truth.n_expanded_reads > 0, "donor"].nunique()
    print(f"{len(reads)} reads, {a.n_donors} donors, {truth.protocol.nunique()} tissues -> {a.out}", file=sys.stderr)
    print(f"expansion tissue: {exp_t}; carrier donors: {n_carriers}", file=sys.stderr)


if __name__ == "__main__":
    main()
