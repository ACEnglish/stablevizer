"""Per-locus tandem repeat instability metrics from a per-read table.

Input: one row per read per TR annotation (cons_start/cons_end), with at least
    locus, region, cons_start, cons_end, period, delta_units, status
Output: one row per (locus, region, cons_start, cons_end).

Design: nothing is normalized across loci or samples. Every scale is either
intrinsic to the locus (its own germline-read noise) or a fixed constant.
"""
import sys
import argparse
import numpy as np
import pandas as pd

GROUP_COLS = ["locus", "region", "cons_start", "cons_end"]

def wilson_lower(k, n, z=1.96):
    """Lower bound of the Wilson score interval for a binomial proportion k/n."""
    if n == 0:
        return 0.0
    p = k / n
    denom = 1 + z**2 / n
    centre = p + z**2 / (2 * n)
    margin = z * np.sqrt(p * (1 - p) / n + z**2 / (4 * n**2))
    return max(0.0, (centre - margin) / denom)


def saturate(x, k=1.0):
    """Fixed bounded map [0, inf) -> [0, 1). x == k gives ~0.63."""
    return 1.0 - np.exp(-np.asarray(x, dtype=float) / k)


def _locus_metrics(g, noise_floor, k_dev, k_spread, z, key_on='is_germ'):
    if key_on == 'is_germ':
        is_som = ~g['is_germ'] #g["status"].isin(["somatic_expansion", "somatic_contraction"])
    else:
        is_som = g['is_somatic']
    som = g.loc[is_som, "delta_units"].to_numpy(dtype=float)
    germ = g.loc[~is_som, "delta_units"].to_numpy(dtype=float)
    n, n_som = len(g), len(som)

    # Locus-intrinsic noise scale (repeat units): robust SD of this locus's own
    # germline reads, floored so loci with perfectly clean germline reads
    # (MAD == 0) don't get an infinite z-score.
    mad_sigma = 1.4826 * np.median(np.abs(germ - np.median(germ))) if len(germ) else 0.0
    sigma = max(noise_floor, mad_sigma)

    vaf = n_som / n if n else np.nan
    vaf_lb = wilson_lower(n_som, n, z)

    if n_som:
        med = float(np.median(som))
        rng = float(som.max() - som.min())
        iqr = float(np.subtract(*np.percentile(som, [75, 25]))) if n_som > 1 else 0.0
        frac_exp = float((som > 0).mean())
    else:
        med = rng = iqr = 0.0
        frac_exp = np.nan

    dev_z = abs(med) / sigma
    spread_z = iqr / sigma
    s_dev = float(saturate(dev_z, k_dev))
    s_spr = float(saturate(spread_z, k_spread))

    # "Noisy-OR": a large shift OR a wide spread each raise the score; a clean
    # single-mode expansion (IQR == 0) is not penalized for having no spread.
    magnitude = 1.0 - (1.0 - s_dev) * (1.0 - s_spr)
    instability = vaf_lb * magnitude

    return pd.Series({
        "n_reads": n,
        "n_somatic": n_som,
        "vaf": vaf,
        "vaf_lb": vaf_lb,
        "frac_expansion": frac_exp,
        "median_delta_units": med,        # signed
        "range_units": rng,               # your original max-min "spread"
        "iqr_units": iqr,                 # robust spread used in the metric
        "noise_sigma_units": sigma,
        "dev_z": dev_z,
        "spread_z": spread_z,
        "burden": vaf_lb * s_dev,         # how much, and how big a shift
        "dispersion": s_spr,              # how many distinct somatic sizes
        "instability": instability,       # single scalar in [0, 1)
        "expected_shift_units": vaf * abs(med),  # raw, interpretable
    })


def instability_metrics(df, noise_floor=1.0, k_dev=1.0, k_spread=1.0, z=1.96,
                        group_cols=GROUP_COLS, key_on='is_germ'):
    """Per-locus instability table.

    noise_floor : minimum germline noise SD, in repeat units
    k_dev, k_spread : saturation scales, in units of noise SD
    z : Wilson interval z (1.96 -> 95%)

    Reads with several non-overlapping TR annotations contribute one row per
    annotation, so each (locus, cons_start, cons_end) is scored independently.
    """
    return (df.groupby(group_cols, sort=False)
              .apply(_locus_metrics, noise_floor, k_dev, k_spread, z, key_on,
                     include_groups=False)
              .reset_index())

def parse_args(args):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("reads", help="read_motifs.tsv")
    parser.add_argument("-o", "--out", default="tr_calls.tsv",
                    help="Output for tsv (%(default)s))")
    args = parser.parse_args(args)
    return args

if __name__ == '__main__':
    args = parse_args(sys.argv[1:])
    try:
        data = pd.read_csv(args.reads, sep='\t')
    except pd.errors.EmptyDataError:
        sys.exit(1)
    result = instability_metrics(data)
    result.to_csv(args.out, sep='\t', index=False)
