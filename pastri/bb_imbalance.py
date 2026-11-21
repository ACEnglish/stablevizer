"""
Allelic imbalance testing with a beta-binomial null.

Input : a table with one row per sample x locus and read counts for each haplotype.
Output: per-locus allelic fraction, beta-binomial p/q-values, an empirical-Bayes
        (shrunken) estimate with credible interval, and an imbalance call.

Requires: numpy, pandas, scipy >= 1.11

Model
-----
Under the null, hap1 reads ~ BetaBinomial(n, mu, rho), where
  mu  = expected hap1 fraction (0.5, or estimated to absorb reference/mapping bias)
  rho = overdispersion (0 = pure binomial; larger = more extra-binomial variance)
Var(fraction) = mu(1-mu) * (1 + (n-1)*rho) / n   -> SD plateaus at ~sqrt(mu(1-mu)*rho)
as depth grows, so high-depth loci are no longer flagged for tiny deviations.

rho (and optionally mu) is estimated by maximum likelihood across all loci within
each sample, so that samples with noisier libraries get their own null.
"""
import numpy as np
import pandas as pd
from scipy import optimize, stats
from scipy.special import expit, logit

RHO_MIN, RHO_MAX = 1e-6, 0.99

def _ab(mu, rho):
    """Convert (mean, overdispersion) to beta shape parameters."""
    s = (1.0 - rho) / rho
    return mu * s, (1.0 - mu) * s


def _two_sided_p(k, n, mu, rho):
    """Two-sided beta-binomial p-value (doubled smaller tail, capped at 1)."""
    a, b = _ab(mu, rho)
    lower = stats.betabinom.cdf(k, n, a, b)         # P(X <= k)
    upper = stats.betabinom.sf(k - 1, n, a, b)      # P(X >= k)
    return np.minimum(1.0, 2.0 * np.minimum(lower, upper))


def fit_betabinom(k, n, fit_mu=False, min_n=10, trim_p=1e-2, trim_rounds=3):
    """
    Estimate (mu, rho) by maximum likelihood.

    fit_mu      : estimate mu too (otherwise fixed at 0.5)
    min_n       : ignore loci below this depth when fitting
    trim_p      : after a first fit, drop loci with p < trim_p and refit, so that
                  genuinely imbalanced loci don't inflate rho. Set to 0 to disable.
                  Tradeoff: too weak a trim leaves imbalanced loci inflating rho (you lose
                  sensitivity); too aggressive a trim truncates the null and deflates rho
                  (you gain false positives). 0.01 with 3 rounds worked well in simulation.
    trim_rounds : number of trim/refit rounds
    """
    k, n = np.asarray(k), np.asarray(n)
    keep = n >= min_n
    mu, rho = 0.5, 0.01

    for rnd in range(trim_rounds + 1):
        kk, nn = k[keep], n[keep]
        if len(kk) < 20:
            raise ValueError(f"Only {len(kk)} loci usable for fitting; need more (>= ~50 recommended).")

        def nll(theta):
            m = expit(theta[0]) if fit_mu else 0.5
            r = np.clip(expit(theta[-1]), RHO_MIN, RHO_MAX)
            a, b = _ab(m, r)
            return -stats.betabinom.logpmf(kk, nn, a, b).sum()

        x0 = [logit(mu), logit(rho)] if fit_mu else [logit(rho)]
        res = optimize.minimize(nll, x0, method="Nelder-Mead")
        mu = expit(res.x[0]) if fit_mu else 0.5
        rho = float(np.clip(expit(res.x[-1]), RHO_MIN, RHO_MAX))

        if trim_p > 0 and rnd < trim_rounds:
            p = _two_sided_p(k, n, mu, rho)
            keep = (n >= min_n) & (p >= trim_p)

    return {"mu": float(mu), "rho": rho, "n_fit": int(keep.sum())}


def analyze(df, hap1="hap1", hap2="hap2", sample="sample",
            delta=0.10, alpha=0.05, fit_mu=False, **fit_kwargs):
    """
    Run the analysis for every sample in `df` (pass sample=None if you only have one).

    delta : minimum meaningful deviation of the hap1 fraction from mu (0.10 = outside 40-60%)
    alpha : FDR threshold (Benjamini-Hochberg, applied within each sample)

    Call rule:  q < alpha  AND  |fraction - mu| > delta
    Extra columns give a shrunken estimate, 95% credible interval and the posterior
    probability that the true fraction lies outside mu +/- delta.
    """
    df = df.copy()
    if sample is None or sample not in df.columns:
        sample = "_all"
        df[sample] = "all"

    out, params = [], []
    for s, g in df.groupby(sample, sort=False):
        k = g[hap1].to_numpy(dtype=int)
        n = k + g[hap2].to_numpy(dtype=int)
        fit = fit_betabinom(k, n, fit_mu=fit_mu, **fit_kwargs)
        mu, rho = fit["mu"], fit["rho"]
        params.append({sample: s, **fit})

        a, b = _ab(mu, rho)
        res = g.copy()
        res["n"] = n
        with np.errstate(invalid="ignore", divide="ignore"):
            res["frac_hap1"] = np.where(n > 0, k / n, np.nan)
        p = np.where(n > 0, _two_sided_p(k, n, mu, rho), np.nan)
        res["p_value"] = p

        q = np.full(len(g), np.nan)
        ok = ~np.isnan(p)
        q[ok] = stats.false_discovery_control(p[ok], method="bh")
        res["q_value"] = q

        # Empirical-Bayes posterior: prior = Beta(a, b) fitted above
        pa, pb = a + k, b + (n - k)
        res["post_mean"] = pa / (pa + pb)
        res["post_ci_lo"] = stats.beta.ppf(0.025, pa, pb)
        res["post_ci_hi"] = stats.beta.ppf(0.975, pa, pb)
        res["post_prob_beyond_delta"] = (stats.beta.cdf(mu - delta, pa, pb)
                                         + stats.beta.sf(mu + delta, pa, pb))

        res["imbalanced"] = (res["q_value"] < alpha) & ((res["frac_hap1"] - mu).abs() > delta)
        res["direction"] = np.where(res["imbalanced"],
                                    np.where(res["frac_hap1"] > mu, "hap1", "hap2"), "")
        out.append(res)

    return pd.concat(out).drop(columns=["_all"], errors="ignore"), pd.DataFrame(params)


# ---------------------------------------------------------------------------
if __name__ == "__main__":
    # Simulated demo: noisy-but-balanced loci plus 5% truly imbalanced ones.
    rng = np.random.default_rng(1)
    # Number of loci to simulate
    L = 500
    depth = np.maximum(rng.lognormal(4.5, 1.0, L).astype(int), 1)

    true_rho = 0.02
    a, b = _ab(0.5, true_rho)
    p_true = rng.beta(a, b, L)
    imb = rng.random(L) < 0.05
    p_true[imb] = rng.choice([0.25, 0.75], imb.sum())
    h1 = rng.binomial(depth, p_true)
    demo = pd.DataFrame({"sample": "S1", "locus": np.arange(L),
                         "hap1": h1, "hap2": depth - h1})
    res, params = analyze(demo)
    print("Fitted parameters (true rho = %.3f):" % true_rho)
    print(params.to_string(index=False), "\n")

    # Compare against a naive binomial test
    naive = np.array([stats.binomtest(int(x), int(x + y), 0.5).pvalue if x + y > 0 else 1
                      for x, y in zip(demo.hap1, demo.hap2)])
    hi = res["n"] > 100
    print(f"Naive binomial, p<0.05 at depth>{L} : {(naive[hi.to_numpy()] < 0.05).mean():.1%} of loci")
    print(f"Beta-binomial call at depth>{L}     : {res.loc[hi, 'imbalanced'].mean():.1%} of loci")
    print(f"(truly imbalanced fraction is {imb.mean():.1%})")
    upper = 50
    d200 = (res["n"] > upper).to_numpy()
    print(f"\nBeta-binomial call, loci with depth > {upper}:")
    print(f"  true imbalanced called : {res.loc[imb & d200, 'imbalanced'].mean():.1%}")
    print(f"  balanced called        : {res.loc[~imb & d200, 'imbalanced'].mean():.1%}")
    print("\nExample rows:")
    print(res.head(5)[["locus", "hap1", "hap2", "frac_hap1", "p_value", "q_value",
                       "post_mean", "post_prob_beyond_delta", "imbalanced"]].round(4).to_string(index=False))
    print(res['frac_hap1'].describe())
