#!/usr/bin/env python

import sys
import argparse
import itertools
from math import exp

import numpy as np
import pandas as pd

from scipy.stats import (
    chi2_contingency,
    binom,
    beta,
    binomtest,
    fisher_exact,
)

from scipy.special import betaln
from statsmodels.stats.multitest import multipletests

from smahtkit.smhtid import SMaHTid
from smahtkit.readcounts import ReadCounts


 #####################
 # Statistical Tests #
 #####################


def single_test(counts: ReadCounts):
    """
    Does only a single column have non-zero alt reads.
    Returns a boolean array and an array of column names
    (or None) for rows where exactly one sample has support.
    Note: This will set_idx('idx')
    """
    tot = counts.alt + counts.ref
    covered_counts = ~np.isnan(tot) & (tot > 0)
    present_counts = (counts.alt != 0) & covered_counts
    is_single = present_counts.sum(axis=1) == 1

    # For each row, get the index of the nonzero column (-1 if not single)
    # argmax returns the first True index; only meaningful where is_single is True
    col_indices = np.where(is_single, present_counts.argmax(axis=1), -1)

    col_names = np.where(is_single, counts.col[col_indices], None)
    ret = pd.DataFrame({"idx": counts.idx,
                        "present": col_names,
                        "n_covered": covered_counts.sum(axis=1),
                        "n_present": present_counts.sum(axis=1),
                        "reject": is_single})
    # How many possible columns there are
    ret['n_col'] = counts.alt.shape[1]

    return ret.set_index('idx')


def outlier_test(counts: ReadCounts, alpha=0.05, vaf_est='loo', **kwargs):
    """
    Run a binomial test of per-col VAF against pooled VAF to identify enrichment/depletion.

    H0: VAF is consistent across columns

    vaf_est options:
        'loo'        : Leave-one-out pooled VAF (default). For each cell (i, j), the null
                       is computed from all columns except j, making it robust to outliers.
        'cumulative' : Global pooled VAF from row-wise total alt / total depth.
    """
    k = np.array(counts.alt)
    n = np.array(counts.alt + counts.ref)

    if vaf_est == 'loo':
        # Sum across all columns, then subtract the focal column — fully vectorized
        alt_row_sum = k.sum(axis=1, keepdims=True)   # (variants, 1)
        dep_row_sum = n.sum(axis=1, keepdims=True)   # (variants, 1)

        loo_alt = alt_row_sum - k   # (variants, cols)
        loo_dep = dep_row_sum - n   # (variants, cols)

        # null_vaf is now a (variants, cols) matrix — one null per cell
        with np.errstate(invalid="ignore", divide="ignore"):
            null_vaf = np.where(loo_dep > 0, loo_alt / loo_dep, np.nan)
    elif vaf_est == 'cumulative':
        null_vaf = (k.sum(axis=1) / n.sum(axis=1))[:, np.newaxis]
    else:
        raise ValueError(f"Unknown vaf_est: {vaf_est!r}. Choose 'loo', 'cumulative', or 'median'.")

    valid = n > 0

    # Vectorized two-sided binomial p-value - null_vaf
    # null_vaf brodcast for scalar-per-row cumulative and per-cell loo shapes
    with np.errstate(invalid="ignore", divide="ignore"):
        p_low  = binom.cdf(k, n, null_vaf)
        p_high = binom.sf(k - 1, n, null_vaf)
        p_vals = np.minimum(2 * np.minimum(p_low, p_high), 1.0)
        p_vals = np.where(valid, p_vals, np.nan)
        col_vaf = np.where(valid, k / n, np.nan)

    # For the output 'null_vaf' column, report the per-cell LOO null or the shared null
    reported_null = null_vaf if vaf_est == 'loo' else np.broadcast_to(null_vaf, k.shape)

    var_idx, tis_idx = np.indices(k.shape)
    data = pd.DataFrame({
        'idx':        np.array(counts.idx)[var_idx.ravel()],
        'col':        np.array(counts.col)[tis_idx.ravel()],
        'vaf':        col_vaf.ravel().round(4),
        'null_vaf':   reported_null.ravel().round(4),
        'p_value':    p_vals.ravel().round(6),
    })

    data['reject'] = ~data['p_value'].isna() & (data['p_value'] < alpha)
    data['finished'] = data['reject']
    # Direction
    gt = data['vaf'] > data['null_vaf']
    valid_mask = ~data['p_value'].isna()
    data['direction'] = np.select(
        [~valid_mask, data['reject'] & gt, data['reject'] & ~gt],
        ['?', '+', '-'],
        default='/'
    )

    return data


def clopper_pearson(k, n, lower=True, alpha=0.05):
    """
    Upper or bound of a Clopper-Pearson CI.
    k = alt reads, n = total reads, alpha = significance level.
    Returns the upper confidence bound on the true VAF.
    """
    if n == 0:
        return 1.0
    if lower:
        return beta.ppf(alpha / 2, k, n - k + 1)
    # Upper bound is the (1 - alpha/2) quantile of Beta(k+1, n-k)
    return beta.ppf(1 - alpha / 2, k + 1, n - k)


def clopper_pearson_vec(k, n, lower=True, alpha=0.05):
    """
    Vectorized Clopper-Pearson CI bound over arrays k, n.
    """
    k = np.asarray(k, dtype=float)
    n = np.asarray(n, dtype=float)

    if lower:
        result = beta.ppf(alpha / 2, k, n - k + 1)
    else:
        result = beta.ppf(1 - alpha / 2, k + 1, n - k)

    return np.where(n == 0, 1.0, result)


def exclusive_test(counts: ReadCounts, alpha=0.05, vaf_est=None, mc=False, **kwargs):
    """
    For SVs with alt support in exactly one sample, test whether absence
    in all other samples is explainable by coverage alone.

    H0: The variant exists at the observed VAF everywhere; we just didn't
        see it in other samples due to low coverage.

    A low p-value means the absences are unlikely under H0 → the variant
    is likely tissue-specific.

    The VAF parameter is estimated as ALT / COV by default.
    Use vaf_est 'upper' or 'lower' for an upper/lower bound estimation of VAF using Clopper-Pearson. Note that the
    returned vaf will be this estimation, not the observed.

    Returns a pandas DataFrame with columns:
        idx: the ReadCounts.idx row that was tested
        col: The significant ReadCounts.col (if any)
        vaf: The variant allele fraction of the column
        p_joint: the joint probabiliy
        adj_pjoint: fbh adjusted joint probability
        reject: boolean by reject

    Note: Every `ReadCounts.idx` row has a test result returned. However, only rows with exactly one column having
    alternate read support are tested. Untested rows will have `result['col'].isna()` and `result['reject'] == False`
    whereas tested rows will record the singly covered col and either True/False in the reject

    For SVs with alt support in a single sample, we computed the joint probability 
    of observing zero alt reads in all remaining samples under the null hypothesis 
    that the variant is present at the Clopper-Pearson lower/upper bound VAF. Low joint 
    probability indicates that coverage alone cannot explain the absence.
    """
    N, M = counts.alt.shape
    alt = counts.alt.astype(float)
    total = (counts.alt + counts.ref).astype(float)

    # --- valid mask: non-NaN and total > 0  [N, M] ---
    valid = ~(np.isnan(alt) | np.isnan(total)) & (total > 0)

    # zero out invalid cells so sums ignore them
    alt_v = np.where(valid, alt,   0.0)
    total_v = np.where(valid, total, 0.0).astype(int)

    # --- "exactly one positive sample" filter ---
    positive = valid & (alt_v > 0)          # [N, M] bool
    n_positive = positive.sum(axis=1)       # [N]
    single_mask = (n_positive == 1)         # [N] rows to test

    # For each row: index of the one positive column
    # argmax on a bool array returns the first True index
    # [N]  (only meaningful where single_mask)
    pos_j = np.argmax(positive, axis=1)

    # Gather positive-sample alt and total using fancy indexing
    row_idx = np.arange(N)
    pos_alt = alt_v[row_idx, pos_j]             # [N]
    pos_tot = total_v[row_idx, pos_j]           # [N]

    # --- VAF computation ---
    if vaf_est is None:
        with np.errstate(invalid="ignore", divide="ignore"):
            vaf = np.where(pos_tot > 0, pos_alt / pos_tot, np.nan)
    else:
        use_lower = (vaf_est == "lower")
        vaf = clopper_pearson_vec(pos_alt, pos_tot, use_lower, alpha)
        vaf = np.where(single_mask, vaf, np.nan)  # blank out irrelevant rows

    # --- joint miss probability ---
    # P(alt==0 | n, vaf) = (1-vaf)^n  →  joint = prod over other samples = (1-vaf)^sum(other_n)
    # sum(other_n) = sum(total_v over valid cells) - pos_tot
    # [N]  sum of all valid totals
    total_valid_sum = total_v.sum(axis=1)
    # [N]  exclude the positive sample
    other_total_sum = total_valid_sum - pos_tot

    with np.errstate(invalid="ignore"):
        p_joint = np.where(
            single_mask,
            (1.0 - vaf) ** other_total_sum,
            np.nan
        )                                                        # [N]

    # --- build DataFrame directly from arrays ---
    idx_arr = np.array(counts.idx)
    col_arr = np.array(counts.col)

    # The positive sample's column name (only meaningful where single_mask)
    pos_sample = col_arr[pos_j]                                  # [N]

    result = pd.DataFrame({
        "idx":     idx_arr,
        "col":     np.where(single_mask, pos_sample, None),
        "vaf":     np.where(single_mask, np.round(vaf,     4), np.nan),
        "p_value": np.where(single_mask, np.round(p_joint, 6), np.nan),
    })

    result["reject"] = False
    valid_p = ~result["p_value"].isna()
    if not valid_p.any():
        return result

    # --- multiple testing correction ---
    if mc:
        result["adj_p_value"] = np.nan
        reject, adj_p, _, _ = multipletests(
            result.loc[valid_p, "p_value"],
            alpha=alpha,
            method="fdr_bh"
        )
        result.loc[valid_p, "adj_p_value"] = adj_p.round(6)
        result.loc[valid_p, "reject"] = reject
    else:
        result.loc[valid_p, "reject"] = result.loc[valid_p, 'p_value'] < alpha

    return result


def batch_effect_p0(alt_A, ref_A, depth_B):
    """
    Probability of seeing zero alt reads in platform B given observations
    from platform A, using a Beta-Binomial model with a uniform prior.

    P(X_B = 0) = B(alpha, beta + depth_B) / B(alpha, beta)

    Parameters
    ----------
    alt_A   : alt read count from platform A
    ref_A   : ref read count from platform A
    depth_B : total depth in platform B (ref_B when alt_B == 0)

    Returns
    -------
    float : p-value; small values indicate a likely batch effect
    """
    alpha = alt_A + 1
    beta  = ref_A + 1
    log_p0 = betaln(alpha, beta + depth_B) - betaln(alpha, beta)
    return exp(log_p0)

def batch_test_full(reads, result):
    """
    Same as batch_test, but instead of only analyzing coverage over the subset of col that we're analyzing, it checks
    over all
    # Parameters
    # reads - original
    # result - The test result

    returns the same result, but with extra columns added about the batch effect tests
    """
    if not result['finished'].any():
        result['batch_p_value'] = np.nan
        return result

    finished_idx = result.loc[result['finished'], 'idx']
    sub = reads[np.isin(reads.idx, finished_idx)]
    by_plat = sub.collapse_by([_.platform_name for _ in reads.col])
    frame = by_plat.to_frame()
    platform_counts = (
            frame.set_index(['idx', 'col'])
            .unstack(level='col')
            .swaplevel(axis=1)
            .sort_index(axis=1)
    )
    platform_counts.columns = [f"{p}_{s}" for p, s in platform_counts.columns]

    n_groups = len(by_plat.idx)

    idx_vals = np.empty(n_groups, dtype=object)
    p_vals   = np.full(n_groups, np.nan)

    for i, (idx, m_df) in enumerate(frame.groupby('idx')):
        idx_vals[i] = idx
        if len(m_df) != 2:
            continue

        m_df = m_df.sort_values(by='alt', ascending=False)
        # No coverage from second plat
        if m_df.iloc[1][['ref', 'alt']].sum() == 0:
            continue
        if m_df.iloc[1]['alt'] == 0: # Only found in one
            p_vals[i] = batch_effect_p0(m_df.iloc[0]['alt'], m_df.iloc[0]['ref'], m_df.iloc[1]['ref'])
        else:
            box = m_df[['ref', 'alt']].values
            _, p_vals[i] = fisher_exact(box)

    p_frame = pd.DataFrame({
        'idx': idx_vals,
        'batch_p_value': np.round(p_vals, 6)
    }).set_index(['idx'])

    # Join p-values and platform breakdown together before merging onto result
    p_frame = p_frame.join(platform_counts)

    result = result.set_index(['idx']).join(p_frame, how='left').reset_index()
    return result

def batch_test(reads, result, level):
    """
    # Parameters
    # reads - original
    # result - The test result
    # level - collapse_level used
        If None, we're comparing all the col, not just a subset

    returns the same result, but with extra columns added about the batch effect tests
    """
    if not result['finished'].any():
        result['batch_p_value'] = np.nan
        return result

    # 1 - Filter reads to finished results only
    finished_idx = result.loc[result['finished'], 'idx']
    sub = reads[np.isin(reads.idx, finished_idx)]

    # 2 - Build new_level with vectorized string join instead of map+lambda
    platforms = [_.platform_name for _ in reads.col]
    if level is not None:
        new_level = [f"{l}:::{p}" for l, p in zip(level, platforms)]
    else:
        new_level = [f"all:::{p}" for p in platforms]

    # 3 - Collapse by new level
    new_collapse = sub.collapse_by(new_level)

    # 4 - Use str.split with expand=True instead of apply(pd.Series) -- much faster
    countable = new_collapse.to_frame()
    countable[['col', 'platform']] = countable['col'].str.split(':::', expand=True)

    # 5 - Use merge instead of MultiIndex isin -- avoids index construction overhead
    result_keys = result[['idx', 'col']]
    care_about = countable.merge(result_keys, on=['idx', 'col'], how='inner')

    # 6 - Collapse read counts
    almost_box = care_about.groupby(['idx', 'col', 'platform'])[['ref', 'alt']].sum()

    # 7 - Vectorize what we can; only loop for the statistical tests themselves
    #     Pre-compute group structure to minimize per-iteration overhead
    grouped = almost_box.groupby(level=[0, 1])
    n_groups = grouped.ngroups

    if n_groups == 0:
        result['batch_p_value'] = np.nan
        return result

    # Pivot platform ref/alt counts into flat columns (replaces d.stack() in loop)
    # Result: idx, col as index; columns like "platform1_ref", "platform1_alt", ...
    platform_counts = (
        almost_box
        .unstack(level='platform')        # platforms become column MultiIndex
        .swaplevel(axis=1)                # (platform, ref/alt) -> (ref/alt, platform)... 
        .sort_index(axis=1)               # ...actually we want (platform, stat) order
    )
    # Flatten MultiIndex columns to "platform_ref", "platform_alt"
    platform_counts.columns = [f"{p}_{s}" for p, s in platform_counts.columns]

    idx_vals = np.empty(n_groups, dtype=object)
    col_vals = np.empty(n_groups, dtype=object)
    p_vals   = np.full(n_groups, np.nan)

    for i, (g, d) in enumerate(grouped):
        idx_vals[i], col_vals[i] = g
        box = d[['ref', 'alt']].values
        n = box.shape[0]
        # All alts are zero or no coverage from one platform
        if np.all(box[:, 1] == 0) or (box.sum(axis=1) == 0).any():
            # Here's where we'd do the single tech test.. No, alt only dang.
            continue

        if n == 2:
            mask = box[:, 1] == 0
            # Use betabinom whenever we have zero alt support from one of the platforms
            if mask.any():
                present = box[~mask][0]
                absent = box[mask][0]
                p_vals[i] = batch_effect_p0(present[1], present[0], absent[0])
            else:
                _, p_vals[i] = fisher_exact(box)
        elif n > 2: # This shouldn't ever happen, only have HiFi and ONT
            _, p_vals[i] = chi2_contingency(box)

    p_frame = pd.DataFrame({
        'idx': idx_vals,
        'col': col_vals,
        'batch_p_value': np.round(p_vals, 6)
    }).set_index(['idx', 'col'])

    # Join p-values and platform breakdown together before merging onto result
    p_frame = p_frame.join(platform_counts)

    result = result.set_index(['idx', 'col']).join(p_frame).reset_index()
    return result

#################
# Cascade Tests #
#################


def core_specific(reads, core_level, tissue_level, **kwargs):
    """
    Perform the Core specific test
    return the annotated DataFrame of passing SVs
    """
    by_core = reads.collapse_by(core_level)

    candidates = single_test(by_core)
    by_core_subset = by_core[candidates['reject']]
    by_core_subset = by_core_subset.min_mask(**kwargs)
    candidates.drop(columns=['reject', 'present'], inplace=True)

    results = []
    for tissue in set(tissue_level):
        mask = np.array(
            [_.split('_')[0] == tissue for _ in by_core_subset.col])

        # Only perform core specific test when something to compare against
        if np.sum(mask) < 2:
            continue
        subset = by_core_subset[:, mask]

        # Only look at the tissues that still have alt coverage
        still_has = single_test(subset)['reject']
        if not still_has.any():
            continue
        subset = subset[still_has]

        result = exclusive_test(subset, **kwargs)
        results.append(result)

    # Empty
    if not results:
        return None
    result = pd.concat(results).set_index('idx').join(candidates).reset_index()
    # We're finished checking this SV because it only has read support in a single core, 
    # so by definition it will only have read support in a single tissue/layer
    result['finished'] = True
    result['annotation'] = "Core Specific"
    result = batch_test(reads, result, core_level)
    return result


def anno_specific(reads, level, blood_idxs=None, fibro_idxs=None, key='Tissue', **kwargs):
    """
    This is generalizable to tissue, tissue*, and layer specific tests, just provide a key for the annotation
    blood_idxs and fibro_idxs are lists of the idxs that have any blood/fibroblast presence

    Parameters to ReadCounts.min_mask and outlier_test can be passed via kwargs
    """
    collapsed = reads.collapse_by(level)
    collapsed = collapsed.min_mask(**kwargs)
    candidates = single_test(collapsed).drop(columns=['reject', 'present'])
    result = exclusive_test(collapsed, **kwargs)
    result['finished'] = result['col'].notna()
    result = result.set_index('idx').join(candidates).reset_index()
    result['has_blood'] = result['idx'].map(blood_idxs) if blood_idxs is not None else None
    result['has_fibro'] = result['idx'].map(fibro_idxs) if fibro_idxs is not None else None
    result['annotation'] = f"{key} Specific"
    result = batch_test(reads, result, level)

    return result


def anno_outlier(reads, level, blood_idxs=None, fibro_idxs=None, key="Tissue", **kwargs):
    """
    This is generalizable to tissue, tissue*, and layer outlier tests, just provide a key for the annotation
    blood_idxs and fibro_idxs are lists of the idxs that have any blood/fibroblast presence

    Parameters to ReadCounts.min_mask and outlier_test can be passed via kwargs
    """
    collapsed = reads.collapse_by(level)
    collapsed = collapsed.min_mask(**kwargs)
    candidates = single_test(collapsed).drop(columns=['reject', 'present'])

    result = outlier_test(collapsed, **kwargs)
    result['keep'] = result['reject']
    result = result.set_index('idx').join(candidates).reset_index()
    result['has_blood'] = result['idx'].map(
        blood_idxs) if blood_idxs is not None else None
    result['has_fibro'] = result['idx'].map(
        fibro_idxs) if fibro_idxs is not None else None
    result['annotation'] = f"{key} Outlier"
    result = batch_test_full(reads, result)

    return result

def calculate_present_present(reads: ReadCounts):
    """
    Simple calculation of what percent of tissues have a call
    """
    frame = single_test(reads)
    frame['present_tissue_pct'] = ((frame['n_present'] / frame['n_covered']) * 100).astype(np.int32)
    return frame[['present_tissue_pct']]

def annotate_origin(orig_reads: ReadCounts, alpha=0.05, min_coverage=10):
    """
    Run all tests and return a DataFrame
    """
    all_results = []

    # Setup the levels
    core_level = [_.proto_core for _ in orig_reads.col]
    tissue_level = [_.protocol for _ in orig_reads.col]
    tissue_dedup_level = [_.dedup_tissue_abv for _ in orig_reads.col]
    
    # Keep this annotation for later
    present_present = calculate_present_present(orig_reads.collapse_by(tissue_level))

    # 0 - initial filtering
    tmp_reads = orig_reads.min_mask(min_coverage=min_coverage, min_col=2)
    filt_svs = np.isnan(tmp_reads.alt).all(axis=1)
    filt_cnt = filt_svs.sum()
    result0 = pd.Series(tmp_reads.idx[filt_svs], name='idx').to_frame()
    # This is counting how many after filtering, so it'll always be zero
    # I'm keeping for notes and to make sure that the e.g. n_col is an int in the output
    candidates = single_test(tmp_reads[filt_svs])
    result0 = result0.set_index('idx').join(candidates).reset_index()
    result0['finished'] = True
    result0['annotation'] = "Low Support"
    result0['level'] = 0
    all_results.append(result0)
    # Initial filtering is performed on original reads so e.g. fibro counts are post filtering
    orig_reads = orig_reads[np.isin(orig_reads.idx, tmp_reads.idx[~filt_svs])]
    print(f"Level 0 - Low Support:\t{filt_cnt}", file=sys.stderr)

    # We're going to be marking these columns
    # We no longer want to check anything that we've already annotated
    blood_col_mask = np.array([_.tissue_abv == 'BLOO' for _ in orig_reads.col])
    blood = orig_reads[:, blood_col_mask]
    blood_idxs = blood.presence_lookup()

    fibro_col_mask = np.array([_.tissue_abv == 'FBRO' for _ in orig_reads.col])
    fibro = orig_reads[:, fibro_col_mask]
    fibro_idxs = fibro.presence_lookup()

    # 1 - Core specific
    result1 = core_specific(orig_reads, core_level,
                            tissue_level, vaf_est='lower',
                            min_coverage=min_coverage)
    if result1 is not None:
        print(f"Level 1 - Core Specific:\t{result1['finished'].sum()}", file=sys.stderr)
        result1['level'] = 1
        all_results.append(result1)
        reads = orig_reads.exclude_finished(result1)
    else:
        print(f"Level 1 - Core Specific:\t0", file=sys.stderr)
        reads = orig_reads

    # 2 - Tissue Specific
    result2 = anno_specific(reads, tissue_level, blood_idxs, fibro_idxs,
                            vaf_est='lower', min_coverage=min_coverage)
    result2['level'] = 2
    print(f"Level 2 - Tissue Specific:\t{result2['finished'].sum()}", file=sys.stderr)
    all_results.append(result2)
    reads = reads.exclude_finished(result2)

    # 3 - Tissue* Specific
    result3 = anno_specific(reads, tissue_dedup_level, blood_idxs, fibro_idxs,
                            key="Tissue*", vaf_est='lower',
                            min_coverage=min_coverage)
    result3['level'] = 3
    print(f"Level 3 - Tissue* Specific:\t{result3['finished'].sum()}", file=sys.stderr)
    all_results.append(result3)
    reads = reads.exclude_finished(result3)

    # 4 - Blood Test
    by_tissue = reads.collapse_by(tissue_level)
    candidates = single_test(by_tissue).drop(columns=['reject', 'present'])

    blood = outlier_test(by_tissue, min_coverage=min_coverage)
    blood['finished'] = (blood['col'] == '3A') & (blood['direction'] == '+')
    blood = blood.set_index('idx').join(candidates).reset_index()
    blood['has_blood'] = True
    blood['has_fibro'] = blood['idx'].isin(fibro_idxs)
    blood['annotation'] = "Blood"
    blood['level'] = 4
    blood = batch_test_full(reads, blood)
    print(f"Level 4 - Blood:\t{blood['finished'].sum()}", file=sys.stderr)
    all_results.append(blood)

    # Now we want to drop blood and fibroblast
    reads = reads.exclude_finished(blood)[:, ~(blood_col_mask | fibro_col_mask)]
    
    # 5 - Repeat #2 and #3
    tissue_level = [_.protocol for _ in reads.col]
    tissue_dedup_level = [_.dedup_tissue_abv for _ in reads.col]

    # Repeat 2 - Tissue Specific
    result5_2 = anno_specific(reads, tissue_level, blood_idxs, fibro_idxs,
                              vaf_est='lower', min_coverage=min_coverage)
    result5_2['level'] = 5.2
    print(f"Level 5.2 - Tissue Specific:\t{result5_2['finished'].sum()}", file=sys.stderr)
    all_results.append(result5_2)
    reads = reads.exclude_finished(result5_2)

    # Repeat 3 - Tissue* Specific
    result5_3 = anno_specific(reads, tissue_dedup_level, blood_idxs, fibro_idxs,
                              key="Tissue*", vaf_est='lower', min_coverage=min_coverage)
    result5_3['level'] = 5.3
    print(f"Level 5.3 - Tissue* Specific:\t{result5_3['finished'].sum()}", file=sys.stderr)
    all_results.append(result5_3)
    reads = reads.exclude_finished(result5_3)

    # 6 - Layer Specific
    layer_level = [_.layer for _ in reads.col]
    result6 = anno_specific(reads, layer_level, blood_idxs, fibro_idxs,
                            min_col=3, min_coverage=min_coverage, key="Layer")
    result6['level'] = 6
    print(f"Level 6 - Layer Specific:\t{result6['finished'].sum()}", file=sys.stderr)
    all_results.append(result6)
    reads = reads.exclude_finished(result6)

    # 7 - Tissue Outlier
    result7 = anno_outlier(reads, tissue_level, blood_idxs, fibro_idxs,
                           min_col=3, min_coverage=min_coverage)
    result7['level'] = 7
    print(f"Level 7 - Tissue Outlier:\t{result7['finished'].sum()}", file=sys.stderr)
    all_results.append(result7)
    reads = reads.exclude_finished(result7)

    # 8 - Tissue* Outlier
    result8 = anno_outlier(reads, tissue_dedup_level, blood_idxs, fibro_idxs,
                           min_col=3, key="Tissue*", min_coverage=min_coverage)
    result8['level'] = 8
    print(f"Level 8 - Tissue* Outlier:\t{result8['finished'].sum()}", file=sys.stderr)
    all_results.append(result8)
    reads = reads.exclude_finished(result8)

    # 9 - Layer Outlier
    result9 = anno_outlier(reads, layer_level, blood_idxs, fibro_idxs,
                           min_col=3, vaf_est='cumulative', key="Layer",
                           min_coverage=min_coverage)
    result9['level'] = 9
    print(f"Level 9 - Layer Outlier:\t{result9['finished'].sum()}", file=sys.stderr)
    all_results.append(result9)
    reads = reads.exclude_finished(result9)

    # 10 - Ubiquitous
    if len(reads.idx) != 0:
        result10 = pd.Series(reads.idx, name='idx').to_frame()
        m_reads = orig_reads[np.isin(orig_reads.idx, reads.idx)]
        candidates = single_test(m_reads).drop(columns=['reject', 'present'])
        result10 = result10.set_index('idx').join(candidates).reset_index()
        alt = np.where(np.isnan(m_reads.alt), 0, m_reads.alt)
        ref = np.where(np.isnan(m_reads.ref), 0, m_reads.ref)
        total = alt + ref
        result10['null_vaf'] = (alt.sum(axis=1) / total.sum(axis=1)).round(4)
        result10['has_blood'] = result10['idx'].map(blood_idxs)
        result10['has_fibro'] = result10['idx'].map(fibro_idxs)
        result10['annotation'] = "Unclassified"
        result10['level'] = 10
        result10['finished'] = True
        # Gather read support by platform -- I might not need to do this any more with the changed batch effect?
        #tmp = reads.dbg_view()
        #block = tmp.groupby(['idx', 'platform_name'])[['ref', 'alt']].sum()
        #row = block.stack()
        #row.index = pd.MultiIndex.from_tuples(
            #[(x[0], f'{x[1]}_{x[2]}') for x in row.index]
        #)
        #row = row.unstack(level=[1])
        #result10 = result10.set_index('idx').join(row).reset_index()
        result10 = batch_test_full(m_reads, result10)
        print(f"Level 10 - Unclassified:\t{result10['finished'].sum()}", file=sys.stderr)
        all_results.append(result10)

    # Organize the results' columns
    column_order = ['idx', 'col', 'annotation', 'level',
                    'n_col', 'n_covered', 'n_present']
    uniq_plat = sorted(list(set([_.platform_name for _ in reads.col])))
    column_order += list(map(lambda x: f'{x[0]}_{x[1]}',
                         itertools.product(uniq_plat, ['ref', 'alt'])))
    column_order += ['batch_p_value', 'vaf', 'null_vaf',
                     'p_value', 'direction', 'reject',
                     'has_blood', 'has_fibro', 'finished']

    result = pd.concat(all_results, ignore_index=True)
    result = result[[c for c in column_order if c in result.columns]]

    # Add in a number of platforms column
    plat_alt = [_ for _ in result.columns if _.endswith('_alt')]
    result['n_platform'] = (result[plat_alt] > 0).sum(axis=1)

    # Multiple testing correction
    mask_core = result['level'].isin([1])
    mask_joint = result['level'].isin([2, 3, 5.2, 5.3, 6])
    mask_binom = result['level'].isin([4, 7, 8, 9])
    has_test = result['p_value'].notna()

    for mask in [mask_core, mask_joint, mask_binom]:
        mask = mask & has_test
        if not mask.any():
            continue
        idx = result.index[mask]
        _, q, _, _ = multipletests(result.loc[idx, 'p_value'], method='fdr_bh')
        result.loc[idx, 'adj_pval'] = q.round(6)
        result.loc[idx, 'adj_reject'] = q < alpha
    
    # Only keep the 'finished' annotations
    result = result[result['finished']].copy().drop(columns='finished')

    # Tack on the present present column
    result = result.set_index('idx').join(present_present, how='left').reset_index()
    print(f"Produced {len(result)} annotations")

    return result

######
# UI #
######

def batch_origin_main(args):
    """
    SMaHT SV Origin Annotator
    
    Provide the `--parquet` from SV pipe or matrices of `--ref` `--alt` read counts as input.

    For `--ref/--alt`, the first column is expected to be the variants' identifiers.
    Columns 1+ SMaHT ids of production samples that are compliant with v2.1 of the SMaHT nomenclature documentation. 
    Uncovered (a.k.a. missing) sites should have blank cells for nulls.
    """
    parser = argparse.ArgumentParser(prog="batch_origin", description=origin_cmd.__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-p", "--parquet", type=str, default=None,
                        help="Input Parquet SV File")
    group.add_argument("-r", "--ref", type=str, default=None,
                        help="Input tsv with reference read matrix")
    parser.add_argument("-a", "--alt", type=str, default=None,
                        help="Input tsv with alternate read matrix")
    parser.add_argument("-o", "--output", metavar="OUT", type=str, default="/dev/stdout",
                        help="Output File (stdout)")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Alpha for reject column (%(default)s)")
    parser.add_argument("--min-coverage", default=10, type=int,
                        help="Minimum coverage of a column to participate in tests (%(default)s)")

    args = parser.parse_args(args)

    if args.ref is not None and args.alt is None:
        parser.error("--ref requires --alt")
    if args.alt is not None and args.ref is None:
        parser.error("--alt requires --ref")
    if args.parquet is not None and (args.ref or args.alt):
        parser.error("--parquet cannot be used with --ref/--alt")

    if args.parquet:
        reads = ReadCounts.from_pq(args.parquet)
    elif args.ref and args.ref:
        reads = ReadCounts.from_tsv(args.ref, args.alt)
    else:
        parser.error("Missing Inputs")

    print(f"Annotating {len(reads.idx)} variants across {len(reads.col)} samples")
    results = annotate_origin(reads, args.alpha, args.min_coverage)

    results.to_csv(args.output, sep='\t', index=False)

def origin_main(args):
    """
    Origin command for read delta tsvs with columns sample (SMaHTid), hap, and is_germ
    Also, I should maybe go back to the qdpi and extractor step and build the documentation from that
    It would be nice to have configs, but I think for now we can assume SMaHT
    """
    parser = argparse.ArgumentParser(prog="origin", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_tsv", type=str,
                        help="Input anno_reads.tsv")
    parser.add_argument("-o", "--output", type=str, default="/dev/stdout",
                        help="Output origin.tsv annotations (%(default)s)")
    parser.add_argument("-s", "--square-reads", default=None,
                        help="Save the per-sample ref/alt read counts to file (off)")
    parser.add_argument("--alpha", type=float, default=0.05,
                        help="Alpha for reject column (%(default)s)")
    parser.add_argument("-m", "--min-sup", type=int, default=3,
                        help="Min number of reads in a sample to allow presence (%(default)s)")
    parser.add_argument("--min-coverage", type=int, default=10,
                        help="Minimum coverage of a column to participate in tests (%(default)s)")
    args = parser.parse_args(args)
    reads = pd.read_csv(args.in_tsv, sep='\t')

    read_parts = []
    parts = []
    for (chrom, start, end, donor, hap), sub_reads in reads.groupby(["chrom", "start", "end", "donor", "hap"]):
        donor_key = f'{donor}.{hap}'
        all_cols = sub_reads['full_name'].unique()

        ref = (sub_reads[sub_reads['is_germ']]
               .groupby(['full_name']).size()
               .reindex(all_cols, fill_value=0))
        ref.name = donor_key
        ref = ref.to_frame().T

        alt = (sub_reads[~sub_reads['is_germ']]
               .groupby(['full_name']).size()
               .reindex(all_cols, fill_value=0))
        # MASKING - fewer than min_sup isn't considered present
        alt[alt < args.min_sup] = 0
        alt.name = donor_key
        alt = alt.to_frame().T
        read_counts = ReadCounts(ref.values, alt.values,
                                 index=ref.index,
                                 columns=[SMaHTid(_) for _ in ref.columns])
        
        if alt.iloc[0].sum() != 0:
            origin = annotate_origin(read_counts, args.alpha, args.min_coverage)
            parts.append(origin)
        
        ref = ref.T[donor_key]
        ref.name = 'ref'

        alt = alt.T[donor_key]
        alt.name = 'alt'
        read_parts.append(pd.concat([ref, alt], axis=1).reset_index())

    output = pd.concat(parts)
    output.to_csv(args.output, index=False, sep='\t')

    if args.square_reads is not None:
        output = pd.concat(read_parts)
        output.to_csv(args.square_reads, sep='\t', index=False)
