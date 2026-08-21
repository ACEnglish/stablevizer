"""
Perform enrichment test for subset of cels in a SOM
"""
import sys
import pickle
import argparse

import pandas as pd
import numpy as np

from scipy.stats import hypergeom
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components
from scipy.spatial import cKDTree
from statsmodels.stats.multitest import multipletests

def square(data, shape):
    ret = np.zeros(shape, dtype=int)
    l = data.groupby(['X', 'Y']).size()
    for pos, v in l.items():
        ret[pos] = v
    return ret


def chi_stat(obs, exp):
    mask = exp > 0
    return np.sum((obs[mask] - exp[mask])**2 / exp[mask])


def flat_idx(x, y):
    return x * som_y + y   # matches ravel()'s C-order used everywhere above


def cluster_masses(A, cand_mask, stat):
    """Connected components of `cand_mask` under adjacency A.
    Returns [(mass, member_unit_indices), ...] with mass = sum(stat) over members."""
    idx = np.where(cand_mask)[0]
    if len(idx) == 0:
        return []
    sub = A[idx][:, idx]
    n_comp, labels = connected_components(sub, directed=False)
    return [(stat[idx[labels == c]].sum(), idx[labels == c]) for c in range(n_comp)]


def som_test(all_count, sub_count, som, ALPHA=0.05, R=2000, CLUSTER_FORMING_ALPHA=0.05, R2=10000, DIST_PRUNE_PERCENTILE=None):
    """
    Performs three tests for enrichment on a som
    alpha - significance
    R - permutation count

    CLUSTER_FORMING_ALPHA = 0.05   # liberal, *uncorrected* per-unit threshold used only
                                    # to decide which units are eligible to join a region;
                                    # actual significance comes from the permutation null
                                    # of the max region mass, below.
    R2 = 10000                     # permutations for tier2 (each iteration does more work
                                    # than tier0/1's, so this defaults smaller; raise it if
                                    # your SOM is small enough to afford it)
    DIST_PRUNE_PERCENTILE = None   # e.g. 95 to refuse to merge two grid-adjacent units
                                    # whose codebook vectors differ more than the 95th
                                    # percentile of all neighbor-pair distances (a "cliff"
                                    # in the U-matrix sense). None = pure grid adjacency,
                                    # the more standard/conservative default.
    """
    ###############################################################################
    # Global omnibus - is there enrichment anywhere?
    ###############################################################################

    shape = som.distance_map().shape
    n = all_count.ravel()
    m = sub_count.ravel()

    rng = np.random.default_rng(0)
    N, M = n.sum(), m.sum()

    expected = M * n / N

    obs_stat = chi_stat(m, expected)

    perm_counts = rng.multivariate_hypergeometric(n, M, size=R)
    perm_stats = np.array([chi_stat(pc, expected) for pc in perm_counts])
    p_global = (1 + (perm_stats >= obs_stat).sum()) / (R + 1)

    print("# Global omnibus - is there enrichment anywhere? : p =", p_global)


    ###############################################################################
    # Per-unit enrichment - which contexts are significant?
    ###############################################################################
    K = len(n)
    pvals = np.array([
        hypergeom.sf(m[k] - 1, N, n[k], M) if n[k] > 0 else np.nan
        for k in range(K)
    ])

    fold_enrichment = (m / M) / (n / N)

    valid = ~np.isnan(pvals)
    reject_compact, qval_compact, _, _ = multipletests(pvals[valid], method='fdr_bh', alpha=ALPHA)

    reject = np.zeros(K, dtype=bool)
    qval = np.full(K, np.nan)
    reject[valid] = reject_compact
    qval[valid] = qval_compact

    x_all, y_all = np.unravel_index(np.arange(K), shape)
    # I don't think I need this

    print("# Per-unit enrichment - which contexts are significant?")
    sig = reject.sum()
    present = (m != 0).sum()

    print(f"# {present} populated clusters with {sig} significant (p < {ALPHA}) Percent: {(sig / present) * 100:.2f}%")

    ###############################################################################
    # Tier 2: topology-aware cluster enrichment
    #
    # Per-unit tests treat each SOM unit independently, so real enrichment spread
    # across a few adjacent (biologically similar) contexts can fail to reach
    # per-unit significance in any single one of them even though the
    # *neighborhood* is clearly non-random. This pools evidence across contiguous
    # suprathreshold units using the SOM's grid topology (and, optionally, the
    # trained codebook distances between neighbors, to avoid merging across a
    # sharp U-matrix boundary), then tests each resulting region against a
    # max-statistic permutation null (Nichols & Holmes 2002-style cluster-mass
    # test), which controls the family-wise error rate across regions.
    ###############################################################################

    # ---- 1. Build the SOM adjacency graph --------------------------------------
    weights = som.get_weights()             # (som_x, som_y, input_len)
    xx, yy = som.get_euclidean_coordinates() # (som_x, som_y) each
    som_x, som_y, _ = weights.shape
    assert (som_x, som_y) == shape, "SOM weight grid and distance_map shape disagree"

    pts = np.column_stack([xx.ravel(), yy.ravel()])   # physical position of unit k
    tree = cKDTree(pts)
    nn_dist = tree.query(pts, k=2)[0][:, 1]           # distance to each unit's closest neighbor
    unit_spacing = np.median(nn_dist)                 # ~1.0 for both hex and rectangular grids
    edges = tree.query_pairs(r=unit_spacing * 1.05, output_type='ndarray')

    ei, ej = edges[:, 0], edges[:, 1]
    w_flat = weights.reshape(K, -1)
    edist = np.linalg.norm(w_flat[ei] - w_flat[ej], axis=1)

    dist_thresh = (np.percentile(edist, DIST_PRUNE_PERCENTILE)
                   if DIST_PRUNE_PERCENTILE is not None else np.inf)
    keep = edist <= dist_thresh

    A = csr_matrix((np.ones(keep.sum(), dtype=bool), (ei[keep], ej[keep])), shape=(K, K))
    A = A.maximum(A.T)  # symmetrize

    print(f"# Tier2 adjacency graph: {keep.sum()}/{len(edist)} SOM neighbor-edges kept"
          + ('' if DIST_PRUNE_PERCENTILE is None else f' (pruned at {DIST_PRUNE_PERCENTILE}th pct dist)'))


    # ---- 2. Observed clusters (reusing tier1's exact hypergeometric p-values) --
    with np.errstate(divide='ignore', invalid='ignore'):
        logp_obs = np.nan_to_num(-np.log10(np.clip(pvals, 1e-300, 1.0)))
    cand_obs = (pvals < CLUSTER_FORMING_ALPHA) & (n > 0)
    obs_clusters = cluster_masses(A, cand_obs, logp_obs)

    # ---- 3. Null distribution of the max cluster mass --------------------------
    #         (same multivariate-hypergeometric null as tier0/tier1)
    perm_counts2 = rng.multivariate_hypergeometric(n, M, size=R2)
    with np.errstate(divide='ignore', invalid='ignore'):
        p_perm = hypergeom.sf(perm_counts2 - 1, N, n[None, :], M)      # (R2, K), vectorized
        logp_perm = np.nan_to_num(-np.log10(np.clip(p_perm, 1e-300, 1.0)))
    cand_perm = (p_perm < CLUSTER_FORMING_ALPHA) & (n[None, :] > 0)

    null_max = np.empty(R2)
    for r in range(R2):
        comps = cluster_masses(A, cand_perm[r], logp_perm[r])
        null_max[r] = max((c[0] for c in comps), default=0.0)

    # ---- 4. Cluster-level (FWER-corrected) significance -------------------------
    cluster_rows = []
    for cid, (mass, members) in enumerate(sorted(obs_clusters, key=lambda c: -c[0])):
        p_clust = (1 + (null_max >= mass).sum()) / (R2 + 1)
        for u in members:
            ux, uy = np.unravel_index(u, shape)
            cluster_rows.append({
                'cluster_id': cid, 'X': ux, 'Y': uy, 'count': m[u],
                'unit_pval': pvals[u], 'unit_qval': qval[u], 'unit_reject_tier1': bool(reject[u]),
                'cluster_size': len(members), 'cluster_mass': mass,
                'cluster_pval': p_clust, 'cluster_reject_tier2': p_clust < ALPHA,
                'fold_enrichment': fold_enrichment[u],
            })

    print(f"# Tier 2: topology-aware cluster enrichment (cluster-forming p < {CLUSTER_FORMING_ALPHA})")
    if not cluster_rows:
        print("# No units passed the cluster-forming threshold - no candidate regions to test.")
        return None

    cluster_df = pd.DataFrame(cluster_rows).sort_values(['cluster_pval', 'cluster_id'])

    n_clusters = cluster_df['cluster_id'].nunique()
    n_sig_clusters = cluster_df.loc[cluster_df['cluster_reject_tier2'], 'cluster_id'].nunique()
    rescued = cluster_df[cluster_df['cluster_reject_tier2'] & ~cluster_df['unit_reject_tier1']]
    print(f"# {n_clusters} candidate region(s), {n_sig_clusters} significant at FWER < {ALPHA}")
    print(f"# {len(rescued)} unit(s) significant via their neighborhood in tier2 but not individually in tier1")

    return cluster_df

def main(args):
    parser = argparse.ArgumentParser(prog="som-test", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-s", "--som", required=True,
                        help="Built SOM")
    parser.add_argument("-n", "--null-map", required=True,
                        help="Null map of loci to SOM")
    parser.add_argument("-m", "--subset-map", required=True,
                        help="Subset map of loci to SOM")
    parser.add_argument("-o", "--output", default="/dev/stdout",
                         help="Output tsv (%(default)s)")
    parser.add_argument("-a", "--alpha", type=float, default=0.05,
                        help="Signifance threshold (%(default)s)")
    parser.add_argument("-R", type=int, default=10000,
                        help="Permutation count for T1 (%(default)s)")
    parser.add_argument("-R2", type=int, default=10000,
                        help="Permutation count for T2 (%(default)s)")
    parser.add_argument("--dist-prune", type=int, 
                        help="Percentile for adjacency matrix pruning (0-100; %(default)s)")
    args = parser.parse_args(args)

    som = pickle.load(open(args.som, 'rb'))

    shape = som['som'].distance_map().shape

    columns = ['chrom', 'start', 'end', 'X', 'Y']
    all_trs = pd.read_csv(args.null_map, sep='\t', names=['X', 'Y'])
    all_tr_count = square(all_trs, shape)

    subset_trs = pd.read_csv(args.subset_map, sep='\t', names=columns)
    subset_tr_count = square(subset_trs, shape)

    clusters = som_test(all_tr_count, subset_tr_count, som['som'],
                        ALPHA=args.alpha, 
                        R=args.R, 
                        CLUSTER_FORMING_ALPHA=0.05, 
                        R2=args.R2, 
                        DIST_PRUNE_PERCENTILE=args.dist_prune)
    
    clusters.to_csv(args.output, sep='\t', index=False, float_format='%.5f')


if __name__ == '__main__':
    main(sys.argv[1:])
