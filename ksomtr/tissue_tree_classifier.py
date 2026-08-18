"""
tissue_tree_classifier.py

Finds, for each locus, the smallest node in a core -> tissue -> layer
hierarchy that "best explains" the set of (donor, core) observations
showing somatic instability, trading off specificity (depth) against
how many observations must be treated as exceptions.

This is the tree-DP version discussed as a lighter-weight alternative to a
full MaxSAT encoding -- it's exact, fast, and easy to audit by hand.

--------------------------------------------------------------------------
HOW TO PLUG IN YOUR DATA
--------------------------------------------------------------------------
You need two inputs:

1. A tree specification: a table mapping each core to its tissue, and each
   tissue to its layer (or however many levels you actually have). See
   `build_tree_from_table()` below -- edit the column names to match yours.

2. An observations table: one row per locus x donor x core that showed
   instability. At minimum needs columns: locus_id, donor_id, core_id.
   (You can derive this from your Table 1 by whatever threshold you're
   using to call "unstable" -- this script starts downstream of that call.)

Search for "EDIT ME" for the places you'll most likely need to adapt.
--------------------------------------------------------------------------
"""
from __future__ import annotations

import os
import sys
import math
from dataclasses import dataclass, field
from typing import Optional
from copy import deepcopy

import numpy as np
import pandas as pd
import smahtkit as sk


# ============================================================
# STAGE 0: Data structures
# ============================================================

@dataclass
class Node:
    id: str
    level: str                      # e.g. "core", "tissue", "layer"
    parent: Optional["Node"] = None
    children: list = field(default_factory=list)
    depth: int = 0                  # 0 = root (layer), increases toward core
    n_leaves: int = 0               # number of core-level descendants (or 1 if leaf)

    # populated per-locus by annotate_coverage()
    covered_obs: set = field(default_factory=set)   # {(donor_id, core_id), ...}

    @property
    def is_leaf(self) -> bool:
        return len(self.children) == 0

    def covered_donors(self) -> set:
        return {d for d, c in self.covered_obs}

    def __repr__(self):
        return f"Node({self.level}:{self.id}, depth={self.depth})"


class Tree:
    def __init__(self):
        self.nodes: dict[str, Node] = {}   # id -> Node
        self.root: Optional[Node] = None

    def add_node(self, node_id: str, level: str, parent_id: Optional[str] = None) -> Node:
        if node_id in self.nodes:
            return self.nodes[node_id]
        node = Node(id=node_id, level=level)
        self.nodes[node_id] = node
        if parent_id is not None:
            parent = self.nodes[parent_id]
            node.parent = parent
            parent.children.append(node)
            node.depth = parent.depth + 1
        else:
            self.root = node
        return node

    def all_nodes(self):
        return list(self.nodes.values())

    def finalize(self):
        """Call once after adding all nodes: computes depth (if root added
        last) and n_leaves for every node. Safe to call multiple times."""
        # recompute depth from root via BFS in case nodes were added out of order
        if self.root is not None:
            self.root.depth = 0
            stack = [self.root]
            while stack:
                n = stack.pop()
                for c in n.children:
                    c.depth = n.depth + 1
                    stack.append(c)

        def count_leaves(node: Node) -> int:
            if node.is_leaf:
                node.n_leaves = 1
                return 1
            total = sum(count_leaves(c) for c in node.children)
            node.n_leaves = max(total, 1)  # avoid 0 for empty branches
            return node.n_leaves

        for n in self.nodes.values():
            if n.parent is None:
                count_leaves(n)


def build_tree_from_table(mapping_df: pd.DataFrame,
                           level_cols: list[str]) -> Tree:
    """
    EDIT ME: adjust level_cols to match your table's column names, ordered
    from coarsest to finest, e.g. ["layer_id", "tissue_id", "core_id"].

    mapping_df: one row per core, with a column for every level giving that
    core's ancestry. E.g.:

        layer_id    tissue_id   core_id
        ectoderm    brain       cortex_A
        ectoderm    brain       cerebellum_A
        mesoderm    adrenal     adrenal_A

    Duplicate ancestry rows (e.g. multiple donors sharing the same core
    label) are fine -- they'll just be deduped when building the tree.
    """
    tree = Tree()
    # add a synthetic super-root above your top level so the DP always has
    # a single top node to consider (cost of this node = "no coherent
    # explanation, everything is an exception")
    tree.add_node("__ROOT__", level="root", parent_id=None)

    for _, row in mapping_df.drop_duplicates(subset=level_cols).iterrows():
        parent_id = "__ROOT__"
        for i, col in enumerate(level_cols):
            level_name = col  # or hardcode nicer names
            node_id = str(row[col])
            tree.add_node(node_id, level=level_name, parent_id=parent_id)
            parent_id = node_id

    tree.finalize()
    return tree


# ============================================================
# STAGE 1: Annotate every node with what it covers, for one locus
# ============================================================

def annotate_coverage(tree: Tree, observations: set[tuple[str, str]]):
    """
    observations: set of (donor_id, core_id) pairs for a single locus.
    Populates node.covered_obs for every node in the tree (bottom-up).
    Must be re-run for each locus (it overwrites covered_obs).
    """
    # clear previous locus's annotation
    for n in tree.nodes.values():
        n.covered_obs = set()

    # seed leaves
    obs_by_core: dict[str, set] = {}
    for donor_id, core_id in observations:
        obs_by_core.setdefault(core_id, set()).add((donor_id, core_id))

    for n in tree.nodes.values():
        if n.is_leaf:
            n.covered_obs = obs_by_core.get(n.id, set())

    # roll up bottom-up; process nodes in decreasing depth order so
    # children are finalized before parents
    for n in sorted(tree.nodes.values(), key=lambda x: -x.depth):
        if not n.is_leaf:
            covered = set()
            for c in n.children:
                covered |= c.covered_obs
            n.covered_obs = covered


# ============================================================
# STAGE 2: Bottom-up DP -- cost of "claiming" each node
# ============================================================

def depth_penalty(node: Node, tree: Tree, mode: str = "linear") -> float:
    """
    Cost of specificity for claiming this node as the explanation.
    Lower depth (closer to root) = higher penalty (less specific).

    mode="linear": penalty = max_depth - node.depth
    mode="mdl":     penalty = log2(total_leaves_in_tree / node.n_leaves)
                    (bits needed to point at this node's subtree --
                     penalizes claiming a big bushy node more than a
                     small one at the same depth)
    """
    if mode == "linear":
        max_depth = max(n.depth for n in tree.nodes.values() if n.covered_obs)
        return max_depth - node.depth
    elif mode == "mdl":
        total_leaves = tree.root.n_leaves if tree.root.n_leaves > 0 else 1
        n_leaves = max(node.n_leaves, 1)
        return math.log2(total_leaves / n_leaves)
    else:
        raise ValueError(f"unknown depth_penalty mode: {mode}")


_DONOR_AGG_FUNCS = {
    "max": max,
    "sum": sum,
    "mean": lambda xs: sum(xs) / len(xs),
}


def cost(node: Node, obs_weights: dict[tuple[str, str], float], lam: float,
         tree: Tree, penalty_mode: str = "linear",
         donor_agg: str = "max") -> tuple[float, set, float]:
    """
    Returns (cost, excluded_obs, excluded_weight).

    obs_weights: {(donor_id, core_id): weight} for every observation at
    this locus. Weight is your confidence in that observation being real
    signal rather than noise -- e.g. VAF, read-support-derived score, or
    1 - p_batch_effect. Defaults to 1.0 per observation if you don't have
    a weight (see `to_obs_weights()` helper below).

    Exclusion cost is: for each donor with >=1 excluded observation,
    aggregate that donor's excluded weights (default: max, so one strong
    off-target hit costs more to exclude than several weak ones from the
    same donor -- switch donor_agg="sum" if you want repeated low-VAF
    hits in the same donor to add up instead).
    """
    all_obs = set(obs_weights.keys())
    excluded = all_obs - node.covered_obs

    donor_excluded_weights: dict[str, list[float]] = {}
    for (d, c) in excluded:
        donor_excluded_weights.setdefault(d, []).append(obs_weights[(d, c)])

    agg_func = _DONOR_AGG_FUNCS[donor_agg]
    excluded_weight = sum(agg_func(ws) for ws in donor_excluded_weights.values())

    penalty = depth_penalty(node, tree, mode=penalty_mode)
    return penalty + lam * excluded_weight, excluded, excluded_weight


def normalize_weight(raw_value: float, saturate_at: float, floor: float = 0.0) -> float:
    """
    EDIT ME: maps a raw signal (VAF, read count, 1-p, etc.) onto a 0-1
    confidence scale so it's on comparable footing with the old integer
    donor-count penalty. `saturate_at` is the value at which you consider
    the signal "fully trustworthy" (weight -> 1.0); values above it are
    clipped to 1.0. `floor` sets a minimum weight for any observation
    that was included in obs_df at all (use 0.0 to let near-zero VAF
    observations count for almost nothing).

    Example: normalize_weight(vaf, saturate_at=0.2) turns VAF >= 0.2 into
    weight 1.0, and scales anything below that linearly toward `floor`.
    """
    if saturate_at <= 0:
        raise ValueError("saturate_at must be > 0")
    scaled = raw_value / saturate_at
    return max(floor, min(scaled, 1.0))


def to_obs_weights(observations: set[tuple[str, str]],
                    weight_lookup: Optional[dict[tuple[str, str], float]] = None
                    ) -> dict[tuple[str, str], float]:
    """Convenience: turn a plain observation set into a weight dict, using
    1.0 for any pair missing from weight_lookup (or all pairs, if no
    lookup given -- reproduces the old unweighted behavior)."""
    weight_lookup = weight_lookup or {}
    return {obs: weight_lookup.get(obs, 1.0) for obs in observations}


# ============================================================
# STAGE 3: Search all nodes for the best single-node explanation
# ============================================================

def best_nodes(tree: Tree, obs_weights: dict, lam: float,
                penalty_mode: str = "linear",
                donor_agg: str = "max") -> tuple[list[Node], float]:
    candidates = [n for n in tree.nodes.values() if len(n.covered_donors()) > 0]
    if not candidates:
        return [], float("inf")

    scored = sorted([(n, cost(n, obs_weights, lam, tree, penalty_mode, donor_agg)[0])
                     for n in candidates], key=lambda x: x[1], reverse=True)
    min_cost = min(s for _, s in scored)
    # tolerate noise...?
    scores = np.array([s for _, s in scored])
    thresh = np.std(scores) / 4 # Need a more rigorous thing here
    # And I need to 'chain', like if there's __ROOT_ and others, go down. If there's ADGL/ADGR, go up?
    winners = [n for n, s in scored if abs(s - min_cost) < thresh]
    return winners, min_cost


# ============================================================
# STAGE 4: Tie-handling / disjunctive fallback
# ============================================================

@dataclass
class LocusResult:
    kind: str                  # "single", "disjunctive", "ambiguous"
    winners: list[Node]
    cost: float
    excluded_obs: set = field(default_factory=set)
    excluded_weight: float = 0.0


def classify_locus(tree: Tree, obs_weights: dict, lam: float,
                    penalty_mode: str = "linear",
                    donor_agg: str = "max") -> LocusResult:
    winners, min_cost = best_nodes(tree, obs_weights, lam, penalty_mode, donor_agg)
    if not winners:
        return LocusResult(kind="no_evidence", winners=[], cost=min_cost)

    if len(winners) == 1:
        _, excl, excl_w = cost(winners[0], obs_weights, lam, tree, penalty_mode, donor_agg)
        return LocusResult("single", winners, min_cost, excl, excl_w)

    parents = {w.parent.id if w.parent else None for w in winners}
    kind = "disjunctive" if len(parents) == 1 else "ambiguous"
    _, excl, excl_w = cost(winners[0], obs_weights, lam, tree, penalty_mode, donor_agg)
    return LocusResult(kind, winners, min_cost, excl, excl_w)


# ============================================================
# STAGE 5: Lambda sweep for stability / confidence
# ============================================================

def classify_locus_with_confidence(tree: Tree, obs_weights: dict,
                                    lam_range: list[float],
                                    lam_default: float,
                                    penalty_mode: str = "linear",
                                    donor_agg: str = "max"):
    """
    Sweeps lambda, tracks which node(s) win at each value, and reports
    the width of the stable interval containing lam_default as a rough
    substitute for a confidence/significance measure.
    """
    results = [(lam, classify_locus(tree, obs_weights, lam, penalty_mode, donor_agg))
               for lam in lam_range]

    def winner_key(res: LocusResult):
        return tuple(sorted(n.id for n in res.winners))

    # find the contiguous interval of lam_range that shares the same
    # winner_key as lam_default's result
    default_result = classify_locus(tree, obs_weights, lam_default, penalty_mode, donor_agg)
    default_key = winner_key(default_result)

    # locate index closest to lam_default
    idx_default = min(range(len(lam_range)), key=lambda i: abs(lam_range[i] - lam_default))

    lo = idx_default
    while lo > 0 and winner_key(results[lo - 1][1]) == default_key:
        lo -= 1
    hi = idx_default
    while hi < len(results) - 1 and winner_key(results[hi + 1][1]) == default_key:
        hi += 1

    interval_width = lam_range[hi] - lam_range[lo]
    total_width = lam_range[-1] - lam_range[0] if lam_range[-1] != lam_range[0] else 1.0
    confidence = interval_width / total_width

    return default_result, confidence, results


# ============================================================
# STAGE 6: Run across all loci
# ============================================================

def run_all_loci(tree: Tree,
                  obs_df: pd.DataFrame,
                  locus_col: str = "locus_id",
                  donor_col: str = "donor_id",
                  core_col: str = "core_id",
                  weight_col: Optional[str] = None,
                  lam_range: list[float] = None,
                  lam_default: float = 1.0,
                  penalty_mode: str = "linear",
                  donor_agg: str = "max") -> pd.DataFrame:
    """
    obs_df: one row per (locus, donor, core) observation of instability.

    weight_col: EDIT ME -- name of a column in obs_df holding your
    per-observation confidence weight (e.g. VAF, or 1 - p_batch_effect).
    If None, every observation is weighted 1.0 (old unweighted behavior).
    If a (locus, donor, core) combo appears more than once in obs_df with
    different weights, the max is taken.

    donor_agg: how a donor's multiple excluded observations combine into
    that donor's contribution to the exclusion penalty. "max" (default)
    means a donor's cost to exclude is driven by their single strongest
    excluded signal; "sum" makes repeated hits from the same donor add up.
    """
    if lam_range is None:
        lam_range = [round(x, 2) for x in _frange(0.1, 5.0, 0.1)]

    rows = []
    for locus_id, group in obs_df.groupby(locus_col):
        if weight_col is not None:
            obs_weights: dict = {}
            for _, r in group.iterrows():
                key = (r[donor_col], r[core_col])
                w = float(r[weight_col])
                obs_weights[key] = max(w, obs_weights.get(key, w))
        else:
            obs = set(zip(group[donor_col], group[core_col]))
            obs_weights = to_obs_weights(obs)

        annotate_coverage(tree, set(obs_weights.keys()))
    
        result, confidence, _ = classify_locus_with_confidence(
            tree, obs_weights, lam_range, lam_default, penalty_mode, donor_agg
        )

        winner_ids = ",".join(n.id for n in result.winners) if result.winners else None
        winner_level = result.winners[0].level if result.winners else None

        rows.append({
            "locus_id": locus_id,
            "n_observations": len(obs_weights),
            "n_donors": len({d for d, c in obs_weights}),
            "total_weight": round(sum(obs_weights.values()), 3),
            "call_kind": result.kind,
            "winning_node(s)": winner_ids,
            "winning_level": winner_level,
            "cost_at_default_lambda": round(result.cost, 3) if result.cost != float("inf") else result.cost,
            "excluded_weight": round(result.excluded_weight, 3),
            "confidence (stable lambda fraction)": round(confidence, 3),
        })

    return pd.DataFrame(rows).sort_values("locus_id").reset_index(drop=True)


def _frange(start, stop, step):
    n = int(round((stop - start) / step))
    return [start + i * step for i in range(n + 1)]


# ============================================================
# EXAMPLE / SELF-TEST with synthetic data
# ============================================================

if __name__ == "__main__":
    # Build Tree from the Protocols
    # layer, dedup_tissue, tissue
    rows = []
    for protocol, val in sk.smhtid.PROTOCOLS.items():
        rows.append([val['layer'], 
                     sk.smhtid.DEDUP_LOOKUP[val['DedupKey']]['tissue_abv'], 
                     val['tissue_abv']])

    levels = ['layer', 'dedup_tissue', 'tissue']
    mapping_df = pd.DataFrame(rows, columns=levels)
    tree = build_tree_from_table(mapping_df, levels)

    rows = []
    for fn in sys.argv[1:]:
        reads = pd.read_csv(fn, sep='\t')
        reads['sample'] = reads['full_name'].apply(sk.SMaHTid)
        reads['donor'] = reads['sample'].apply(lambda x: x.donor)
        reads['tissue'] = reads['sample'].apply(lambda x: x.tissue_abv)
        # Use unstable to figure out which samples to actually deal with
        other = os.path.join(os.path.dirname(fn), "output.unstable.tsv")
        other = pd.read_csv(other, sep='\t')
        other['tissue'] = other['protocol'].apply(lambda x: sk.smhtid.PROTOCOLS[x]['tissue_abv'])

        other['mag'] = ((other['delta_mid']**2 + other['spread']**2) ** 0.5) * other['vaf']
        scaled = other['mag'] / np.percentile(other['mag'], 90)
        scaled[scaled > 1.0] = 1.0
        other['weight'] = scaled

        other.set_index(['donor', 'tissue', 'hap'], inplace=True)
        reads = reads.groupby(['donor', 'tissue', 'hap', 'is_germ']).size().unstack()
        reads = reads[reads.index.isin(other.index)]
        reads['weight'] = other['weight']
        #reads[False].fillna(0) / reads.sum(axis=1)
        # Assuming a per-locus output
        locus = os.path.dirname(fn).split('/')[-1]
        for (donor, tissue, hap), row in reads.iterrows():
            rows.append((locus, f'{donor}.{hap}', tissue, row['weight']))
    obs = pd.DataFrame(rows, columns=['locus_id', 'donor_id', 'tissue', 'weight'])

    #obs["vaf_norm"] = obs["vaf"].apply(lambda v: normalize_weight(v, saturate_at=50))
    # Should maybe put lam_default to some weight of how many tissues are present?
    summary_weighted_norm = run_all_loci(
        tree, obs, weight_col="weight", core_col='tissue', donor_agg='max', lam_default=0.5,
        penalty_mode='linear'
    )
    summary_weighted_norm.to_csv("tissue_classifier.tsv", sep='\t', index=False)
