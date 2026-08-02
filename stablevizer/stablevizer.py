"""
Tandem Repeat Instability Describer

Given a tsv of read lenghths, identify donor/tissues with instability and plot.
"""
import sys
import argparse
import itertools

import kmedoids
import numpy as np
import pandas as pd
import seaborn as sb

from tqdm import tqdm
import matplotlib.pyplot as plt
from sklearn.cluster import MeanShift, estimate_bandwidth
from sklearn.metrics import pairwise_distances

from stablevizer.protocols import PROTOCOLS

# Custom color palettes
PROTOCOL_PALETTE = {k:v['color'] for k,v in PROTOCOLS.items()}
TISSUE_PALETTE = {_['tissue_abv']:_['color'] for _ in PROTOCOLS.values()}
HAP_PALETTE = dict(zip([0, 1, 2], sb.color_palette("Set2", 3)))
THIRD_PALETTE = {False: 'gray', True: 'black'}

def parse_args(args):
    parser = argparse.ArgumentParser(prog="stablevizer", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("in_tsv", type=str,
                        help="Read lengths tsv")
    parser.add_argument("-o", "--output", default="output",
                        help="Output prefix (%(default)s)")
    parser.add_argument("-r", "--min-reads-donor", type=int, default=3,
                        help="Minimum number of non-germ reads required in donor (%(default)s)")
    parser.add_argument("-R", "--min-reads-tissue", type=int, default=3,
                        help="Minimum number of non-germ reads required in tissue (%(default)s)")
    parser.add_argument("-b", "--min-bandwidth", type=int, default=10,
                        help="Minimum MeanShift clustering bandwidth")
    parser.add_argument("-g", "--germ-vaf", type=float, default=0.8,
                        help="Minimum fraction of reads to collect germline cluster (%(default)s)")
    parser.add_argument("-q", "--germ-q", type=float, default=0.05,
                        help="Germline length interval for filtering haplotagging errors (%(default)s)")
    parser.add_argument("--rehaplotype", action="store_true",
                        help="Perform naive kmedoid clustering of readlengths to reassign haplotypes (%(default)s)")
    parser.add_argument("--abs-delta", action="store_true",
                        help="Calculate abs(∆) (%(default)s)")
    parser.add_argument("-t", "--title", default=None,
                        help="Locus detail to put in plot titles")
    parser.add_argument("--detail", type=str, default=None,
                        help="Produce individual donor detail plots (comma separated)")
    parser.add_argument("--all-detail", action="store_true",
                        help="Produce a large plot with all unstable donor's read details (%(default)s)")
    parser.add_argument("--ALL-detail", action="store_true",
                        help="Produce a large plot with ALL donor's read details (%(default)s)")
    parser.add_argument("--subset", type=str, default=None,
                        help="Subset reads to donors (comma separated)")
    parser.add_argument("--pdf", action="store_true",
                        help="Save instability plot as pdf (%(default)s)")
    parser.add_argument("--save-reads", action="store_true",
                        help="Save the annotated reads to `.anno_reads.tsv`")

    args = parser.parse_args(args)
    if (args.all_detail or args.ALL_detail) and (args.all_detail == args.ALL_detail):
        print("Only provide one of --all or --ALL detail", file=sys.stderr)
        sys.exit(1)

    if args.detail:
        args.detail = args.detail.split(',')
    else:
        args.detail = []

    return args


def locus_viz(data, donor="Sample", third='auto', fig=None, germ=None):
    """
    Make a per-sample 3-panel swarmplots of read length distributions
    
    When `third == 'auto'`, the last panel defaults to `is_{protocol}` where protocol is the most common ~is_germ protocol
    When `third == None`, the last panel is skipped.
    Otherwise, provide a pd.Series with `.name` of e.g. `is_skin` and boolan values for custom third plot

    Provide an existing fig (ax) to fill-in to an existing plot. Otherwise a new `pls.subplots` is made.
    """
    n_col = 2
    if isinstance(third, str) and third == 'auto':
        most = data[~data['is_germ']]['protocol'].value_counts().idxmax()
        third = data['protocol'] == most
        third.name = f'is_{most}'
        n_col += 1
        
    if fig is None:
        fig, ax = plt.subplots(ncols=n_col, sharey=True)
    else:
        ax = fig.subplots(ncols=n_col, sharey=True)

    np.random.seed(123)
    sb.stripplot(data=data, x='is_germ', y='length',
                hue='hap', 
                 ax=ax[0],
                 palette=HAP_PALETTE,
                 zorder=1)
    
    np.random.seed(123)
    sb.stripplot(data=data, x='is_germ', y='length',
                hue='protocol',
                 ax=ax[1],
                 palette=PROTOCOL_PALETTE,
                 zorder=1)
    ax[1].legend().remove()

    np.random.seed(123)
    if third is not None:
        sb.stripplot(data=data, x='is_germ', y='length',
                     hue=third, 
                     ax=ax[2], 
                     palette=THIRD_PALETTE,
                     zorder=1,
                     dodge=True,
                     hue_order=[False, True])
    fig.suptitle(f'{donor} TR View')


    if germ is not None:
        for i in ax:
            for _, row in germ[germ['hap'] != 0].iterrows():
                color = HAP_PALETTE[row['hap']]
                i.axhline(row['lower'], color=color, linestyle='--', zorder=20)
                i.axhline(row['upper'], color=color, linestyle='--', zorder=20)

    return (fig, ax)

def perform_clustering(data, min_bandwidth=10, germ_vaf=0.80, germ_q=0.05, absolute=False):
    """
    Updated data in place with new columns of:
      - `is_germ` boolean if the read is germline
      - `∆` the length difference of this read from its predicted germline length
    returns a summary DataFrame with per-haplotype clustering information of
      - `donor` being clustered
      - `hap` being clustered
      - `nclust` produced by MeanShift
      - `ngerm` clusters produced
      - `nassign` reads in germline clusters
      - `total` reads
      - `bandwidth` parameter used by MeanShift
      - `pct_germ` of all reads

    Note that these summary numbers are not updated by the `fix_haplotypes`, so are unreliable
    """
    result_parts = []
    germ_parts = []
    summary = []
    print(f"Clustering {len(data):,} reads across {data['donor'].nunique()} donors", file=sys.stderr)
    for _, sub in tqdm(data.groupby(['donor', 'hap'])):
        X = sub['length'].values.reshape((-1,1))

        # No cluster variability
        if len(X) == 1 or X.std() == 0:
            continue
            
        # Was playing with something
        #centers = sub.groupby('hap')['length'].describe().sort_values(by='mean')
        #if len(centers) > 2 and 1 in centers.index and 2 in centers.index:
        #    lower = centers[centers.index != 0].iloc[0]['75%']
        #    upper = centers[centers.index != 0].iloc[1]['25%']
        #    bw = (upper - lower - 2) / 2
        #if bw is None or bw < 0:
        
        bw = max(estimate_bandwidth(X, quantile=0.3, n_samples=500), min_bandwidth)
            
        m = MeanShift(bandwidth=bw).fit(X)
        labels = pd.Series(m.labels_, name='is_germ', index=sub.index)
        clusters = labels.value_counts()
        clusters.name = 'read_count'
        clusters = clusters.to_frame()
        
        germ_centroids = sub.groupby(['hap'])['length'].mean().values
        dists = np.abs(m.cluster_centers_.ravel()[:, None] - germ_centroids[None, :])
        clusters['distance'] = dists.min(axis=1)

        # Rank based sorting
        clusters['read_rank'] = clusters['read_count'].rank(ascending=False)
        clusters['dist_rank'] = clusters['distance'].rank(ascending=True)
        clusters['combined_rank'] = clusters[['read_rank', 'dist_rank']].sum(axis=1)
        clusters.sort_values(by='combined_rank', inplace=True)
        
        tot_reads = len(sub)
        assigned_reads = 0
        n_germ = - 1
        germ_labels = []
        
        clusters.sort_values(by='distance', inplace=True)
        while n_germ < len(clusters) and assigned_reads / tot_reads < germ_vaf:
            n_germ += 1
            assigned_reads += clusters.iloc[n_germ]['read_count']
            germ_labels.append(clusters.index[n_germ])

        summary.append([_[0],
                        _[1],
                        len(clusters),
                        len(germ_labels),
                        assigned_reads,
                        tot_reads, 
                        bw])

        assign = labels.isin(germ_labels).to_frame()
        germ_length = sub[assign['is_germ']]['length']
        germ_spread = pd.Series([_[1],
                                 germ_length.mean(),
                                 germ_length.quantile(germ_q),
                                 germ_length.quantile(1 - germ_q)
                                 ],
                                index=['hap', 'mean', 'lower', 'upper'],
                                name=_[0])
        assign['germ_length'] = germ_spread['mean']
        germ_spread.name = _[0]
        germ_spread['hap'] = _[1]
        result_parts.append(assign)
        germ_parts.append(germ_spread)
        
    summary = pd.DataFrame(summary, columns=['donor', 'hap', 'nclust', 'ngerm', 'nassign', 'total', 'bandwidth'])
    summary['pct_germ'] = summary['nassign'] / summary['total']

    result = pd.concat(result_parts)
    germ = pd.DataFrame(germ_parts).round().astype(int)
    germ.index.name = 'donor'
    data = data.join(result)
    data['is_germ'] = data['is_germ'].infer_objects(copy=False).fillna(True).astype(bool)
    
    fix_soma_haplotypes(data, germ)

    data['delta'] = data['length'] - data['germ_length']
    if absolute:
        data['delta'] = data['delta'].abs()

    return data, summary, germ

def fix_soma_haplotypes(data, germ_lookup):
    """
    The haplotype assignment might be bad, so for each donor we look at the
    germ_length distribution per haplotype (germ_lookup) and reassign any
    non-germ read whose length falls within [lower, upper] of that
    haplotype's distribution back to germ. We also flip its assigned
    haplotype (assumes hap in {1, 2}) and reset its expected germ_length.

    Mutates `data` in place.
    """
    read_cnt = 0
    donor_cnt = 0

    for donor, sub in data.groupby('donor'):
        # Sometimes there is no germline, I guess
        if donor not in germ_lookup.index:
            continue

        to_check = sub.loc[~sub['is_germ']].copy()
        if to_check.empty:
            continue

        # Gotta [[ to ensure a frame is returned
        donor_lookup = germ_lookup.loc[[donor]]

        change = pd.Series(False, index=to_check.index)
        for _, spread in donor_lookup.iterrows():
            if spread['hap'] == 0:
                continue
            change |= to_check['length'].between(spread['lower'], spread['upper'])

        if not change.any():
            continue

        to_check['hap'] = to_check['hap'].where(~change, to_check['hap'] % 2 + 1)
        to_check['is_germ'] = change
        to_check['germ_length'] = to_check['hap'].map(donor_lookup.set_index('hap')['mean'])

        data.loc[to_check.index, ['is_germ', 'hap', 'germ_length']] = (
            to_check[['is_germ', 'hap', 'germ_length']]
        )

        read_cnt += change.sum()
        donor_cnt += 1

    print(f"Masked {read_cnt} oddly haplotyped reads in {donor_cnt} donors within germ-q", file=sys.stderr)


def rehaplotype(data):
    """
    rehaplotype based on length
    """
    print("Rehaplotyping", file=sys.stderr)
    #reads_changed = 0
    #donors_changed = 0
    for _, sub in tqdm(data.groupby(['donor'])):
        dist = pairwise_distances(sub['length'].values.reshape((-1, 1)))
        med = kmedoids.KMedoids(2).fit(dist).labels_ + 1
        # Doesn't work - because there's no 1/2 matching
        #t = (sub['hap'] != med).sum()
        #reads_changed += t
        #if t:
            #donors_changed += 1
        data.loc[sub.index, 'hap'] = med
    #print("Changed {reads_changed} read haplotypes acros {donors_changed}", file=sys.stderr)
    return data

def instability_plot(filt_view, title=None, absolute=False):
    # Instability Plot
    fig, ax = plt.subplots(dpi=180)

    filt_view['tissue'] = filt_view['protocol'].apply(lambda x: PROTOCOLS[x]['tissue_abv'])

    order = filt_view['donor'].unique()
    marker_shapes = ['o', 's', '^', 'D', 'v', 'P', 'X', '*', 'h', '<', '>', 'p']
    # Ensure there's always enough with cycle
    markers = dict(zip(order, itertools.cycle(marker_shapes)))

    p = sb.scatterplot(data=filt_view, x='delta_mid', 
                       y='spread', 
                       hue='tissue', 
                       size='vaf', 
                       style='donor',
                       markers=markers,
                       palette=TISSUE_PALETTE,
                       edgecolors='black',
                       ax=ax,
                       zorder=2,
                       alpha=1)

    if not title:
        title = 'TR Instability Plot'

    xlim = list(p.get_xlim())
    ylim = list(p.get_ylim())
    if absolute:
        xlim[0] = 0
        ylim[0] = 0

    p.set(xlabel='Median somatic read-length ∆ from germline (bp)',
          ylabel='Max - Min somatic read-length (bp)',
          xlim=xlim,
          ylim=ylim,
          title=title)
    plt.grid(zorder=1)
    _ = plt.legend(bbox_to_anchor=(1, 1.02), fontsize=8)
    _ = plt.tight_layout()

    counts = filt_view['tissue'].value_counts().to_dict()
    counts.update(filt_view['donor'].value_counts().to_dict())
    
    header_map = {
        'tissue': 'Tissue (#donor)',
        'donor': 'Donor (#tissue)',
        'vaf': 'VAF'
    }

    for text in p.legend_.get_texts():
        label = text.get_text()
        if label in header_map:
            text.set_text(header_map[label])
            text.set_weight('bold')
        elif label in counts:
            text.set_text(f"{label} ({counts[label]})")
    
    return fig


def run_stablevizer(args):
    args = parse_args(args)
    data = pd.read_csv(args.in_tsv, sep='\t')
    expected = ['donor', 'protocol', 'hap', 'length']
    all_good = True
    for i in expected:
        if i not in data.columns:
            print("Expected column {i} not in {args.in_tsv}", file=sys.stderr)
            all_good = False

    avail_donor = data['donor'].unique()
    for i in args.detail:
        if not i in avail_donor:
            print("--detail {i} not in args.in_tsv}")
            all_good = False
            
    if not all_good:
        print("Fix input / args", file=sys.stderr)
        sys.exit(1)

    if args.subset:
        s = args.subset.split(',')
        print(f"Subsetting to {len(s)} donors", file=sys.stderr)
        data.drop(data.index[~data['donor'].isin(s)], inplace=True)
    
    if not len(data):
        print("No reads to cluster!!!", file=sys.stderr)
        exit(1)

    if args.rehaplotype:
        data = rehaplotype(data)
    else:
        unphased = data['hap'] == 0
        pct = (unphased.mean() * 100)
        print(f"Dropping {unphased.sum()} ({pct:.1f}%) unphased reads", file=sys.stderr)
        data.drop(data.index[unphased], inplace=True)

    data, clusters, germ = perform_clustering(data,
                                              min_bandwidth=args.min_bandwidth,
                                              germ_vaf=args.germ_vaf,
                                              germ_q=args.germ_q,
                                              absolute=args.abs_delta)
    
    idx = ['donor', 'hap']
    a = clusters.set_index(idx)
    b = germ.reset_index().set_index(idx)
    out = a.join(b)
    out.to_csv(f'{args.output}.germline.tsv', sep='\t')
    if args.save_reads:
        data.to_csv(f'{args.output}.anno_reads.tsv', sep='\t', index=False, float_format="%.1f")

    # We want at least a few non-germline reads per-donor - and for now we'll ignore hap=0 clusters
    m_filt = ~data['is_germ']
    alt_coverage = data[m_filt].groupby('donor').size()
    keep_donors = alt_coverage[alt_coverage >= args.min_reads_donor].index
    keep_donors = set(keep_donors)

    # Now we want to group the reads within donor by tissue 
    # You could be more clever here. Instead of all reads, we should give an opportunity
    # for each haplotype to have the instability.
    # This opens the door for potential reassignment of reads to the other haplotype
    grp = data[m_filt & data['donor'].isin(keep_donors)].groupby(['donor', 'protocol'])
    # Record the upper/lower ∆
    lower = grp['delta'].min()
    mean = grp['delta'].mean()
    upper = grp['delta'].max()
    spread = grp['length'].max() - grp['length'].min()
    alt_reads = grp.size()

    view = pd.concat([
        lower.round(1),
        mean.round(1),
        grp['delta'].quantile(.5).round(1),
        upper.round(1),
        spread.round(1),
        alt_reads],
        axis=1)

    view.columns = ['delta_min', 'delta_mean', 'delta_mid', 'delta_max', 'spread', 'alt_reads']

    # Record the per-protocol coverage for calculating VAF
    tot_reads = data.groupby(['donor', 'protocol']).size()
    tot_reads.name = 'coverage'
    view = view.join(tot_reads, how='left')

    view = view.reset_index()
    view['vaf'] = (view['alt_reads'] / view['coverage']).round(4)

    # Now subset to only donor/tissue with minimum somatic read support
    mask = (view['alt_reads'] >= args.min_reads_tissue) # & (view['mean'] > 20)
    filt_view = view[mask].copy()
    # Filtered Summary TSV
    print(f"Identified {filt_view['donor'].nunique()} donor / {filt_view['protocol'].nunique()} protocols with instability", file=sys.stderr)
    filt_view.to_csv(f"{args.output}.unstable.tsv", sep='\t', index=False)
    
    if len(filt_view):
        m_fig = instability_plot(filt_view, args.title, args.abs_delta)
        fmt = 'pdf' if args.pdf else 'png'
        plt.rcParams['pdf.fonttype'] = 42
        m_fig.savefig(f'{args.output}.instability.{fmt}', format=fmt, bbox_inches='tight')
    else:
        print(f"No instability detected. Skipping main plot", file=sys.stderr)
    
    if args.detail:
        print(f"Detailing {len(args.detail)} donors' TR", file=sys.stderr)
        avail_donors = data['donor'].unique()
        for donor in args.detail:
            if not donor in avail_donors:
                print(f"Donor {donor} not in reads. Cannot plot `--detail`", file=sys.stderr)
                continue
            m_fig, _ = locus_viz(data[data['donor'] == donor], donor, germ=germ.loc[[donor]])
            m_fig.savefig(f"{args.output}.{donor}_detail.png", bbox_inches='tight')
    
    if not (args.all_detail or args.ALL_detail):
        sys.exit(0)
    
    has_soma = filt_view['donor'].unique()
    uniq_donor = has_soma if args.all_detail else data['donor'].unique()
    print(f"Making all detail plot for {len(uniq_donor)} donors", file=sys.stderr)
    side  = int(np.ceil(np.sqrt(uniq_donor.size)))

    parent_fig = plt.figure(figsize=(8 * side, 4 * side))
    subfigs = parent_fig.subfigures(side, side)
    
    for subfig, donor in zip(subfigs.ravel(), sorted(uniq_donor)):
        third = 'auto' if donor in has_soma else None
        sub = data[data['donor'] == donor]
        locus_viz(sub, donor, third=third, fig=subfig, germ=germ.loc[[donor]])

    parent_fig.savefig(f"{args.output}.all_donor_detail.png", bbox_inches='tight')
