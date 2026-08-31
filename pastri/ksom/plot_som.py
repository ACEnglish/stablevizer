"""
Plot the SOM
"""
import pickle
import argparse

import minisom
import numpy as np
import pandas as pd
from tqdm import tqdm
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
from matplotlib.patches import RegularPolygon, Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm, colorbar
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize
from shapely.geometry import Polygon
from shapely.ops import unary_union


def som_plot(som,
             heatmap=None,
             heatmap_label="UMatrix",
             color_map=cm.Blues,
             norm=None,
             enrichment_df=None):
    """
    Main plotting function of SOM
    """
    xx, yy = som.get_euclidean_coordinates()
    weights = som.get_weights()

    if heatmap is None:
        heatmap = som.distance_map()
        
    if norm is None:
        norm = Normalize(vmin=np.nanmin(heatmap), vmax=np.nanmax(heatmap))
        
    f = plt.figure(figsize=(10,10), dpi=180)
    ax = f.add_subplot(111)

    ax.set_aspect('equal')

    if enrichment_df is not None and not enrichment_df.empty:
        to_highlight = enrichment_df.set_index(['X', 'Y']).index
    else:
        to_highlight = None

    # iteratively add hexagons
    all_hex = {}
    for i in range(weights.shape[0]):
        for j in range(weights.shape[1]):
            unit_sig = to_highlight is None or (i,j) in to_highlight
            wy = yy[(i, j)] * np.sqrt(3) / 2
            cx, cy = xx[(i, j)], wy
            r = 0.87 / np.sqrt(3) if unit_sig else 0.70 / np.sqrt(3)
            # Scale X relative to origin, then translate to (cx, cy)
            hex_trans = (
                mtransforms.Affine2D().scale(2 / np.sqrt(3), 1.0).translate(cx, cy)
                + ax.transData
            )
            hex_patch = RegularPolygon((0, 0),
                                       numVertices=6,
                                       radius=r,
                                       transform=hex_trans,
                                       facecolor=color_map(norm(heatmap[i, j])),
                                       alpha=0.6 - (int(not unit_sig) * 0.2),
                                       edgecolor='black' if unit_sig else 'gray')
            ax.add_patch(hex_patch)
            all_hex[(i, j)] = hex_patch

    # Draw Cluster Outer Boundaries
    if enrichment_df is not None and not enrichment_df.empty:
        uniq_clusters = enrichment_df['cluster_id'].unique()
        r_union = 1.0 / np.sqrt(3)  # Exact radius (no margin gap) to ensure proper polygon union
        # 0.5 for a turn
        angles = np.linspace(0.5, 2 * np.pi + 0.5, 6, endpoint=False)
        for cluster_id, group in enrichment_df.groupby('cluster_id'):
            if not group['cluster_reject'].any():
                continue
            polys = []

            for _, row in group.iterrows():
                i, j = int(row['X']), int(row['Y'])
                cx = xx[(i, j)]
                cy = yy[(i, j)] * np.sqrt(3) / 2
                verts = [(cx + r_union * np.cos(a), cy + r_union * np.sin(a)) for a in angles]
                polys.append(Polygon(verts))
                # Unset the edge
                all_hex[(i, j)].set_edgecolor(None)#"#000000")

            if polys:
                merged = unary_union(polys)
                geoms = merged.geoms if merged.geom_type == 'MultiPolygon' else [merged]

                for geom in geoms:
                    # Draw outer perimeter
                    x, y = geom.exterior.xy
                    ax.plot(x, y,
                            color='black',
                            linewidth=2.5,
                            linestyle='-')

                    # Draw inner holes (if any)
                    for interior in geom.interiors:
                        hx, hy = interior.xy
                        ax.plot(hx, hy,
                                color='black',
                                linewidth=2.5,
                                linestyle='-')


    plt.xticks([])
    plt.yticks([])
    buff = 0.55
    ax.set_xlim(xx.min() - buff, xx.max() + buff)
    ax.set_ylim(yy.min()*np.sqrt(3)/2 - buff, yy.max()*np.sqrt(3)/2 + buff)
    
    # Heatmap
    divider = make_axes_locatable(plt.gca())
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    cb1 = colorbar.ColorbarBase(ax_cb, cmap=color_map, norm=norm,
                                orientation='vertical', alpha=.4)
    cb1.ax.get_yaxis().labelpad = 16
    cb1.ax.set_ylabel(heatmap_label,
                      rotation=270, fontsize=16)
    plt.gcf().add_axes(ax_cb)
    plt.tight_layout()
    return f, ax


def plot_som(args):
    """
    Main entrypoint
    """
    parser = argparse.ArgumentParser(prog="plot-som", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    # TODO - this should be optional
    parser.add_argument("-b", "--bed", required=True,
                        help="Input `.map.bed`")
    parser.add_argument("-s", "--som", required=True,
                        help="Input `som.pkl`")
    parser.add_argument("-o", "--output", required=True,
                        help="Output figure `.png`")
    parser.add_argument("-t", "--title", default="SOM Plot",
                        help="Figure title (%(default)s)")
    parser.add_argument("-T", "--bar-title", default=None,
                        help="Color bar title (default `--metric`))")
    parser.add_argument("--metric", default="count", choices=['count', 'mean', 'sum', 'median'],
                        help="Default heatmap count, with other choices, expect sixth column")
    parser.add_argument("-e", "--enrichment", default=None,
                        help="Enrichment test result")
    parser.add_argument("--reverse", action="store_true",
                        help="Reverse color map")
    # `code/laytr/laytr/somplot.py`
    # TODO: Enable XYO Markers
    #parser.add_argument("-X", dtype=str)
    #parser.add_argument("-Y", dtype=str)
    #parser.add_argument("-O", dtype=str)

    args = parser.parse_args(args)
    
    if args.bar_title is None:
        args.bar_title = args.metric

    if args.enrichment:
        args.enrichment = pd.read_csv(args.enrichment, sep='\t')

    som = pickle.load(open(args.som, 'rb'))

    columns = ['chrom', 'start', 'end', 'X', 'Y']
    if args.metric != 'count':
        columns.append("sixth")

    wxy = pd.read_csv(args.bed, sep='\t',
                      names=columns
                     ).set_index(['chrom', 'start', 'end'])
    
    if args.metric == 'count':
        wxy['sixth'] = 1

    cells =  np.zeros(som['som'].distance_map().shape)
    l = wxy.groupby(['X', 'Y'])['sixth'].agg(args.metric)
    for pos, v in l.items():
        cells[pos] = v

    fig, ax = som_plot(som['som'],
                       heatmap=cells, 
                       heatmap_label=args.bar_title, 
                       color_map=cm.RdYlBu_r if args.reverse else cm.RdYlBu,
                       enrichment_df=args.enrichment)
    fig.suptitle(args.title)
    plt.savefig(args.output)

