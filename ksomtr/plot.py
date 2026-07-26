import numpy as np
import pickle
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import RegularPolygon, Ellipse
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import cm, colorbar
from matplotlib.lines import Line2D
from matplotlib.colors import Normalize

def start_plot(som,  heatmap=None, heatmap_label="UMatrix", color_map=cm.Blues, norm=None):
    xx, yy = som.get_euclidean_coordinates()
    weights = som.get_weights()

    if heatmap is None:
        heatmap = som.distance_map()
        
    if norm is None:
        norm = Normalize(vmin=np.nanmin(heatmap), vmax=np.nanmax(heatmap))
        
    f = plt.figure(figsize=(10,10))
    ax = f.add_subplot(111)

    ax.set_aspect('equal')

    # iteratively add hexagons
    for i in range(weights.shape[0]):
        for j in range(weights.shape[1]):
            wy = yy[(i, j)] * np.sqrt(3) / 2
            hex = RegularPolygon((xx[(i, j)], wy),
                                 numVertices=6,
                                 radius=.95 / np.sqrt(3),
                                 facecolor=color_map(norm(heatmap[i, j])),
                                 alpha=.4,
                                 edgecolor='gray')
            ax.add_patch(hex)
    xrange = np.arange(weights.shape[0])
    yrange = np.arange(weights.shape[1])
    plt.xticks(xrange, xrange)
    plt.yticks(yrange * np.sqrt(3) / 2, yrange)
    buff = 0.55
    ax.set_xlim(xx.min() - buff, xx.max() + buff)
    ax.set_ylim(yy.min()*np.sqrt(3)/2 - buff, yy.max()*np.sqrt(3)/2 + buff)
    
    # Heatmap
    divider = make_axes_locatable(plt.gca())
    ax_cb = divider.new_horizontal(size="5%", pad=0.05)
    cb1 = colorbar.ColorbarBase(ax_cb, cmap=color_map, norm=norm,  # <-- pass norm here too
                                orientation='vertical', alpha=.4)
    cb1.ax.get_yaxis().labelpad = 16
    cb1.ax.set_ylabel(heatmap_label,
                      rotation=270, fontsize=16)
    plt.gcf().add_axes(ax_cb)
    plt.tight_layout()
    return f, ax

if __name__ == '__main__':

    model = "models/adotto.v2.norm4.som.pkl"
    mapping = "models/adotto.full_mapping.bed"
    seqstat = "models/adotto.v2.seqstat.bed"
    output_pfx = "models/adotto.v2.norm4"

    som = pickle.load(open(model, 'rb'))
    wxy = pd.read_csv(mapping, sep='\t',
                      names=['chrom', 'start', 'end', 'X', 'Y']
                     ).set_index(['chrom', 'start', 'end'])
    l = wxy.groupby(['X', 'Y']).size()
    counts =  np.zeros(som['som'].distance_map().shape)
    for pos, v in l.items():
        counts[pos] = v
        
    start_plot(som['som'], heatmap=counts, 
               heatmap_label="TR Count", 
               color_map=cm.RdYlBu)
    plt.savefig(output_pfx + '.count.png')

    trs = pd.read_csv(seqstat, sep='\t', 
                      names=['chrom', 'start', 'end', 'gc', 'hom', 'ent'],
                     ).set_index(['chrom', 'start', 'end'])

    trs = trs.join(wxy)

    l = trs.groupby(['X', 'Y'])['gc'].mean()
    gcpct =  np.zeros(som['som'].distance_map().shape)
    for pos, v in l.items():
        gcpct[pos] = v
        
    l = trs.groupby(['X', 'Y'])['hom'].mean()
    hompct =  np.zeros(som['som'].distance_map().shape)
    for pos, v in l.items():
        hompct[pos] = v
        
    l = trs.groupby(['X', 'Y'])['ent'].mean()
    entpct =  np.zeros(som['som'].distance_map().shape)
    for pos, v in l.items():
        entpct[pos] = v

    start_plot(som['som'], heatmap=gcpct, 
               heatmap_label="Mean GC Pct",
               color_map=cm.RdYlBu)
    plt.savefig(output_pfx + '.gc.png')

    start_plot(som['som'], heatmap=hompct, 
               heatmap_label="Percent Homopolymer",
               color_map=cm.RdYlBu)
    plt.savefig(output_pfx + '.hom.png')

    start_plot(som['som'], heatmap=entpct, 
               heatmap_label="Entropy",
               color_map=cm.RdYlBu)
    plt.savefig(output_pfx + '.ent.png')
    # Gotta do it the hard way for patho and stuff
    # And I still don't have markers
    # Or you do the work to make a better environment


