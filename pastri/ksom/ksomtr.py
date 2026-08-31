import sys
import pickle
import argparse

import pysam
import minisom
from truvari.annotations.lcr import sequence_entropy
import numpy as np
import pandas as pd
from tqdm import tqdm

def sample_farthest_point(df, cols, n_samples, random_state=None):
    rng = np.random.default_rng(random_state)
    points = df[cols].to_numpy()
    n = len(points)

    mask = np.ones(n, dtype=bool)
    min_dists = np.full(n, np.inf)

    start = rng.integers(n)
    selected_idx = [start]
    mask[start] = False
    tot = n_samples - 1
    for _ in tqdm(range(tot), total=tot):
        last_point = points[selected_idx[-1]]
        dists = np.linalg.norm(points - last_point, axis=1)
        min_dists = np.minimum(min_dists, dists)

        masked_dists = np.where(mask, min_dists, -np.inf)
        next_idx = np.argmax(masked_dists)

        selected_idx.append(next_idx)
        mask[next_idx] = False

    return df.iloc[selected_idx]

def downsample(args):
    parser = argparse.ArgumentParser(prog="downsample")
    parser.add_argument("-b", "--bed-fn", required=True)
    parser.add_argument("-n", "--num", required=True)
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="%(default)s")
    args = parser.parse_args(args)
    
    # Should be asserting that this is the shape..
    df = pd.read_csv(in_fn, sep='\t', names=['chrom', 'start', 'end', 'a', 'b', 'c'])
    fps_sample = sample_farthest_point(df, ['a', 'b', 'c'], n_samples=NSAMP, random_state=8811)
    fps_sample.to_csv(args.output, sep='\t', index=False, header=False)



def ksom_main(args):

    parser = argparse.ArgumentParser(prog="ksomtr", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str, default=None,
                        help="Command to execute: "  + "\n".join(TOOLS.keys()))
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")
    args = parser.parse_args(args)
    TOOLS[args.cmd](args.options)
