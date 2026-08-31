"""
Map a kvec to a som. Outputs a tsv of `X\tY` for each entry in the kvec
"""
import pickle
import argparse
import numpy as np
from tqdm import tqdm

def map_kvec(args):
    """
    Map a kvec onto a SOM. Just outputs list of coordinates
    """
    parser = argparse.ArgumentParser(prog="map-kvec", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-k", "--kvec", required=True,
                        help="Input kmer vectors")
    parser.add_argument("-s", "--som", required=True,
                        help="Input SOM")
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="Output `.map` %(default)s")
    args = parser.parse_args(args)

    data = np.load(args.kvec)
    som = pickle.load(open(args.som, 'rb'))
    fout = open(args.output, 'w')
    for i in tqdm(data['kmers']):
        print(*som['som'].winner(i), sep='\t', file=fout)
    fout.close()
