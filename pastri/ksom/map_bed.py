"""
Create kmer vectors of reference fasta bed regions and map them to a SOM
"""
import pickle
import argparse

import pysam
from tqdm import tqdm

from pastri.ksom.build_kvec import parse_bed_regions, seq_to_kvec

def map_bed(args):
    """
    Map a bed file onto a som
    """
    parser = argparse.ArgumentParser(prog="map-bed", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("-b", "--bed-fn", required=True,
                        help="Input bed file")
    parser.add_argument("-r", "--ref-fn", required=True,
                        help="Input reference fasta")
    parser.add_argument("-s", "--som-fn", required=True,
                        help="Input SOM")
    parser.add_argument("-o", "--output", default="/dev/stdout",
                        help="Output `.map.bed` %(default)s")

    args = parser.parse_args(args)

    regions = parse_bed_regions(args.bed_fn)
    ref = pysam.FastaFile(args.ref_fn)
    som = pickle.load(open(args.som_fn, 'rb'))
    fout = open(args.output, 'w')

    for region in tqdm(regions):
        seq = ref.fetch(*region)
        vec = seq_to_kvec(seq, som['k'], count=som['count'])
        w = som['som'].winner(vec)
        print(*region, *w, sep='\t', file=fout)
    fout.close()


