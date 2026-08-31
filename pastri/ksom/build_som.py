"""
Build a SOM
"""
import pickle
import argparse
import minisom
import numpy as np

def build_som(args):
    """
    Train the SOM on the kmers
    """
    parser = argparse.ArgumentParser(prog="build_som", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("kmers_fn",
                        help="Input kmer vectors.npz from `kvec` command")
    parser.add_argument("-o", "--out-fn", required=True,
                        help="Output SOM.pkl file")
    parser.add_argument("-s", "--sigma", type=float, default=1.5,
                        help="SOM sigma parameter %(default)s")
    parser.add_argument("-l", "--learning-rate", type=float, default=1,
                        help="SOM learning rate parameter %(default)s")
    parser.add_argument("-i", "--iters", type=int, default=1_000_000,
                        help="SOM training iterations %(default)s")
    parser.add_argument("-p", "--pca-init", action="store_true",
                        help="PCA initialization of neuron placement (%(default)s)")
    parser.add_argument("--seed", default=None)
    args = parser.parse_args(args)

    data = np.load(args.kmers_fn)
    # TODO: all of these could be parameters so different
    # SOMs could be tried more easily
    som = minisom.MiniSom(25, 25, data['kmers'].shape[1],
                          sigma=args.sigma,
                          learning_rate=args.learning_rate,
                          topology='hexagonal',
                          neighborhood_function='gaussian',
                          activation_distance='euclidean',
                          random_seed=args.seed
                         )

    if args.pca_init:
        som.pca_weights_init(data['kmers'])

    som.train_batch(data['kmers'], args.iters, verbose=True)
    
    output = {'som': som,
              'k': int(data['k']),
              'count': bool(data['count']),
              }

    with open(args.out_fn, 'wb') as fout:
        pickle.dump(output, fout)
