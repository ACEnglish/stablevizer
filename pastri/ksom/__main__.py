import argparse
import importlib


# Lazy import 
CMDS = [
    'build-kvec',
    'build-som',
    'map-kvec',
    'map-bed',
    'seqstat-bed',
    'plot-som',
    'downsample',
]

def main(args):
    """
    Main
    """
    desc = "Kmer based SOM scripts\nAvailable commands:\n\t" + "\n\t".join(CMDS)
    parser = argparse.ArgumentParser(prog="svpipe", description=desc,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("command", metavar="CMD", choices=CMDS, type=str,
                        help="Command to run")
    parser.add_argument("options", metavar="OPTS", nargs=argparse.REMAINDER,
                        help="Options to pass to command")
    args = parser.parse_args(args)
    func = args.command.replace('-', '_')
    path = f'pastri.ksom.{func}'
    cmd = getattr(importlib.import_module(path), func)

    cmd(args.options)

if __name__ == '__main__':
    main()
