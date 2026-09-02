import sys
import argparse
from importlib.metadata import version

import pastri
from pastri import __version__
from pastri.msa import msa_main
from pastri.plume import plume_main
#from pastri.plump import plump_main
from pastri.origin import origin_main
from pastri.locus_setup import locus_setup_main
from pastri.ksom.__main__ import main as ksom_main
from pastri.smahtid import smahtid_main, quick_smahtid_main

def flat_version(args):
    """Print the version"""
    if len(args) and args[0].count("-v"):
        print(f"Pastri {version('pastri')}")
    else:
        print(f"Pastri v{__version__}")


TOOLS = {       
         "tr": locus_setup_main,
         "plume": plume_main,
         #"plump": plump_main,
         "origin": origin_main,
         "msa": msa_main,
         "ksom": ksom_main,
         "smahtid": smahtid_main,
         "qsmahtid": quick_smahtid_main,
         "version": flat_version,
}

USAGE = f"""\
Pastri v{__version__} - Pattern Analysis of Somatic TR Instability

    tr       Extract locus info from qdpi
    plume    Analyze plume patterns
    plump    Analyze plump patterns
    origin   Tissue origin and batch effect tests
    ksom     Kmer based self-organizing map
    smhtid   SMaHT nomenclature deconvolution
    version  Print the version and exit
"""

def main():
    """
    Main entrypoint for pastri tools
    """
    parser = argparse.ArgumentParser(prog="pastri", description=USAGE,
                            formatter_class=argparse.RawDescriptionHelpFormatter)

    parser.add_argument("cmd", metavar="CMD", choices=TOOLS.keys(), type=str, default=None,
                        help="Command to execute")
    parser.add_argument("options", metavar="OPTIONS", nargs=argparse.REMAINDER,
                        help="Options to pass to the command")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit()
    args = parser.parse_args()

    TOOLS[args.cmd](args.options)

if __name__ == '__main__':
    main()
