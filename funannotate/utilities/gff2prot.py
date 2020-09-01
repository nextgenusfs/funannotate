#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from natsort import natsorted
import funannotate.library as lib


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='gff2prot.py',
                                     description='''Script to convert GFF3 and FASTA proteins.''',
                                     epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-g', '--gff3', required=True,
                        help='Genome annotation GFF3 format')
    parser.add_argument('-f', '--fasta', required=True,
                        help='Genome in FASTA format')
    parser.add_argument('--no_stop', action='store_true',
                        help='Dont print stop codon')
    args = parser.parse_args(args)

    # translate GFF3 to proteins
    # load into dictionary
    Genes = {}
    Genes = lib.gff2dict(args.gff3, args.fasta, Genes)

    for k, v in natsorted(list(Genes.items())):
        if v['type'] == 'mRNA':
            for i, x in enumerate(v['ids']):
                if args.no_stop:
                    Prot = v['protein'][i].rstrip('*')
                else:
                    Prot = v['protein'][i]
                sys.stdout.write('>%s %s\n%s\n' % (x, k, lib.softwrap(Prot)))


if __name__ == "__main__":
    main(sys.argv[1:])
