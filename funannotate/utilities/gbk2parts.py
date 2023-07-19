#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import funannotate.library as lib


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='gbk2parts.py',
                                     description='''Script to convert GBK file to its components.''',
                                     epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-g', '--gbk', required=True,
                        help='Genome in GenBank format')
    parser.add_argument('-o', '--output', required=True,
                        help='Output basename')
    args = parser.parse_args(args)

    # setup output files
    tblout = f'{args.output}.tbl'
    gffout = f'{args.output}.gff3'
    protout = f'{args.output}.proteins.fa'
    transout = f'{args.output}.mrna-transcripts.fa'
    cdsout = f'{args.output}.cds-transcripts.fa'
    dnaout = f'{args.output}.scaffolds.fa'
    lib.gb2parts(args.gbk, tblout, gffout, protout, transout, cdsout, dnaout)

if __name__ == "__main__":
    main(sys.argv[1:])
