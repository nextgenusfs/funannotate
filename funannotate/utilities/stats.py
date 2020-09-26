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
    parser = argparse.ArgumentParser(prog='stats.py',
                                     description='''Script to run some simple genome annotation stats''',
                                     epilog="""Written by Jon Palmer (2020) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-f', '--fasta', required=True,
                        help='Genome in FASTA format')
    parser.add_argument('-o', '--out', required=True,
                        help='JSON output stats file')
    parser.add_argument('-g', '--gff',
                        help='Genome annotation in GFF3 format')
    parser.add_argument('-t', '--tbl',
                        help='Genome annotation in TBL format')
    parser.add_argument('--transcript_alignments',
                        help='transcript alignments in GFF3 format')
    parser.add_argument('--protein_alignments',
                        help='protein alignments in GFF3 format')
    args = parser.parse_args(args)


    if not args.gff and not args.tbl:
        print('Warning: no genome annotation passed (-t or -g), will only output genome assembly stats')
    elif args.tbl:
        print('Generating stats from Genome FASTA file and TBL annotation')
        lib.annotation_summary(args.fasta, args.out, tbl=args.tbl,
                               transcripts=args.transcript_alignments,
                               proteins=args.protein_alignments)
    elif args.gff:
        print('Generating stats from Genome FASTA file and GFF3 annotation')
        lib.annotation_summary(args.fasta, args.out, gff=args.gff,
                               transcripts=args.transcript_alignments,
                               proteins=args.protein_alignments)
    print('Finished writing JSON stats file: {}'.format(args.out))


if __name__ == "__main__":
    main(sys.argv[1:])
