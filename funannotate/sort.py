#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)

import sys
import argparse
from Bio import SeqIO
from funannotate.library import countfasta


def SortRenameHeaders(input, basename, output, minlen=False):
    # sort records and write temp file
    with open(output, 'w') as output:
        with open(input, 'r') as input:
            records = list(SeqIO.parse(input, 'fasta'))
            records.sort(cmp=lambda x, y: cmp(len(y), len(x)))
            counter = 1
            for rec in records:
                rec.name = ''
                rec.description = ''
                rec.id = basename + '_' + str(counter)
                if len(rec.id) > 16:
                    print(
                        'Error. {:} fasta header too long.  Choose a different --base name. NCBI/GenBank max is 16 characters'.format(rec.id))
                    sys.exit(1)
                if minlen:
                    if len(rec.seq) < int(minlen):
                        continue
                counter += 1
                SeqIO.write(rec, output, 'fasta')


def main(args):
        # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='sort_rename.py', usage="%(prog)s [options] -i genome.fa -o sorted.fa",
                                     description='''Script that sorts input by length and then renames contig headers.''',
                                     epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help='Multi-fasta genome file')
    parser.add_argument('-o', '--out', required=True,
                        help='Cleaned output (FASTA)')
    parser.add_argument('-b', '--base', default='scaffold',
                        help='Basename of contig header')
    parser.add_argument(
        '-m', '--minlen', help='Contigs shorter than threshold are discarded')
    args = parser.parse_args(args)

    print('{:,} contigs records loaded'.format(countfasta(args.input)))
    print("Sorting and renaming contig headers")
    if args.minlen:
        print("Removing contigs less than {:} bp".format(args.minlen))
    SortRenameHeaders(args.input, args.base, args.out, minlen=args.minlen)
    print('{:,} contigs saved to file'.format(countfasta(args.out)))


if __name__ == "__main__":
    main(sys.argv[1:])
