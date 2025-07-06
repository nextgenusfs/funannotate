#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import absolute_import, division, print_function, unicode_literals

import sys
import argparse
from Bio.SeqIO.FastaIO import SimpleFastaParser
from funannotate.library import countfasta, softwrap


def SortRenameHeaders(input, basename, output, minlen=0, simplify=False):
    Seqs = []
    with open(input, "r") as infile:
        for header, sequence in SimpleFastaParser(infile):
            Seqs.append((header, len(sequence), sequence))
    # sort by length
    sortedSeqs = sorted(Seqs, key=lambda x: x[1], reverse=True)
    # loop through and return contigs and keepers
    counter = 1
    with open(output, "w") as outfile:
        for name, length, seq in sortedSeqs:
            if simplify:  # try to just split at first space
                if " " in name:
                    newName = name.split(" ")[0]
                else:
                    newName = name
            else:
                newName = f"{basename}_{counter}"
            if len(newName) > 16:
                print(
                    f"Error. {newName} fasta header too long.",
                    "Choose a different --base name.",
                    "NCBI/GenBank max is 16 characters.",
                )
                raise SystemExit(1)
            if minlen > 0:
                if length >= minlen:
                    # ony write if length
                    outfile.write(">{:}\n{:}\n".format(newName, softwrap(seq)))
            else:
                # always write if we aren't filtering by length
                outfile.write(">{:}\n{:}\n".format(newName, softwrap(seq)))
            counter += 1


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)

    parser = argparse.ArgumentParser(
        prog="sort_rename.py",
        usage="%(prog)s [options] -i genome.fa -o sorted.fa",
        description="Script that sorts input by length and then renames contig headers.",
        epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
        formatter_class=MyFormatter,
    )
    parser.add_argument("-i", "--input", required=True, help="Multi-fasta genome file")
    parser.add_argument("-o", "--out", required=True, help="Cleaned output (FASTA)")
    parser.add_argument(
        "-b", "--base", default="scaffold", help="Basename of contig header"
    )
    parser.add_argument(
        "-s",
        "--simplify",
        action="store_true",
        help="Try to simplify headers, split at first space",
    )
    parser.add_argument(
        "-m", "--minlen", type=int, default=0, help="Contigs shorter than threshold are discarded"
    )
    args = parser.parse_args(args)

    print(("{:,} contigs records loaded".format(countfasta(args.input))))
    print("Sorting and renaming contig headers")
    if args.minlen:
        print(("Removing contigs less than {:} bp".format(args.minlen)))
    SortRenameHeaders(
        args.input, args.base, args.out, minlen=args.minlen, simplify=args.simplify
    )
    print(("{:,} contigs saved to file".format(countfasta(args.out))))


if __name__ == "__main__":
    main(sys.argv[1:])
