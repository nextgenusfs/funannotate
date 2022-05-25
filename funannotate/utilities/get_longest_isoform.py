#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys, re, os, gzip, argparse
from Bio import SeqIO

def main(inargs):
        # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='get_longest_isoform',
                                     description='''Script to extract longest isoform of protein or transcript file from funannotate or where gene is tagged in header.''',
                                     epilog="""Written by Jason Stajich (2022) @hyphaltip""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help='fasta formatted transcript or protein file')
    parser.add_argument('-o', '--output', help='Output basename')

    parser.add_argument('-v', '--verbose', help='Extra verbose output',dest='verbose', default=False, action='store_true')

    args = parser.parse_args(inargs)
    genes = {}
    if not args.output:
        args.output = args.input + ".longest"
    transmatch = re.compile(r'\-T\d+$')
    genematch  = re.compile(r'gene[:=](\S+)')
    recCount = 0
    handle = args.input
    if args.input.endswith('.gz'):
        handle =  gzip.open(args.input,"rt")
    for rec in SeqIO.parse(handle, "fasta"):
        id = rec.id
        description = rec.description
        geneid = id
        m = transmatch.search(id)
        if m:
            geneid = description.split()[1]
        else:
            m = genematch.search(description)
            if m:
                geneid = m.group(1)
        if geneid == id:
            if args.verbose:
                print("Warning: could not parse gene name from header '{}' '{}'".format(id,description))
        if geneid not in genes or len(rec) > len(genes[geneid]):
            genes[geneid] = rec
        recCount += 1

    print("{} genes and {} total sequences (isoforms) seen".format(len(genes),recCount))
    SeqIO.write(genes.values(),args.output,'fasta')

if __name__ == "__main__":
    main(sys.argv[1:])
