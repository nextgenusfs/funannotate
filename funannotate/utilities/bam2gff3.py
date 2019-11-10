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
    parser = argparse.ArgumentParser(prog='bam2gff3.py',
                                     description='''Script to convert BAM to GFF3.''',
                                     epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--bam', required=True, help='input BAM')
    parser.add_argument('-o', '--output', required=True, help='Output GFF3')
    args = parser.parse_args(args)

    # convert BAM to gff3
    lib.bam2gff3(args.bam, args.output)


if __name__ == "__main__":
    main(sys.argv[1:])
