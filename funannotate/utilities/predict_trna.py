#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import funannotate.library as lib


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(
        prog='predict_trna.py',
        description='''Run tRNAscan-SE on a genome FASTA file and produce a
        non-overlapping tRNA GFF3 file, the same trnascan.no-overlaps.gff3
        style result produced in predict_misc by funannotate predict.''',
        epilog="""Written by Jon Palmer / Jason Stajich nextgenusfs@gmail.com""",
        formatter_class=MyFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help='Genome FASTA file (should be soft-masked)')
    parser.add_argument('-o', '--out', required=True,
                        help='Output non-overlapping tRNA GFF3 file (trnascan.no-overlaps.gff3)')
    parser.add_argument('--genes',
                        help='Existing gene models GFF3 to drop overlapping tRNA models. Default: none (keep all)')
    parser.add_argument('--gaps',
                        help='Assembly gaps GFF3/BED to drop overlapping tRNA models. Default: none')
    parser.add_argument('--trnascan',
                        help='Pre-computed tRNAscan-SE tabular output, skips running tRNAscan-SE')
    parser.add_argument('--cpus', default=1, type=int,
                        help='Number of CPUs to use for tRNAscan-SE')
    parser.add_argument('--tmpdir',
                        help='Directory to write intermediate files. Default: directory of --out')
    parser.add_argument('--keep_tmp', action='store_true',
                        help='Keep intermediate/temporary files')
    parser.add_argument('--logfile',
                        help='Logfile output file. Default: predict-tRNA.log in --tmpdir')
    args = parser.parse_args(args)

    args.input = os.path.abspath(args.input)
    args.out = os.path.abspath(args.out)
    outdir = os.path.dirname(args.out) or '.'
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    tmpdir = os.path.abspath(args.tmpdir) if args.tmpdir else outdir
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)

    log_name = args.logfile if args.logfile else os.path.join(tmpdir, 'predict-tRNA.log')
    if os.path.isfile(log_name):
        os.remove(log_name)
    lib.setupLogging(log_name)

    if not lib.checkannotations(args.input):
        lib.log.error('Input genome FASTA not found: {}'.format(args.input))
        sys.exit(1)

    # Step 1: run tRNAscan-SE (or copy pre-calculated results), length-filter,
    # and convert to GFF3 -- same as lib.runtRNAscan() used in predict.py
    tRNAscanGFF3 = os.path.join(tmpdir, 'trnascan.gff3')
    if args.trnascan:
        lib.log.info('Using pre-computed tRNAscan-SE results: {}'.format(args.trnascan))
    else:
        lib.log.info('Predicting tRNAs with tRNAscan-SE')
    result = lib.runtRNAscan(
        args.input, tmpdir, tRNAscanGFF3, cpus=args.cpus, precalc=args.trnascan)

    # lib.runtRNAscan() returns False both when tRNAscan-SE errors out AND
    # when it runs fine but finds zero tRNAs (empty output) -- the latter is
    # a normal result, not a failure, so match predict.py's behavior of
    # writing an empty valid GFF3 rather than treating it as fatal.
    if not result:
        lib.log.info('No tRNAs found (or tRNAscan-SE failed, check logfile: {})'.format(log_name))
        with open(args.out, 'w') as outfile:
            outfile.write('##gff-version 3\n')
        return

    # Step 2: drop tRNA models that overlap existing gene models/assembly gaps
    # -- same as lib.validate_tRNA() used in predict.py. validate_tRNA always
    # requires a "genes" GFF3 to intersect against, so fall back to an empty
    # placeholder (keeps all tRNA models) when --genes isn't provided.
    genesGFF3 = args.genes
    placeholderGenes = False
    if not genesGFF3:
        genesGFF3 = os.path.join(tmpdir, 'no-genes.placeholder.gff3')
        with open(genesGFF3, 'w') as outfile:
            outfile.write('##gff-version 3\n')
        placeholderGenes = True

    lib.validate_tRNA(tRNAscanGFF3, genesGFF3, args.gaps, args.out)
    lib.log.info(
        '{:,} tRNAscan models are valid (non-overlapping)'.format(
            lib.countGFFgenes(args.out)))

    if placeholderGenes and not args.keep_tmp:
        os.remove(genesGFF3)


if __name__ == "__main__":
    main(sys.argv[1:])
