#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import uuid
import argparse
import shutil
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
import funannotate.library as lib


def runTanTan(input, output):
    # this is the simplest masking solution, although certainly not ideal
    cmd = ['tantan', '-c', input]
    lib.runSubprocess2(cmd, '.', lib.log, output)


def RepeatModelMask(input, cpus, tmpdir, output, repeatlib, debug):
    lib.log.info("Loading sequences and soft-masking genome")
    outdir = os.path.join(tmpdir, 'RepeatModeler')
    input = os.path.abspath(input)
    output = os.path.abspath(output)
    # lets run RepeatModeler here to get repeat library
    if os.path.exists(outdir):
        shutil.rmtree(outdir)
    os.makedirs(outdir)
    lib.log.info("Soft-masking: building RepeatModeler database")
    with open(debug, 'a') as debug_log:
        subprocess.call(['BuildDatabase', '-name', 'Repeats', input],
                        cwd=outdir, stdout=debug_log, stderr=debug_log)
    lib.log.info("Soft-masking: generating repeat library using RepeatModeler")
    with open(debug, 'a') as debug_log:
        subprocess.call(['RepeatModeler', '-e', 'ncbi', '-database', 'Repeats', '-pa', str(cpus)],
                        cwd=outdir, stdout=debug_log, stderr=debug_log)
    # find name of folder
    RP_folder = '.'
    for i in os.listdir(outdir):
        if i.startswith('RM_'):
            RP_folder = i
    library = os.path.abspath(repeatlib)
    if lib.checkannotations(os.path.join(outdir, RP_folder, 'consensi.fa.classified')):
        shutil.copyfile(os.path.join(outdir, RP_folder,
                                     'consensi.fa.classified'), library)
    # now soft-mask the genome for gene predictors
    outdir2 = os.path.join(tmpdir, 'RepeatMasker')
    if os.path.isdir(outdir2):
        shutil.rmtree(outdir2)
    os.makedirs(outdir2)
    if not os.path.isfile(library):
        lib.log.info(
            "Soft-masking: running RepeatMasker with default library (RepeatModeler found 0 models)")
        with open(debug, 'a') as debug_log:
            subprocess.call(['RepeatMasker', '-e', 'ncbi', '-gff', '-species', 'fungi', '-pa', str(cpus), '-xsmall', '-dir', '.', input],
                            cwd=outdir2, stdout=debug_log, stderr=debug_log)
    else:
        lib.log.info("Soft-masking: running RepeatMasker with custom library")
        with open(debug, 'a') as debug_log:
            subprocess.call(['RepeatMasker', '-e', 'ncbi', '-gff', '-lib', library, '-pa', str(cpus), '-xsmall', '-dir', '.', input],
                            cwd=outdir2, stdout=debug_log, stderr=debug_log)
    for file in os.listdir(outdir2):
        if file.endswith('.masked'):
            shutil.copyfile(os.path.join(outdir2, file), output)


def RepeatMask(input, library, cpus, tmpdir, output, debug):
    input = os.path.abspath(input)
    output = os.path.abspath(output)
    outdir = os.path.join(tmpdir, 'RepeatMasker')
    # now soft-mask the genome for gene predictors
    lib.log.info("Soft-masking: running RepeatMasker with custom library")
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    with open(debug, 'a') as debug_log:
        subprocess.call(['RepeatMasker', '-e', 'ncbi', '-lib', os.path.abspath(library), '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input],
                        stderr=debug_log, stdout=debug_log, cwd=tmpdir)
    for file in os.listdir(outdir):
        if file.endswith('.masked'):
            os.rename(os.path.join(outdir, file), output)


def RepeatMaskSpecies(input, species, cpus, tmpdir, output, debug):
    input = os.path.abspath(input)
    output = os.path.abspath(output)
    outdir = os.path.join(tmpdir, 'RepeatMasker')
    # now soft-mask the genome for gene predictors
    lib.log.info(
        "Soft-masking: running RepeatMasker using %s species" % species)
    if not os.path.isdir(outdir):
        os.makedirs(outdir)
    with open(debug, 'a') as debug_log:
        subprocess.call(['RepeatMasker', '-e', 'ncbi', '-species', species, '-pa', str(cpus), '-xsmall', '-dir', 'RepeatMasker', input],
                        stderr=debug_log, stdout=debug_log, cwd=tmpdir)
    for file in os.listdir(outdir):
        if file.endswith('.masked'):
            os.rename(os.path.join(outdir, file), output)


def main(args):
    # setup menu with argparse
    class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
        def __init__(self, prog):
            super(MyFormatter, self).__init__(prog, max_help_position=48)
    parser = argparse.ArgumentParser(prog='funannotate-mask.py',
                                     description='''Wrapper for RepeatModeler/RepeatMasker''',
                                     epilog="""Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
                                     formatter_class=MyFormatter)
    parser.add_argument('-i', '--input', required=True,
                        help='genome assembly FASTA format')
    parser.add_argument('-o', '--out', required=True,
                        help='Output softmasked FASTA file')
    parser.add_argument('--debug', action='store_true',
                        help='Keep intermediate files')
    parser.add_argument('-m', '--method', default='tantan', choices=[
                        'repeatmodeler', 'repeatmasker', 'tantan'], help='Method to mask repeats with')
    parser.add_argument('-s', '--repeatmasker_species',
                        help='RepeatMasker species, will skip repeatmodeler')
    parser.add_argument('-l', '--repeatmodeler_lib',
                        help='Pre-computed RepeatModeler (or other) repetitive elements')
    parser.add_argument('--cpus', default=2, type=int,
                        help='Number of CPUs to use')
    args = parser.parse_args(args)

    # create log file for Repeats(capture stderr)
    log_name = 'funannotate-mask.log'
    if os.path.isfile(log_name):
        os.remove(log_name)

    # initialize script, log system info and cmd issue at runtime
    lib.setupLogging(log_name)
    cmd_args = " ".join(sys.argv)+'\n'
    lib.log.debug(cmd_args)
    print("-------------------------------------------------------")
    lib.SystemInfo()

    # get version of funannotate
    version = lib.get_version()
    lib.log.info("Running funanotate v{:}".format(version))

    repeats = None
    tmpdir = None
    if args.method == 'tantan':
        programs = ['tantan']
        lib.CheckDependencies(programs)
        lib.log.info('Soft-masking simple repeats with tantan')
        runTanTan(args.input, args.out)
    else:
        programs = ['RepeatMasker']
        if args.method == 'repeatmodeler':
            programs += ['BuildDatabase', 'RepeatModeler']
        lib.CheckDependencies(programs)

        # create tmpdir
        pid = uuid.uuid4()
        tmpdir = 'mask_'+str(pid)
        os.makedirs(tmpdir)

        # parse options which dictates how repeatmodeler/masker are run
        if not args.repeatmodeler_lib:  # no fasta file given, so
            if not args.repeatmasker_species:  # no species given, so run entire repeatmodler + repeat masker
                repeats = 'repeatmodeler-library.'+str(pid)+'.fasta'
                RepeatModelMask(args.input, args.cpus, tmpdir,
                                args.out, repeats, log_name)
            else:
                RepeatMaskSpecies(
                    args.input, args.repeatmasker_species, args.cpus, tmpdir, args.out, log_name)
        else:
            if lib.checkannotations(args.repeatmodeler_lib):
                RepeatMask(args.input, args.repeatmodeler_lib,
                           args.cpus, tmpdir, args.out, log_name)
            else:
                lib.log.error('ERROR: repeat library is not a valid file: {:}'.format(
                    args.repeatmodeler_lib))
                sys.exit(1)

    # output some stats on %reads masked.
    scaffolds = 0
    maskedSize = 0
    GenomeLength = 0
    with open(args.out, 'r') as input:
        for rec, Seq in SimpleFastaParser(input):
            scaffolds += 1
            GenomeLength += len(Seq)
            maskedSize += lib.n_lower_chars(Seq)

    percentMask = maskedSize / float(GenomeLength)
    lib.log.info('Repeat soft-masking finished: \nMasked genome: {:}\nnum scaffolds: {:,}\nassembly size: {:,} bp\nmasked repeats: {:,} bp ({:.2f}%)'.format(
        os.path.abspath(args.out), scaffolds, GenomeLength, maskedSize, percentMask*100))
    if repeats:
        lib.log.info('RepeatModeler library: {:}'.format(repeats))
    # clean up
    if not args.debug:
        if tmpdir:
            lib.SafeRemove(tmpdir)
    print("-------------------------------------------------------")


if __name__ == "__main__":
    main(sys.argv[1:])
