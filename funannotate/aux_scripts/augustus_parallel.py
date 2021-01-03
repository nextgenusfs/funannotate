#!/usr/bin/env python

import sys
import subprocess
import os
import uuid
import shutil
import argparse
from Bio import SeqIO
import funannotate.library as lib

# setup menu with argparse


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)


parser = argparse.ArgumentParser(prog='augustus_parallel.py',
                                 usage="%(prog)s [options] -i genome.fasta -s botrytis_cinera -o prediction_output_base",
                                 description='''Script runs augustus in parallel to use multiple processors''',
                                 epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
                                 formatter_class=MyFormatter)
parser.add_argument('-i', '--input', required=True,
                    help='Genome in FASTA format')
parser.add_argument('-o', '--out', required=True,
                    help='Basename of output files')
parser.add_argument('-s', '--species', required=True,
                    help='Augustus species name')
parser.add_argument('--hints', help='Hints file (PE)')
parser.add_argument('--cpus', default=2, type=int,
                    help='Number of CPUs to run')
parser.add_argument('-v', '--debug', action='store_true',
                    help='Keep intermediate files')
parser.add_argument('--logfile', default='augustus-parallel.log',
                    help='logfile')
parser.add_argument('--local_augustus')
parser.add_argument('--AUGUSTUS_CONFIG_PATH')
parser.add_argument('-e', '--extrinsic', help='augustus extrinsic file')
parser.add_argument('--no-progress', dest='progress', action='store_false',
                    help='no progress on multiprocessing')
args = parser.parse_args()

# check for augustus installation
if args.AUGUSTUS_CONFIG_PATH:
    AUGUSTUS = args.AUGUSTUS_CONFIG_PATH
else:
    try:
        AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
    except KeyError:
        print("$AUGUSTUS_CONFIG_PATH environmental variable not found, Augustus is not properly configured")
        sys.exit(1)

if AUGUSTUS.endswith('config'):
    AUGUSTUS_BASE = AUGUSTUS.replace('config', '')
elif AUGUSTUS.endswith('config'+os.sep):
    AUGUSTUS_BASE = AUGUSTUS.replace('config'+os.sep, '')
else:
    AUGUSTUS_BASE = AUGUSTUS

# see if local species passed
if args.local_augustus:
    LOCALAUGUSTUS = args.local_augustus
else:
    LOCALAUGUSTUS = AUGUSTUS

# setup hints and extrinic input, hard coded for protein and transcript alignments from funannotate
extrinsic = '--extrinsicCfgFile={:}'.format(args.extrinsic)


def countGFFgenes(input):
    count = 0
    with open(input, 'r') as f:
        for line in f:
            if "\tgene\t" in line:
                count += 1
    return count


def runAugustus(Input):
    if '_part' in Input:
        chr = Input.split('_part')[0]
    else:
        chr = Input
    species = '--species='+args.species
    hints_input = '--hintsfile='+args.hints
    aug_out = os.path.join(tmpdir, Input+'.augustus.gff3')
    core_cmd = ['augustus', species, '--AUGUSTUS_CONFIG_PATH={:}'.format(LOCALAUGUSTUS), '--softmasking=1',
                '--gff3=on', '--UTR=off', '--stopCodonExcludedFromCDS=False', os.path.join(tmpdir, chr+'.fa')]
    if args.hints:
        core_cmd.insert(2, extrinsic)
        core_cmd.insert(3, hints_input)
    if Input in ranges:
        start = ranges.get(Input)[0]
        end = ranges.get(Input)[1]
        core_cmd.insert(2, '--predictionStart='+str(start))
        core_cmd.insert(3, '--predictionEnd='+str(end))
    # try using library module
    lib.runSubprocess2(core_cmd, '.', lib.log, aug_out)


log_name = args.logfile
if os.path.isfile(log_name):
    os.remove(log_name)

# initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

lib.log.debug('AUGUSTUS_CONFIG_PATH={:}'.format(AUGUSTUS))
lib.log.debug('Augustus Base directory={:}'.format(AUGUSTUS_BASE))
lib.log.debug('Local Augustus path={:}'.format(LOCALAUGUSTUS))

# first step is to split input fasta file into individual files in tmp folder
lib.log.debug("Splitting contigs and hints files")
tmpdir = 'augustus_tmp_'+str(uuid.uuid4())
os.makedirs(tmpdir)
scaffolds = []
global ranges
ranges = {}
with open(args.input, 'r') as InputFasta:
    for record in SeqIO.parse(InputFasta, 'fasta'):
        contiglength = len(record.seq)
        if contiglength > 500000:  # split large contigs
            num_parts = contiglength / 500000 + 1
            chunks = contiglength / num_parts
            for i in range(0, int(num_parts)):
                name = str(record.id)+'_part'+str(i+1)
                scaffolds.append(name)
                outputfile = os.path.join(tmpdir, str(record.id)+'.fa')
                if i == 0:  # this is first record
                    start = 1
                    end = chunks + 10000
                else:
                    start = end - 10000
                    end = start + chunks + 10000
                if end > contiglength:
                    end = contiglength
                if not name in ranges:
                    ranges[name] = (start, end)
                with open(outputfile, 'w') as output:
                    SeqIO.write(record, output, 'fasta')
        else:
            name = str(record.id)
            scaffolds.append(name)
            outputfile = os.path.join(tmpdir, name+'.fa')
            with open(outputfile, 'w') as output:
                SeqIO.write(record, output, 'fasta')

# now loop through each scaffold running augustus
if args.cpus > len(scaffolds):
    num = len(scaffolds)
else:
    num = args.cpus
lib.log.debug("Running Augustus on %i chunks, using %i CPUs" %
              (len(scaffolds), num))
lib.runMultiProgress(runAugustus, scaffolds, num, progress=args.progress)


lib.log.debug("Augustus prediction is finished, now concatenating results")
with open(os.path.join(tmpdir, 'augustus_all.gff3'), 'w') as output:
    for file in scaffolds:
        file = os.path.join(tmpdir, file+'.augustus.gff3')
        with open(file) as input:
            output.write(input.read())

if lib.checkannotations(os.path.join(tmpdir, 'augustus_all.gff3')):
    lib.log.debug('Augustus finished, now joining results')
if lib.which_path('join_aug_pred.pl'):
    join_script = 'join_aug_pred.pl'
else:
    join_script = os.path.join(AUGUSTUS_BASE, 'scripts', 'join_aug_pred.pl')

cmd = '{:} < {:} > {:}'.format(join_script, os.path.join(
    tmpdir, 'augustus_all.gff3'), args.out)
lib.log.debug(cmd)

with open(args.out, 'w') as finalout:
    with open(os.path.join(tmpdir, 'augustus_all.gff3'), 'r') as infile:
        subprocess.call([join_script], stdin=infile, stdout=finalout)

if not args.debug:
    shutil.rmtree(tmpdir)
lib.log.info('{:,} predictions from Augustus'.format(countGFFgenes(args.out)))
