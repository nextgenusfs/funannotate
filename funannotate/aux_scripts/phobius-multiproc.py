#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import uuid
import time
import multiprocessing
import argparse
import shutil
import funannotate.library as lib

# setup menu with argparse


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)


parser = argparse.ArgumentParser(
    prog='phobius-multiproc.py',
    usage="%(prog)s [options] -i proteome.fasta",
    description='''Script that runs phobius remotely.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)
parser.add_argument('-i', '--input', required=True, help='whole proteome')
parser.add_argument('-o', '--out', required=True, help='Phobius results')
parser.add_argument('-e', '--email', help='Email address for IPRSCAN server')
parser.add_argument('-l', '--logfile',
                    default='phobius-multiproc.log', help='Logfile')
parser.add_argument('--debug', action='store_true',
                    help='Keep intermediate files')
args = parser.parse_args()


def runPhobiusRemote(Input):
    base = Input.split('/')[-1]
    base = base.split('.fa')[0]
    OUTPATH = os.path.join(TMPDIR, base)
    cmd = ['perl', os.path.join(parentdir, 'phobius-remote.pl'),
           '--email', args.email, '-f', 'short', '--outfile', base, Input]
    lib.runSubprocess(cmd, TMPDIR, lib.log)
    time.sleep(1)  # make sure there is time for all files to show up
    os.rename(OUTPATH+'.out.txt', OUTPATH+'.phobius')
    os.remove(OUTPATH+'.sequence.txt')


def runPhobiusLocal(Input):
    base = Input.split('/')[-1]
    base = base.split('.fa')[0]
    OUTPATH = os.path.join(TMPDIR, base+'.phobius')
    cmd = ['phobius.pl', '-short', Input]
    lib.runSubprocess2(cmd, TMPDIR, lib.log, OUTPATH)


global parentdir
parentdir = os.path.join(os.path.dirname(__file__))

# create log file
log_name = args.logfile
if os.path.isfile(log_name):
    os.remove(log_name)

# initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

# create tmpdir to store fasta files and output files
TMPDIR = 'phobius_' + str(uuid.uuid4())

# split fasta
lib.splitFASTA(args.input, TMPDIR)

# now get list of files in tmpdir
proteins = []
for file in os.listdir(TMPDIR):
    if file.endswith('.fa'):
        proteins.append(file)

# now run the script
if lib.which('phobius.pl'):
    lib.runMultiProgress(runPhobiusLocal, proteins,
                         multiprocessing.cpu_count())
else:
    lib.runMultiProgress(runPhobiusRemote, proteins,
                         29)  # max is 30 jobs at a time

# collect all results
phobius = []
for file in os.listdir(TMPDIR):
    if file.endswith('.phobius'):
        phobius.append(os.path.join(TMPDIR, file))

# write output
TMdomain = 0
SigPep = 0
with open(args.out, 'w') as output:
    output.write("%s\t%s\t%s\t%s\n" % ('ID', 'TM', 'SP', 'Prediction'))
    for x in phobius:
        with open(x, 'r') as input:
            line = input.readlines()
            try:
                result = line[1].split(' ')
                result = [x for x in result if x]
                if result[1] == 'prediction':
                    continue
                if int(result[1]) > 0:
                    TMdomain += 1
                if result[2] == 'Y':
                    SigPep += 1
                output.write("%s\t%s\t%s\t%s\n" % (
                    result[0], result[1],
                    result[2], result[3].replace('\n', '')))
            except IndexError:
                pass

# clean
if not args.debug:
    shutil.rmtree(TMPDIR)
lib.log.debug("%i total proteins, %i TMdomain, %i Signal Peptide" %
              (len(phobius), TMdomain, SigPep))
