#!/usr/bin/env python

import sys
import os
import argparse
import shutil
import inspect
from Bio.SeqIO.FastaIO import SimpleFastaParser
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='funannotate-mask.py',
    description = '''Wrapper for RepeatModeler/RepeatMasker''',
    epilog = """Written by Jon Palmer (2018) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--input', required=True, help='genome assembly FASTA format')
parser.add_argument('-o', '--out', required=True, help='Output softmasked FASTA file')
parser.add_argument('--debug', action='store_true', help='Keep intermediate files')
parser.add_argument('-s', '--repeatmasker_species', help='RepeatMasker species, will skip repeatmodeler')
parser.add_argument('-l', '--repeatmodeler_lib', help='Pre-computed RepeatModeler (or other) repetitive elements')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to use')
args=parser.parse_args()
         
#create log file for Repeats(capture stderr)
log_name = 'funannotate-mask.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print("-------------------------------------------------------")
lib.SystemInfo()

#get version of funannotate
version = lib.get_version()
lib.log.info("Running %s" % version)

programs = ['RepeatMasker']
if not args.repeatmodeler_lib or not args.repeatmasker_species:
    programs = programs + ['BuildDatabase', 'RepeatModeler']
lib.CheckDependencies(programs)

#create tmpdir
pid = os.getpid()
tmpdir = 'mask_'+str(pid)
os.makedirs(tmpdir)
repeats = None
#parse options which dictates how repeatmodeler/masker are run
if not args.repeatmodeler_lib: #no fasta file given, so
    if not args.repeatmasker_species: #no species given, so run entire repeatmodler + repeat masker
        repeats = 'repeatmodeler-library.'+str(pid)+'.fasta'
        lib.RepeatModelMask(args.input, args.cpus, tmpdir, args.out, repeats, log_name)
    else:
        lib.RepeatMaskSpecies(args.input, args.repeatmasker_species, args.cpus, tmpdir, args.out, log_name)
else:
    if lib.checkannotations(args.repeatmodeler_lib):
        lib.RepeatMask(args.input, args.repeatmodeler_lib, args.cpus, tmpdir, args.out, log_name)
    else:
        lib.log.error('ERROR: repeat library is not a valid file: {:}'.format(args.repeatmodeler_lib))
        sys.exit(1)

#output some stats on %reads masked.
scaffolds = 0
maskedSize = 0
GenomeLength = 0
with open(args.out, 'rU') as input:
    for rec, Seq in SimpleFastaParser(input):
        scaffolds += 1
        GenomeLength += len(Seq)
        maskedSize += lib.n_lower_chars(Seq)

percentMask = maskedSize / float(GenomeLength)
lib.log.info('Repeatmasking finished: \nMasked genome: {:}\nnum scaffolds: {:,}\nassembly size: {:,} bp\nmasked repeats: {:,} bp ({:.2f}%)'.format(os.path.abspath(args.out), scaffolds, GenomeLength, maskedSize, percentMask*100))
if repeats:
    lib.log.info('RepeatModeler library: {:}'.format(repeats))
#clean up
if not args.debug:
    lib.SafeRemove(tmpdir)
print("-------------------------------------------------------")