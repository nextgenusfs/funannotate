#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import uuid
import argparse
import fnmatch
import subprocess
import funannotate.library as lib


def worker(input):
    logfile = os.path.join(tmpdir, 'Trinity-gg.log')
    cmd = input.split(' ')
    # make sure no empty items
    cmd = [x for x in cmd if x]
    with open(logfile, 'a') as output:
        subprocess.call(cmd, cwd=os.path.join(
            tmpdir, 'trinity_gg'), stdout=output, stderr=output)


def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        worker(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def find_files(directory, pattern):
    for root, dirs, files in os.walk(directory):
        for basename in files:
            if fnmatch.fnmatch(basename, pattern):
                filename = os.path.join(root, basename)
                yield filename


def runTrinityGG(genome, readTuple, longReads, shortBAM, output, args=False):
    '''
    function will run genome guided Trinity. First step will be to run hisat2 to align reads
    to the genome, then pass that BAM file to Trinity to generate assemblies
    '''
    if not lib.checkannotations(shortBAM):
        # build hisat2 index, using exons and splice sites
        lib.log.info("Building Hisat2 genome index")
        cmd = ['hisat2-build', '-p',
               str(args.cpus), genome, os.path.join(tmpdir, 'hisat2.genome')]
        lib.runSubprocess4(cmd, '.', lib.log)
        # align reads using hisat2
        lib.log.info("Aligning reads to genome using Hisat2")
        # use bash wrapper for samtools piping for SAM -> BAM -> sortedBAM
        # use half number of threads for bam compression threads
        bamthreads = (args.cpus + 2 // 2) // 2
        if args.stranded != 'no' and not readTuple[2]:
            hisat2cmd = ['hisat2', '-p', str(args.cpus), '--max-intronlen',
                         str(args.max_intronlen), '--dta',
                         '-x', os.path.join(tmpdir, 'hisat2.genome'),
                         '--rna-strandness', args.stranded]
        else:
            hisat2cmd = ['hisat2', '-p', str(args.cpus), '--max-intronlen',
                         str(args.max_intronlen), '--dta',
                         '-x', os.path.join(tmpdir, 'hisat2.genome')]
        if readTuple[0] and readTuple[1]:
            hisat2cmd = hisat2cmd + ['-1', readTuple[0], '-2', readTuple[1]]
        if readTuple[2]:
            hisat2cmd = hisat2cmd + ['-U', readTuple[2]]
        samtools_cmd = ['samtools', 'sort', '--reference', genome,
                        '-@', str(bamthreads), '-o', shortBAM, '-']
        lib.log.debug('{} | {}'.format(' '.join(hisat2cmd), ' '. join(samtools_cmd)))
        p1 = subprocess.Popen(hisat2cmd, stdout=subprocess.PIPE, stderr=FNULL)
        p2 = subprocess.Popen(samtools_cmd, stdout=subprocess.PIPE, stderr=FNULL, stdin=p1.stdout)
        p1.stdout.close()
        p2.communicate()
    else:
        lib.log.info('Existig Hisat2 alignments found: {:}'.format(shortBAM))

    # now launch Trinity genome guided
    TrinityLog = os.path.join(tmpdir, 'Trinity-gg.log')
    lib.log.info("Running genome-guided Trinity, logfile: %s" % TrinityLog)
    lib.log.info(
        "Clustering of reads from BAM and preparing assembly commands")
    jaccard_clip = []
    if args.jaccard_clip:
        jaccard_clip = ['--jaccard_clip']
    if args.stranded != 'no':
        cmd = ['Trinity', '--SS_lib_type', args.stranded, '--no_distributed_trinity_exec',
               '--genome_guided_bam', shortBAM, '--genome_guided_max_intron', str(
                   args.max_intronlen),
               '--CPU', str(args.cpus), '--max_memory', args.memory, '--output', os.path.join(tmpdir, 'trinity_gg')]
    else:
        cmd = ['Trinity', '--no_distributed_trinity_exec', '--genome_guided_bam', shortBAM,
               '--genome_guided_max_intron', str(
                   args.max_intronlen), '--CPU', str(args.cpus),
               '--max_memory', args.memory, '--output', os.path.join(tmpdir, 'trinity_gg')]
    cmd = cmd + jaccard_clip
    if longReads and lib.checkannotations(longReads):
        cmd = cmd + ['--long_reads', os.path.realpath(longReads)]
    lib.runSubprocess2(cmd, '.', lib.log, TrinityLog)
    commands = os.path.join(tmpdir, 'trinity_gg', 'trinity_GG.cmds')

    # this will create all the Trinity commands, will now run these in parallel using multiprocessing
    # in Python (seems to be much faster than Parafly on my system)
    file_list = []
    with open(commands, 'r') as cmdFile:
        for line in cmdFile:
            line = line.replace('\n', '')
            # don't think this should be appended to every command....
            line = line.replace('--no_distributed_trinity_exec', '')
            line = line.replace('"', '')  # don't need these double quotes
            file_list.append(line)
    lib.log.info("Assembling "+"{0:,}".format(len(file_list)) +
                 " Trinity clusters using %i CPUs" % (args.cpus-1))
    lib.runMultiProgress(safe_run, file_list, args.cpus-1, progress=args.progress)

    # collected output files and clean
    outputfiles = os.path.join(
        tmpdir, 'trinity_gg', 'trinity_output_files.txt')
    with open(outputfiles, 'w') as fileout:
        for filename in find_files(os.path.join(tmpdir, 'trinity_gg'), '*inity.fasta'):
            fileout.write('%s\n' % filename)
    # now grab them all using Trinity script
    cmd = ['perl', os.path.abspath(os.path.join(
        TRINITY, 'util', 'support_scripts', 'GG_partitioned_trinity_aggregator.pl')), 'Trinity_GG']
    lib.runSubprocess5(cmd, '.', lib.log, outputfiles, output)
    lib.log.info('{:,} transcripts derived from Trinity'.format(
        lib.countfasta(output)))


# setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)


parser = argparse.ArgumentParser(prog='trinity.py',
                                 description='''Script is a wrapper for multiprocessing trinity.''',
                                 epilog="""Written by Jon Palmer (2017-2018) nextgenusfs@gmail.com""",
                                 formatter_class=MyFormatter)
parser.add_argument('-f', '--fasta', required=True,
                    help='Genome in FASTA format')
parser.add_argument('-l', '--left', help='Left (R1) FASTQ Reads')
parser.add_argument('-r', '--right', help='Right (R2) FASTQ Reads')
parser.add_argument('-s', '--single', help='Single ended FASTQ Reads')
parser.add_argument('-t', '--tmpdir', help='temoporary directory')
parser.add_argument('--long', help='Long reads')
parser.add_argument('-b', '--bam', help='BAM file')
parser.add_argument('--logfile', help='logfile')
parser.add_argument('-o', '--out', required=True, help='Trinity transcripts')
parser.add_argument('--memory', default='50G',
                    help='RAM to use for Jellyfish/Trinity')
parser.add_argument('--stranded', default='no',
                    choices=['RF', 'FR', 'F', 'R', 'no'], help='RNA seq strandedness')
parser.add_argument('--jaccard_clip', action='store_true',
                    help='Turn on jaccard_clip for dense genomes')
parser.add_argument('--max_intronlen', default=3000,
                    help='Maximum intron length for gene models')
parser.add_argument('--cpus', default=2, type=int,
                    help='Number of CPUs to use')
parser.add_argument('--TRINITYHOME', '--TRINITY_HOME', dest='TRINITYHOME', required=True,
                    help='Path to Trinity config directory, $TRINITYHOME')
parser.add_argument('--no-progress', dest='progress', action='store_false',
                    help='no progress on multiprocessing')
args = parser.parse_args()

# start logic here
global tmpdir, parentdir, TRINITY, FNULL
TRINITY = args.TRINITYHOME
parentdir = os.path.join(os.path.dirname(__file__))
if args.tmpdir:
    tmpdir = args.tmpdir
else:
    tmpdir = 'trinity_GG_'+str(uuid.uuid4())
if args.logfile:
    log_name = args.logfile
else:
    log_name = 'funannotate-trinity.log'

if os.path.isfile(log_name):
    os.remove(log_name)

if args.left and args.right:
    all_reads = (os.path.abspath(args.left), os.path.abspath(args.right), None)
elif args.single:
    all_reads = (None, None, os.path.abspath(args.single))
else:
    print('Error: specify either PE reads via -l and -r or SE reads via -s')
    sys.exit(1)
# initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

# run trinity
runTrinityGG(args.fasta, all_reads, args.long, args.bam, args.out, args=args)
