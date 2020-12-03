#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import multiprocessing
import shutil
import argparse
import subprocess
from Bio.SeqIO.FastaIO import SimpleFastaParser
import funannotate.library as lib


def list_slice(S, step):
    return [S[i::step] for i in range(step)]

# script to run tbl2asn in parallel for larger genomes


def split_tbl2asn(folder):
    '''
    function to chunk the genome and annotation files into parts if > 10,000 contigs to
    conform to NCBI recommendations and avoid the 2GB threshold of sequin files
    '''
    numSeqs = 0
    genomeSize = 0
    with open(os.path.join(folder, 'genome.fsa'), 'r') as fastain:
        for Header, Seq in SimpleFastaParser(fastain):
            numSeqs += 1
            genomeSize += len(Seq)
    # if less than 10,000 contigs and less than 100 MB, then don't split and just run it
    if numSeqs < 10000 and genomeSize < int(100e6):
        # move to subfolder for multiprocessing to work correctly
        if os.path.isdir(os.path.join(folder, '1')):
            lib.SafeRemove(os.path.join(folder, '1'))
        os.makedirs(os.path.join(folder, '1'))
        shutil.copyfile(os.path.join(folder, 'genome.fsa'),
                        os.path.join(folder, '1', 'genome.fsa'))
        shutil.copyfile(os.path.join(folder, 'genome.tbl'),
                        os.path.join(folder, '1', 'genome.tbl'))
    else:
        # rounded_up = -(-numerator // denominator) #nice trick to round up
        if genomeSize > int(100e6):
            chunks = -(-genomeSize // int(100e6))  # split into 100 MB chunks
        else:
            chunks = -(-numSeqs // 10000)
        Records = []
        with open(os.path.join(folder, 'genome.fsa'), 'r') as fastain:
            for tup in SimpleFastaParser(fastain):
                Records.append(tup)
        # sort the fasta tuples by size
        Records = sorted(Records, key=lambda x: len(x[1]), reverse=True)
        # shuffle them into lists like dealing playing cards then all chunks have similar sizes
        sliced_records = list_slice(Records, chunks)
        # loop through and add headers to dictionary for tbl splitting lookup
        headers = {}
        for i, x in enumerate(sliced_records):
            if os.path.isdir(os.path.join(folder, str(i+1))):
                lib.SafeRemove(os.path.join(folder, str(i+1)))
            os.makedirs(os.path.join(folder, str(i+1)))
            with open(os.path.join(folder, str(i+1), 'genome'+str(i+1)+'.fsa'), 'w') as outfile:
                for seq in x:
                    outfile.write('>{:}\n{:}\n'.format(seq[0], seq[1]))
                    headers[seq[0]] = i+1
        # now parse tbl file and split in same way as fasta files
        with open(os.path.join(folder, 'genome.tbl'), 'r') as tblin:
            for contig in lib.readBlocks(tblin, '>Feature'):
                ID = contig[0].split(' ')[-1].rstrip()
                filenum = None
                if ID in headers:
                    filenum = headers.get(ID)
                if filenum:
                    with open(os.path.join(folder, str(filenum), 'genome'+str(filenum)+'.tbl'), 'a') as tblout:
                        tblout.write(''.join(contig))


def tbl2asn_safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        tbl2asn_runner(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def tbl2asn_runner(cmd, dir):
    cmd = cmd + ['-Z', os.path.join(dir, 'discrepency.report.txt'), '-p', dir]
    FNULL = open(os.path.devnull, 'w')
    subprocess.call(cmd, stdout=FNULL, stderr=FNULL)


def runtbl2asn_parallel(folder, template, discrepency, organism, isolate, strain, parameters, version, cpus):
    '''
    function to run NCBI tbl2asn
    '''
    # make sure ouput that will be appended to is not there
    for file in [os.path.join(folder, 'genome.val'), os.path.join(folder, 'errorsummary.val'), os.path.join(folder, 'genome.gbf'), discrepency]:
        lib.SafeRemove(file)
    # get funannotate version
    fun_version = lib.get_version()
    # input should be a folder
    if not os.path.isdir(folder):
        lib.log.error("tbl2asn error: %s is not a directory, exiting" % folder)
        sys.exit(1)
    # based on organism, isolate, strain, construct meta info for -j flag
    if not organism:
        lib.log.error("tbl2asn error: organism not specified")
        sys.exit(1)
    meta = "[organism=" + organism + "]"
    if isolate:
        isolate_meta = "[isolate=" + isolate + "]"
        meta = meta + " " + isolate_meta
    if strain:
        strain_meta = "[strain=" + strain + "]"
        meta = meta + " " + strain_meta
    cmd = ['tbl2asn', '-y', '"Annotated using '+fun_version+'"', '-N', str(version),
           '-t', template, '-M', 'n', '-j', '"'+meta+'"', '-V', 'b',
           '-c', 'f', '-T', '-a', 'r10u']
    # check for custom parameters
    if parameters:
        params = parameters.split(' ')
        cmd = cmd + params
    # check for folders in the input folder, if present, run tbl2asn on each folder and then combine
    multiple = []
    for file in os.listdir(folder):
        if os.path.isdir(os.path.join(folder, file)):
            multiple.append(os.path.join(folder, file))
    if len(multiple) == 0:
        multiple.append(folder)
    p = multiprocessing.Pool(cpus)
    results = []
    for i in multiple:
        results.append(p.apply_async(tbl2asn_safe_run, (cmd, i)))
    p.close()
    p.join()
    # now collect the results make in main folder
    # first delete any of the outputs you might be appending to
    with open(os.path.join(folder, 'genome.val'), 'a') as validation:
        with open(discrepency, 'a') as discrep:
            with open(os.path.join(folder, 'errorsummary.val'), 'a') as summary:
                with open(os.path.join(folder, 'genome.gbf'), 'a') as genbank:
                    for dirName, subdirList, fileList in os.walk(folder, topdown=False):
                        if len(subdirList) > 0:
                            continue
                        for f in fileList:
                            if f == 'errorsummary.val':
                                with open(os.path.join(dirName, f)) as infile:
                                    summary.write(infile.read())
                            elif f.endswith('.val'):
                                with open(os.path.join(dirName, f)) as infile:
                                    validation.write(infile.read())
                            elif f.endswith('.gbf'):
                                with open(os.path.join(dirName, f)) as infile:
                                    genbank.write(infile.read())
                            elif f.endswith('.tbl'):
                                shutil.copyfile(os.path.join(
                                    dirName, f), os.path.join(folder, f))
                            elif f.endswith('.sqn'):
                                shutil.copyfile(os.path.join(
                                    dirName, f), os.path.join(folder, f))
                            elif f == 'discrepency.report.txt':
                                with open(os.path.join(dirName, f)) as infile:
                                    discrep.write(infile.read())

# setup menu with argparse


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)


parser = argparse.ArgumentParser(prog='tbl2asn_parallel.py',
                                 description='''Run tbl2asn multipthreaded.''',
                                 epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
                                 formatter_class=MyFormatter)
parser.add_argument('-i', '--input', required=True,
                    help='Genome in TLB format')
parser.add_argument('-f', '--fasta', required=True,
                    help='Genome in FASTA format')
parser.add_argument('-s', '--species', required=True,
                    help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('-o', '--out', required=True, help='output directory')
parser.add_argument('--sbt', required=True, help='SBT file')
parser.add_argument('--isolate', help='Isolate name (e.g. Af293)')
parser.add_argument('--strain', help='Strain name (e.g. CEA10)')
parser.add_argument('-c', '--cpus', default=1, type=int, help='Number of CPUs')
parser.add_argument('-d', '--discrep', required=True,
                    help='discprepency report')
parser.add_argument('-t', '--tbl2asn', default='-l paired-ends',
                    help='Parameters for tbl2asn, linkage and gap info')
parser.add_argument('-v', '--version', default=1,
                    type=int, help='Genome version')
args = parser.parse_args()

if not os.path.isdir(args.out):
    os.makedirs(args.out)
# copy input files into directory
if args.input != os.path.join(args.out, 'genome.tbl'):
    shutil.copyfile(args.input, os.path.join(args.out, 'genome.tbl'))
if args.fasta != os.path.join(args.out, 'genome.fsa'):
    shutil.copyfile(args.fasta, os.path.join(args.out, 'genome.fsa'))
# now split input if genome is large enough
split_tbl2asn(args.out)
# now run in parallel
runtbl2asn_parallel(args.out, args.sbt, args.discrep, args.species,
                    args.isolate, args.strain, args.tbl2asn, args.version, args.cpus)
