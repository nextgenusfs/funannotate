#!/usr/bin/env python

import sys
import multiprocessing
import subprocess
import os, re
import time
import shutil
import argparse
import math
from Bio import SeqIO
from natsort import natsorted
import funannotate.library as lib


def worker(inputList):
    output = inputList[-2]
    logfile = inputList[-1]
    cmd = inputList[:-2]
    with open(output, 'w') as outfile:
        with open(logfile, 'w') as log:
            subprocess.call(cmd, stdout=outfile, stderr=log)


def safe_run(*args, **kwargs):
    """Call run(), catch exceptions."""
    try:
        worker(*args, **kwargs)
    except Exception as e:
        print(("error: %s run(*%r, **%r)" % (e, args, kwargs)))


def create_partitions(fasta, genes, partition_list, proteins=False,
                      transcripts=False, repeats=False, num=50,
                      tmpdir='.', interval=2000, partitions=True):
    # function to create EVM partition intervals that do not split genes
    if not os.path.isdir(tmpdir):
        os.makedirs(tmpdir)
    SeqRecords = SeqIO.index(fasta, 'fasta')
    PID = os.getpid()
    bedGenes = os.path.join(tmpdir, 'genes.{}.bed'.format(PID))
    Results = []
    with open(genes, 'r') as infile:
        for line in infile:
            if line.startswith('#') or line.startswith('\n'):
                continue
            line = line.rstrip()
            cols = line.split('\t')
            if cols[2] == 'gene':
                Results.append([cols[0], int(cols[3]), int(cols[4]), cols[8],
                                cols[5], cols[6]])
    # sort the results by contig and position
    ChrGeneCounts = {}
    sortedResults = natsorted(Results, key=lambda x: (x[0], x[1]))
    with open(bedGenes, 'w') as outfile:
        for x in sortedResults:
            outfile.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(x[0], x[1], x[2],
                                                            x[3], x[4], x[5]))
            if not x[0] in ChrGeneCounts:
                ChrGeneCounts[x[0]] = 1
            else:
                ChrGeneCounts[x[0]] += 1
    ChrNoGenes = len(SeqRecords) - len(ChrGeneCounts)
    lib.log.debug('{:,} total contigs; skipping {:,} contigs with no genes'.format(len(SeqRecords), ChrNoGenes))
    if partitions:
        # now merge overlaping genes [strand] to get conservative locus boundaries
        cmd = ['bedtools', 'merge', '-s', '-i', bedGenes]
        merged = {}
        for line in lib.execute(cmd):
            line = line.rstrip()
            chr, start, end = line.split('\t')
            if chr not in merged:
                merged[chr] = [(int(start), int(end), -1)]
            else:
                diff = int(start) - merged[chr][-1][1]
                merged[chr].append((int(start), int(end), diff))
        # parse Results and get coordinates to partitions
        Partitions = {}
        Commands = {}
        for k, v in natsorted(merged.items()):
            if not k in ChrGeneCounts:  # no genes, so can safely skip
                continue
            Partitions[k] = []
            if len(v) > num:
                chunks = math.ceil(len(v) / num)
                num_genes = int(round(len(v) / chunks))
                chunks = int(chunks)
                for i in range(chunks):
                    if k in Commands:
                        continue
                    i = i + 1
                    if i == 1:
                        start = 1
                    else:
                        phase = int(round(interval / 4))
                        if len(Partitions[k]) > 0:
                            start = Partitions[k][-1][1] - phase
                        else:
                            start = 1
                    loc = i * num_genes
                    if i == chunks:
                        end = len(SeqRecords[k])
                    else:
                        if loc >= len(v):
                            end = len(SeqRecords[k])
                        else:
                            end = getBreakPoint(v, loc, direction='reverse',
                                                gap=interval)
                            if not end:
                                end = getBreakPoint(v, loc, direction='forward',
                                                    gap=interval)
                    if not end:
                        Commands[k] = {'n': len(v)}
                        lib.log.debug('{} --> {} bp'.format(
                            k, len(SeqRecords[k])))
                    else:
                        partLen = end - start
                        if partLen < 10000:
                            continue
                        Partitions[k].append((start, end))
                        partName = '{}_{}-{}'.format(k, start, end)
                        Commands[partName] = {'n': num_genes, 'chr': k}
                        lib.log.debug('{} --> {} bp'.format(
                            partName, partLen))
            else:
                Commands[k] = {'n': len(v)}
                lib.log.debug('{} --> {} bp'.format(
                    k, len(SeqRecords[k])))
        # now loop through partitions and write files for EVM
        with open(partition_list, 'w') as partout:
            for chr, p in natsorted(Partitions.items()):
                chrDir = os.path.join(tmpdir, chr)
                if not os.path.isdir(chrDir):
                    os.makedirs(chrDir)
                if len(p) == 0:
                    partout.write('{}\t{}\t{}\n'.format(chr,
                                                        os.path.abspath(
                                                            chrDir),
                                                        'N'))
                    chrFasta = os.path.join(chrDir, os.path.basename(fasta))
                    with open(chrFasta, 'w') as fastaout:
                        fastaout.write('>{}\n{}\n'.format(
                            chr, lib.softwrap(str(SeqRecords[chr].seq))))
                    genePred = os.path.join(chrDir, os.path.basename(genes))
                    RangeFinder(genes, chr, 1, len(SeqRecords[chr]),
                                genePred)
                    if proteins:
                        protPred = os.path.join(
                            chrDir, os.path.basename(proteins))
                        RangeFinder(proteins, chr, 1, len(SeqRecords[chr]),
                                    protPred)
                    if transcripts:
                        tranPred = os.path.join(chrDir,
                                                os.path.basename(transcripts))
                        RangeFinder(transcripts, chr, 1, len(SeqRecords[chr]),
                                    tranPred)
                    if repeats:
                        repPred = os.path.join(
                            chrDir, os.path.basename(repeats))
                        RangeFinder(repeats, chr, 1, len(SeqRecords[chr]),
                                    repPred)
                else:
                    for coords in p:
                        partDir = os.path.join(
                            chrDir, '{}_{}-{}'.format(chr, coords[0], coords[1]))
                        if not os.path.isdir(partDir):
                            os.makedirs(partDir)
                        partout.write('{}\t{}\t{}\t{}\n'.format(
                            chr, os.path.abspath(chrDir), 'Y',
                            os.path.abspath(partDir)))
                        partFasta = os.path.join(
                            partDir, os.path.basename(fasta))
                        with open(partFasta, 'w') as fastaout:
                            fastaout.write('>{}\n{}\n'.format(
                                chr, lib.softwrap(
                                    str(SeqRecords[chr].seq[coords[0] - 1:coords[1]])
                                    )))
                        # split genes GFF3
                        genePred = os.path.join(
                            partDir, 'gene_predictions.gff3')
                        RangeFinder(genes, chr, coords[0], coords[1],
                                    genePred)
                        if proteins:
                            protPred = os.path.join(partDir,
                                                    os.path.basename(proteins))
                            RangeFinder(proteins, chr, coords[0], coords[1],
                                        protPred)
                        if transcripts:
                            tranPred = os.path.join(partDir,
                                                    os.path.basename(transcripts))
                            RangeFinder(transcripts, chr, coords[0], coords[1],
                                        tranPred)
                        if repeats:
                            repPred = os.path.join(partDir,
                                                   os.path.basename(repeats))
                            RangeFinder(repeats, chr, coords[0], coords[1],
                                        repPred)
    else:
        Commands = {}
        with open(partition_list, 'w') as partout:
            for chr in SeqRecords:
                if not chr in ChrGeneCounts:  # no genes so skip
                    continue
                Commands[chr] = {'n': len(SeqRecords[chr])}
                chrDir = os.path.join(tmpdir, chr)
                if not os.path.isdir(chrDir):
                    os.makedirs(chrDir)
                partout.write('{}\t{}\t{}\n'.format(chr,
                                                    os.path.abspath(chrDir),
                                                    'N'))
                chrFasta = os.path.join(chrDir, os.path.basename(fasta))
                with open(chrFasta, 'w') as fastaout:
                    fastaout.write('>{}\n{}\n'.format(
                        chr, lib.softwrap(str(SeqRecords[chr].seq))))
                genePred = os.path.join(chrDir, os.path.basename(genes))
                RangeFinder(genes, chr, 1, len(SeqRecords[chr]),
                            genePred)
                if proteins:
                    protPred = os.path.join(chrDir, os.path.basename(proteins))
                    RangeFinder(proteins, chr, 1, len(SeqRecords[chr]),
                                protPred)
                if transcripts:
                    tranPred = os.path.join(chrDir,
                                            os.path.basename(transcripts))
                    RangeFinder(transcripts, chr, 1, len(SeqRecords[chr]),
                                tranPred)
                if repeats:
                    repPred = os.path.join(chrDir, os.path.basename(repeats))
                    RangeFinder(repeats, chr, 1, len(SeqRecords[chr]),
                                repPred)

    return Commands


def RangeFinder(input, chrom, start, end, output, EVM=False):
#    if not EVM:
#        EVM = os.environ['EVM_HOME']
#    RFScript = os.path.join(EVM, 'EvmUtils', 'gff_range_retriever.pl')
#    cmd = [RFScript, chrom, str(start), str(end), 'ADJUST_TO_ONE']
#    with open(output, 'w') as outfile:
#        with open(os.path.abspath(input)) as infile:
#            p = subprocess.Popen(cmd, stdin=infile, stdout=outfile)
#            p.wait()

    adjust_coord = start - 1;
#    print("end is %d, adjust_coord is %d chrom is %s"%(chrom,end,adjust_coord))
    got_spacer = False
    adjust_to_one = True  # not sure when this would bbe false
    with open(output, 'w') as outfile:
        with open(os.path.abspath(input)) as infile:
            for line in infile:
                if not re.search(r'\w',line):
                    if not got_spacer:
                        outfile.write(line)
                    got_spacer = True
                elif line.startswith('#'):
                    outfile.write(line)
                else:
                    row = line.split("\t")
                    row[3] = int(row[3])
                    row[4] = int(row[4])
                    (contig, lend, rend) = (row[0], row[3], row[4])
                    if contig == chrom and lend >= start and rend <= end:
                        if adjust_to_one:
                            row[3] -= adjust_coord
                            row[4] -= adjust_coord
                        row[3] = "%d"%(row[3])
                        row[4] = "%d"%(row[4])
                        outfile.write("\t".join(row))
                        got_spacer = False
        outfile.flush()


def getBreakPoint(tupList, idx, direction='reverse', gap=2000):
    # takes list of tuples of coords and a starting index (idx). finds closest
    # break point in between tuple coordSorted
    solution = False
    while not solution:
        try:
            start, end, diff = tupList[idx]
        except IndexError:
            return False
        if diff >= gap:
            phase = int(round(diff / 2))
            solution = end + phase
        else:
            if direction == 'reverse':
                idx -= 1
            else:
                idx += 1
    return solution


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)


parser = argparse.ArgumentParser(
    prog='funannotate-runEVM.py',
    description='''Funannotate script to run evidencemodeler.''',
    epilog="""Written by Jon Palmer (2020) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)
parser.add_argument('-f', '--fasta', required=True, help='genome fasta')
parser.add_argument('-g', '--genes', required=True, help='gene predictions')
parser.add_argument('-w', '--weights', required=True, help='weights file')
parser.add_argument('-o', '--out', required=True, help='output GFF3')
parser.add_argument('-p', '--proteins', help='protein alignments')
parser.add_argument('-t', '--transcripts', help='transcript alignments')
parser.add_argument('-r', '--repeats', help='repeats GFF3')
parser.add_argument('-i', '--interval', type=int, default=2000,
                    help='length of gene free interval to use')
parser.add_argument('-n', '--num-gene-partition', dest='gene_partition',
                    type=int, default=35,
                    help='number of features per partition')
parser.add_argument('-m', '--min-intron-len', dest='min_intron', type=int,
                    default=10, help='minimum intron')
parser.add_argument('--no-partitions', dest='no_partitions',
                    action='store_false', help='no contig splits')
parser.add_argument('-c', '--cpus', type=int, default=2, help='num cpus')
parser.add_argument('-d', '--dir', default='.', help='directory')
parser.add_argument('-l', '--logfile', help='logfile')
parser.add_argument('--EVM_HOME',
                    help='Path to Evidence Modeler home directory, $EVM_HOME')
args = parser.parse_args()

# initialize script, log system info and cmd issue at runtime
PID = os.getpid()
if args.logfile:
    log_name = args.logfile
else:
    log_name = 'funannotate-evm.{}.log'.format(PID)
if os.path.isfile(log_name):
    os.remove(log_name)
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv) + '\n'
lib.log.debug(cmd_args)

# create output directory
tmpdir = args.dir
if os.path.exists(tmpdir):
    shutil.rmtree(tmpdir)
os.makedirs(tmpdir)

# set some EVM script locations
perl = 'perl'
if args.EVM_HOME:
    EVM = args.EVM_HOME
else:
    try:
        EVM = os.environ['EVM_HOME']
    except:
        lib.log.error('Could not find EVM_HOME environmental variable')
        sys.exit(1)

Combine = os.path.join(EVM, 'EvmUtils', 'recombine_EVM_partial_outputs.pl')
Convert = os.path.join(EVM, 'EvmUtils', 'convert_EVM_outputs_to_GFF3.pl')

# split partitions
partitions = os.path.join(tmpdir, 'partitions_list.out')
if args.no_partitions:
    lib.log.info('EVM: partitioning input to ~ {} genes per partition'.format(
        args.gene_partition))
else:
    lib.log.info('EVM: partitioning each contig separately')
cmdinfo = create_partitions(args.fasta, args.genes, partitions,
                            proteins=args.proteins,
                            transcripts=args.transcripts,
                            repeats=args.repeats,
                            tmpdir=tmpdir, num=args.gene_partition,
                            interval=args.interval,
                            partitions=args.no_partitions)
lib.log.debug('Finished partitioning, generating command list')
# sort the cmdinfo by number of putative genes and distribute into sub files
num_workers = args.cpus - 1
if num_workers < 1:
    num_workers = 1
c = 0
file_list = []
tasks = 0
for s in sorted(cmdinfo.items(), key=lambda x: x[1]['n'], reverse=True):
    key, d = s
    if 'chr' in d:
        outputDir = os.path.abspath(os.path.join(tmpdir, d['chr'], key))
    else:
        outputDir = os.path.abspath(os.path.join(tmpdir, key))
    cmd = [os.path.join(EVM, 'evidence_modeler.pl'),
           '-G', os.path.basename(args.fasta),
           '-g', os.path.basename(args.genes),
           '-w', os.path.abspath(args.weights),
           '--min_intron_length', str(args.min_intron),
           '--exec_dir', outputDir]
    if args.proteins:
        cmd += ['-p', os.path.basename(args.proteins)]
    if args.transcripts:
        cmd += ['-e', os.path.basename(args.transcripts)]
    if args.repeats:
        cmd += ['--repeats', os.path.basename(args.repeats)]
    cmd += [os.path.join(outputDir, 'evm.out'),
            os.path.join(outputDir, 'evm.out.log')]
    file_list.append(cmd)

# run runMultiProgress
lib.runMultiProgress(safe_run, file_list, num_workers)

# now combine the paritions
cmd4 = [perl, Combine, '--partitions', os.path.basename(partitions),
        '--output_file_name', 'evm.out']
lib.runSubprocess(cmd4, tmpdir, lib.log)

# now convert to GFF3
cmd5 = [perl, Convert, '--partitions', os.path.basename(partitions),
        '--output', 'evm.out',
        '--genome', os.path.abspath(args.fasta)]
lib.runSubprocess(cmd5, tmpdir, lib.log)

# now concatenate all GFF3 files together for a genome then
lib.log.info("Converting to GFF3 and collecting all EVM results")
with open(args.out, 'w') as out:
    for root, dirs, files in os.walk(tmpdir):
        for file in files:
            if file == 'evm.out.gff3':
                filename = os.path.join(root, file)
                with open(filename, 'r') as readfile:
                    shutil.copyfileobj(readfile, out)
