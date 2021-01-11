#!/usr/bin/env python

import sys
import os
import uuid
import subprocess
import shutil
import itertools
import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
import funannotate.library as lib

# setup menu with argparse


class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)


parser = argparse.ArgumentParser(
    prog='funannotate-p2g.py',
    description='''Funannotate script to run tblastn/exonerate protein2genome.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class=MyFormatter)
parser.add_argument('-p', '--proteins', required=True,
                    help='Protein multi-fasta input')
parser.add_argument('-g', '--genome', required=True,
                    help='Genome multi-fasta input')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs')
parser.add_argument('--tblastn', help='tBLASTN run previously')
parser.add_argument('-o', '--out', required=True,
                    help='Final exonerate output file')
parser.add_argument('-t', '--tblastn_out', help='Save tblastn output')
parser.add_argument('--maxintron', default=3000, help='Maximum intron size')
parser.add_argument('--exonerate_pident', default=80,
                    help='Exonerate pct identity')
parser.add_argument('--logfile', default='funannotate-p2g.log', help='logfile')
parser.add_argument('--ploidy', default=1, type=int, help='Ploidy of assembly')
parser.add_argument('--debug', action='store_true',
                    help='Keep intermediate folders if error detected')
parser.add_argument('-f', '--filter', default='diamond', choices=[
                    'diamond', 'tblastn'],
                    help='Method to use for pre-filter for exonerate')
parser.add_argument('-d', '--filter_db',
                    help='Premade diamond protein database')
parser.add_argument('--EVM_HOME',
                    help='Path to Evidence Modeler home directory, $EVM_HOME')
parser.add_argument('--no-progress', dest='progress', action='store_false',
                    help='no progress on multiprocessing')
args = parser.parse_args()

# do some checks and balances
if args.EVM_HOME:
    EVM = args.EVM_HOME
else:
    try:
        EVM = os.environ["EVM_HOME"]
    except KeyError:
        lib.log.error("$EVM_HOME environmental variable not found, Evidence Modeler is not properly configured.  You can use the --EVM_HOME argument to specifiy a path at runtime")
        sys.exit(1)

ExoConverter = os.path.join(EVM, 'EvmUtils', 'misc',
                            'exonerate_gff_to_alignment_gff3.pl')

log_name = args.logfile
if os.path.isfile(log_name):
    os.remove(log_name)

# initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

# get version of programs
exo_version = subprocess.Popen(
    ['exonerate', '--version'], stdout=subprocess.PIPE,
    universal_newlines=True).communicate()[0].split('\n')[0]
exo_version = exo_version.split('version ')[-1]
blast_version = subprocess.Popen(
    ['tblastn', '-version'], stdout=subprocess.PIPE,
    universal_newlines=True).communicate()[0].split('\n')[0]
blast_version = blast_version.split(': ')[-1]
if args.filter == 'diamond':
    diamond_version = lib.getDiamondVersion()


def runDiamond(input, query, cpus, output, premade_db=None):
    # create DB of protein sequences
    # check diamond version
    if lib.getDiamondVersion() >= '2.0.5':
        # run in frameshift mode
        cmd = ['diamond', 'blastx', '--threads', str(cpus), '-q', input,
               '--db', 'diamond', '-o', 'diamond.matches.tab', '-e', '1e-10',
               '-k', '0', '--more-sensitive', '--unal', '0', '-c', '1', '-F', '15',
               '-f', '6', 'sseqid', 'slen',
               'sstart', 'send', 'qseqid', 'qlen', 'qstart', 'qend', 'pident',
               'length', 'evalue', 'score', 'qcovhsp', 'qframe']
    else:
        if int(cpus) > 8:
            cpus = 8
        cmd = ['diamond', 'blastx', '--threads', str(cpus), '-q', input,
               '--db', 'diamond', '-o', 'diamond.matches.tab', '-e', '1e-10',
               '-k', '0', '--more-sensitive', '-f', '6', 'sseqid', 'slen',
               'sstart', 'send', 'qseqid', 'qlen', 'qstart', 'qend', 'pident',
               'length', 'evalue', 'score', 'qcovhsp', 'qframe']
    if premade_db is None:
        db_cmd = ['diamond', 'makedb', '--threads',
                  str(cpus), '--in', query, '--db', 'diamond']
        lib.runSubprocess4(db_cmd, output, lib.log)
    else:
        lib.log.debug('Using premade Diamond database: {}'.format(premade_db))
        os.symlink(os.path.abspath(premade_db),
                   os.path.join(output, 'diamond.dmnd'))
    # now run search
    lib.runSubprocess4(cmd, output, lib.log)


def runtblastn(input, query, cpus, output, maxhits):
    # start by formatting blast db/dustmasker filtered format
    cmd = ['dustmasker', '-in', input, '-infmt', 'fasta', '-parse_seqids',
           '-outfmt', 'maskinfo_asn1_bin', '-out', 'genome_dust.asnb']
    lib.runSubprocess(cmd, output, lib.log)
    cmd = ['makeblastdb', '-in', input, '-dbtype', 'nucl',
           '-parse_seqids', '-mask_data', 'genome_dust.asnb', '-out', 'genome']
    lib.runSubprocess(cmd, output, lib.log)
    cmd = ['tblastn', '-num_threads', str(cpus), '-db', 'genome',
           '-query', query, '-max_target_seqs', str(maxhits),
           '-db_soft_mask', '11',
           '-threshold', '999', '-max_intron_length', str(args.maxintron),
           '-evalue', '1e-10', '-outfmt', '6', '-out', 'filter.tblastn.tab']
    lib.runSubprocess(cmd, output, lib.log)


def parseDiamond(blastresult):
    Results = {}
    with open(blastresult, 'r') as input:
        for line in input:
            cols = line.rstrip().split('\t')
            hit = cols[0] + ':::' + cols[4]
            coords = [int(cols[6]), int(cols[7])]
            start_extend = (int(cols[2])*3) - 3
            end_extend = (int(cols[1]) - int(cols[3]))*3
            qframe = int(cols[13])
            if qframe > 0:
                start = int(cols[6]) - start_extend
                end = int(cols[7]) + end_extend
            else:
                start = int(cols[7]) - end_extend
                end = int(cols[6]) + start_extend
            # make sure not below zero
            if start < 1:
                start = 1
            # make sure not greater than length of contig
            if end > int(cols[5]):
                end = int(cols[5])
            if end > start:
                Results[hit] = (start, end)
            else:
                lib.log.debug('P2G Error in coords: {:} start={:} stop={:} | {:}'.format(
                    coords, start, end, cols))
    # convert Dictionary to a list that has  hit:::scaffold:::start:::stop
    HitList = []
    for k, v in list(Results.items()):
        finalhit = k+':::'+str(v[0])+':::'+str(v[1])
        HitList.append(finalhit)
    return HitList


def parseBlast(blastresult):
    Results = {}
    with open(blastresult, 'r') as input:
        for line in input:
            cols = line.split('\t')
            hit = cols[0] + ':::' + cols[1]
            if int(cols[8]) < int(cols[9]):
                start = cols[8]
                end = cols[9]
            else:
                start = cols[9]
                end = cols[8]
            if not hit in Results:
                Results[hit] = (start, end)
            else:
                # get old start stop
                old = Results.get(hit)
                if int(start) < int(old[0]):
                    newstart = start
                else:
                    newstart = old[0]
                if int(end) > int(old[1]):
                    newstop = end
                else:
                    newstop = old[1]
                Results[hit] = (newstart, newstop)
    # convert Dictionary to a list that has  hit:::scaffold:::start:::stop
    HitList = []
    for k, v in list(Results.items()):
        finalhit = k+':::'+str(v[0])+':::'+str(v[1])
        HitList.append(finalhit)
    return HitList


def runExonerate(input):
    s = input.split(':::')
    ProtID = s[0]
    ScaffID = s[1]
    ScaffStart = int(s[2])
    ScaffEnd = int(s[3])
    # get the protein model
    query = os.path.join(tmpdir, ProtID+'.'+str(uuid.uuid4())+'.fa')
    with open(query, 'w') as output:
        SeqIO.write(protein_dict[ProtID], output, 'fasta')
    # now get the genome region, use different variable names for SeqRecords to avoid collision
    scaffold = os.path.join(tmpdir, ScaffID+'.'+ProtID +
                            '.'+str(ScaffStart)+'-'+str(ScaffEnd)+'.fa')
    with open(scaffold, 'w') as output2:
        with open(os.path.join(tmpdir, 'scaffolds', ScaffID+'.fa'), 'r') as fullscaff:
            for header, Sequence in SimpleFastaParser(fullscaff):
                # grab a 3 kb cushion on either side of hit region, careful of scaffold ends
                start = ScaffStart - 3000
                if start < 1:
                    start = 1
                end = ScaffEnd + 3000
                if end > len(Sequence):
                    end = len(Sequence)
                output2.write('>%s\n%s\n' % (header, Sequence[start:end]))
    exoname = ProtID+'.'+ScaffID+'__'+str(start)+'__'
    # check that input files are created and valid
    exonerate_out = os.path.join(tmpdir, 'exonerate.' + exoname + '.out')
    ryo = "AveragePercentIdentity: %pi\n"
    cmd = ['exonerate', '--model', 'p2g', '--showvulgar', 'no',
           '--showalignment', 'no', '--showquerygff', 'no',
           '--showtargetgff', 'yes', '--maxintron', str(args.maxintron),
           '--percent', str(args.exonerate_pident),
           '--ryo', ryo, query, scaffold]
    if lib.checkannotations(query) and lib.checkannotations(scaffold):
        # run exonerate, capture errors
        with open(exonerate_out, 'w') as output3:
            proc = subprocess.Popen(
                cmd, stdout=output3, stderr=subprocess.PIPE,
                universal_newlines=True)
        stderr = proc.communicate()
        if 'WARNING' in stderr[1]:
            lib.log.debug('Error in input:{:}'.format(input))
            lib.log.debug('%s, Len=%i, %i-%i; %i-%i' %
                          (header, len(Sequence), ScaffStart,
                           ScaffEnd, start, end))
            os.rename(query, os.path.join(
                tmpdir, 'failed', os.path.basename(query)))
            os.rename(scaffold, os.path.join(
                tmpdir, 'failed', os.path.basename(scaffold)))
        else:
            for y in [query, scaffold]:
                try:
                    lib.SafeRemove(y)
                except OSError:
                    lib.log.debug("Error removing %s" % (y))
        # check filesize of exonerate output, no hits still have some output
        # data in them, should be safe dropping anything smaller than 500 bytes
        if lib.getSize(exonerate_out) < 500:
            os.remove(exonerate_out)
    else:
        lib.log.debug('Error in query or scaffold:{:}'.format(input))
        lib.SafeRemove(query)
        lib.SafeRemove(scaffold)


# count number of proteins to look for
total = lib.countfasta(args.proteins)
lib.log.info('Mapping {:,} proteins to genome using {:} and exonerate'.format(
    total, args.filter))

# make tmpdir
tmpdir = 'p2g_' + str(uuid.uuid4())
if not os.path.isdir(tmpdir):
    os.makedirs(tmpdir)
    os.makedirs(os.path.join(tmpdir, 'failed'))
    os.makedirs(os.path.join(tmpdir, 'scaffolds'))

if args.filter == 'tblastn':
    lib.log.debug("BLAST v%s; Exonerate v%s" % (blast_version, exo_version))
    # check for tblastn input
    if args.tblastn:
        lib.log.info("Using pre-calculated tBLASTN result")
        BlastResult = args.tblastn
    else:
        BlastResult = os.path.join(tmpdir, 'filter.tblastn.tab')
        runtblastn(os.path.abspath(args.genome), os.path.abspath(
            args.proteins), args.cpus, tmpdir, args.ploidy*5)  # 2X ploidy for tBLASTn filter
    # parse results
    Hits = parseBlast(BlastResult)
else:
    lib.log.debug("Diamond v%s; Exonerate v%s" %
                  (diamond_version, exo_version))
    # run Diamond
    BlastResult = os.path.join(tmpdir, 'diamond.matches.tab')
    runDiamond(os.path.abspath(args.genome), os.path.abspath(
        args.proteins), args.cpus, tmpdir, premade_db=args.filter_db)
    Hits = parseDiamond(BlastResult)

lib.log.info('Found {0:,}'.format(len(Hits)) +
             ' preliminary alignments --> aligning with exonerate')

# index the genome and proteins
# do index here in case memory problems?
protein_dict = SeqIO.index(os.path.abspath(args.proteins), 'fasta')

# split genome fasta into individual scaffolds
with open(os.path.abspath(args.genome), 'r') as input:
    for record in SeqIO.parse(input, "fasta"):
        SeqIO.write(record, os.path.join(
            tmpdir, 'scaffolds', record.id + ".fa"), "fasta")

# run multiprocessing exonerate
lib.runMultiProgress(runExonerate, Hits, args.cpus,
                     progress=args.progress)

# now need to loop through and offset exonerate predictions back to whole scaffolds
exonerate_raw = os.path.join(tmpdir, 'exonerate.out.combined')
with open(exonerate_raw, 'w') as output:
    for file in os.listdir(tmpdir):
        if file.endswith('.out'):
            with open(os.path.join(tmpdir, file), 'r') as exoresult:
                offset = int(file.split('__')[1])
                for line in itertools.islice(exoresult, 3, None):
                    if line.startswith('#') or line.startswith('Average') or line.startswith('-- completed'):
                        output.write(line)
                    else:
                        cols = line.split('\t')
                        cols[3] = str(int(cols[3])+offset)
                        cols[4] = str(int(cols[4])+offset)
                        output.write('\t'.join(cols))

# convert to GFF3 using ExoConverter from EVM
with open(args.out, 'w') as output:
    subprocess.call([ExoConverter, exonerate_raw], stdout=output, stderr=FNULL)

# output some quick summary of exonerate alignments that you found
Found = lib.countGFFgenes(exonerate_raw)
lib.log.info('Exonerate finished: found {:,} alignments'.format(Found))

# check for saving output of tblastn
if args.tblastn_out:
    shutil.copyfile(BlastResult, args.tblastn_out)

# finally clean-up your mess if failed is empty
if args.debug:
    try:
        os.rmdir(os.path.join(tmpdir, 'failed'))
        empty = True
    except OSError:
        empty = False
    if empty:
        lib.SafeRemove(tmpdir)
    else:
        lib.log.error("Failed exonerate alignments found, see files in %s" %
                      os.path.join(tmpdir, 'failed'))
else:
    if os.path.isdir(tmpdir):
        lib.SafeRemove(tmpdir)
