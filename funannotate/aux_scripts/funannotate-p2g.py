#!/usr/bin/env python

import sys
import os
import uuid
import subprocess
import shutil
import itertools
import argparse
import timeit
import datetime
from Bio import SeqIO
import funannotate.library as lib
from pkg_resources import parse_version

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
parser.add_argument('--tmpdir', default='/tmp', help='volume to write tmp files')
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
    if parse_version(lib.getDiamondVersion()) >= parse_version('2.0.5'):
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
        lib.runSubprocess(db_cmd, output, lib.log, only_failed=True)
    else:
        lib.log.debug('Using premade Diamond database: {}'.format(premade_db))
        os.symlink(os.path.abspath(premade_db),
                   os.path.join(output, 'diamond.dmnd'))
    # now run search
    lib.runSubprocess(cmd, output, lib.log, only_failed=True)


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
           '-evalue', '1e-10',
           '-outfmt', '6 sseqid slen sstart send qseqid qlen qstart qend pident length evalue score',
           '-out', 'filter.tblastn.tab']
    lib.runSubprocess(cmd, output, lib.log)


def parseDiamond(blastresult):
    Results = {}
    with open(blastresult, 'r') as input:
        for line in input:
            cols = line.rstrip().split('\t')
            #hit = cols[0] + ':::' + cols[4]
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
                if cols[0] not in Results:
                    Results[cols[0]] = [(cols[4], start, end, int(cols[5]), float(cols[8]), float(cols[10]))]
                else:
                    Results[cols[0]].append((cols[4], start, end, int(cols[5]), float(cols[8]), float(cols[10])))
            else:
                lib.log.debug('P2G Error in coords: {:} start={:} stop={:} | {:}'.format(
                    coords, start, end, cols))
    return Results


def parseBlast(blastresult):
    # 6 'qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
    # '6 sseqid slen sstart send qseqid qlen qstart qend pident length evalue score'
    Results = {}
    with open(blastresult, 'r') as input:
        for line in input:
            cols = line.split('\t')
            hit = cols[4] + ':::' + cols[0]
            if int(cols[2]) < int(cols[3]):
                start = cols[2]
                end = cols[3]
            else:
                start = cols[3]
                end = cols[2]
            if not hit in Results:
                Results[hit] = (int(start), int(end), int(cols[1]), float(cols[8]), float(cols[10]))
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
                Results[hit] = (int(newstart), int(newstop), old[2], old[3], old[4])
    # convert Dictionary to a list that has  hit:::scaffold:::start:::stop
    filtered = {}
    for k, v in list(Results.items()):
        qseqid, sseqid = k.split(':::')
        if not qseqid in filtered:
            filtered[qseqid] = [(sseqid, v[0], v[1], v[2], v[3], v[4])]
        else:
            filtered[qseqid].append((sseqid, v[0], v[1], v[2], v[3], v[4]))
    return filtered


def spawn(cmd, **kwargs):
    # last item in cmd list is output file
    run_cmd = cmd[:-1]
    outfile = cmd[-1]
    with open(outfile, 'w') as output:
        p = subprocess.Popen(run_cmd, stdout=output, stderr=subprocess.PIPE, **kwargs)
    stderr = p.communicate()
    if p.returncode != 0:
        print('ERROR: {}'.format(' '.join(run_cmd)))
        if stderr:
            print(stderr.decode("utf-8"))
        return 1
    else:
        return 0

# write contig slices based on the diamond/tblastn parsed locations for staging exonerate runs
def contig_writer(scaff):
    # input here is dict items, so (contig, [(filename, start, end)])
    # just need to index the scaffold, then write all the files
    # old debugging
    # PID = os.getpid()
    #    lib.log.debug('PID {} working on {} splits of {}'.format(PID, len(scaff[1]), scaff[0]))
    record = SeqIO.index(os.path.join(tmpdir, 'scaffolds', '{}.fa'.format(scaff[0])),'fasta')
    recstr = record[scaff[0]] # move this outside the loop reduces the number of lookups required
                              # just reuse this object to create slices
    for z in scaff[1]:
        # I wonder if this would be faster if we just kept a set in memory and
        # skipped the isfile call and just checked to see if a filename was in our in-memory set?
        if not os.path.isfile(z[0]):
            with open(z[0], 'w') as outfile:
                outfile.write('>{}\n{}\n'.format(scaff[0], lib.softwrap(str(recstr[z[1]:z[2]].seq))))
#    lib.log.debug('PID {} for {} finished'.format(PID, scaff[0]))

# count number of proteins to look for
total = lib.countfasta(args.proteins)
lib.log.info('Mapping {:,} proteins to genome using {:} and exonerate'.format(
    total, args.filter))

# make tmpdir
tmpdir = os.path.join(args.tmpdir, 'p2g_' + str(uuid.uuid4()))
if not os.path.isdir(tmpdir):
    os.makedirs(tmpdir)
    os.makedirs(os.path.join(tmpdir, 'failed'))
    os.makedirs(os.path.join(tmpdir, 'scaffolds'))
    os.makedirs(os.path.join(tmpdir, 'proteins'))

filter_tic = timeit.default_timer()
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
    if args.tblastn:
        lib.log.info("Using pre-calculated Diamond result")
        BlastResult = args.tblastn
    else:
        BlastResult = os.path.join(tmpdir, 'diamond.matches.tab')
        runDiamond(os.path.abspath(args.genome), os.path.abspath(
            args.proteins), args.cpus, tmpdir, premade_db=args.filter_db)
    # parse results
    Hits = parseDiamond(BlastResult)
filter_toc = timeit.default_timer()
filter_runtime = filter_toc - filter_tic

# now generate the fasta files and build exonerate commands
# index the genome and proteins
fasta_tic = timeit.default_timer()
protein_dict = SeqIO.index(os.path.abspath(args.proteins), 'fasta')
genome_dict = SeqIO.index(os.path.abspath(args.genome), 'fasta')
exo_cmds = []
scaffold_splits = {}
# first write the proteins on single loop through
# collect the scaffold splits by each scaffold so can run this in parallel
# there seems to be a race condition if trying to access the same index at the same time from diff processes

for k, v in Hits.items():
    protfile = os.path.join(tmpdir, 'proteins', k.replace('|', '_')+'.fasta')
    if not os.path.isfile(protfile):
        # I tested and there was no signif difference with using Biopython SeqIO.write
        # but leaving this string -> softwrap
        protSeq = str(protein_dict[k].seq)
        with open(protfile, 'w') as outfile:
            outfile.write('>{}\n{}\n'.format(k, lib.softwrap(protSeq)))
    for i, x in enumerate(sorted(v, key=lambda y: y[4], reverse=True)):
        if i <= args.ploidy*2:  # look for 2 hits for each copy
            # lets add 3kb in each direction to make sure exonerate can find it
            start = x[1] - 3000
            if start < 0:
                start = 0
            end = x[2] + 3000
            if end > x[3]:
                end = x[3]
            exoname = k+'.'+x[0]+'__'+str(start)+'__'
            # check that input files are created and valid
            exonerate_out = os.path.join(tmpdir, 'exonerate.' + exoname + '.out')
            dnafile = os.path.join(tmpdir, 'scaffolds', '{}_{}-{}.fasta'.format(x[0], start, end))
            if not x[0] in scaffold_splits:
                scaffold_splits[x[0]] = [(dnafile, start, end)]
                # write the whole scaffold to file to use later
                with open(os.path.join(tmpdir, 'scaffolds', '{}.fa'.format(x[0])), 'w') as outfile:
                    outfile.write('>{}\n{}\n'.format(x[0], lib.softwrap(str(genome_dict[x[0]].seq))))
            else:
                scaffold_splits[x[0]].append((dnafile, start, end))
            # now build exonerate command and add to list
            cmd = ['exonerate', '--model', 'p2g', '--showvulgar', 'no',
                '--showalignment', 'no', '--showquerygff', 'no',
                '--showtargetgff', 'yes', '--maxintron', str(args.maxintron),
                '--percent', str(args.exonerate_pident),
                '--ryo', "AveragePercentIdentity: %pi\n", protfile, dnafile, exonerate_out]
            exo_cmds.append(cmd)

# clear memory for these objects as needed
genome_dict = None
protein_dict = None
# now we want to run the contig writer across multiple processes
lib.log.debug('Writing {} contig splits for exonerate'.format(len(exo_cmds)))
lib.runMultiProgress(contig_writer, scaffold_splits.items(), args.cpus, progress=False)

fasta_toc = timeit.default_timer()
fasta_runtime = fasta_toc - fasta_tic
lib.log.info('Found {:,} preliminary alignments with {:} in {:} --> generated FASTA files for exonerate in {:}'.format(
    len(exo_cmds), args.filter, str(datetime.timedelta(seconds=filter_runtime)).rsplit('.', 1)[0],
    str(datetime.timedelta(seconds=fasta_runtime)).rsplit('.', 1)[0]))

# run multiprocessing exonerate
exo_tic = timeit.default_timer()
lib.runMultiProgress(spawn, exo_cmds, args.cpus,
                     progress=args.progress)
exo_toc = timeit.default_timer()
exo_runtime = exo_toc - exo_tic
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
lib.log.info('Exonerate finished in {:}: found {:,} alignments'.format(
    str(datetime.timedelta(seconds=exo_runtime)).rsplit('.', 1)[0], Found))

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
