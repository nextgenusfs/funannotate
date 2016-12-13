#!/usr/bin/env python

import sys, os, subprocess, csv, shutil, inspect, itertools, argparse
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-p2g.py',
    description='''Funannotate script to run tblastn/exonerate protein2genome.''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-p','--proteins', required=True, help='Protein multi-fasta input')
parser.add_argument('-g','--genome', required=True, help='Genome multi-fasta input')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs')
parser.add_argument('--tblastn', help='tBLASTN run previously')
parser.add_argument('-o','--out', required=True, help='Final exonerate output file')
parser.add_argument('--maxintron', default = 3000, help='Maximum intron size')
parser.add_argument('--logfile', default ='funannotate-p2g.log', help='logfile')
args=parser.parse_args() 


log_name = args.logfile
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)

#get version of programs
exo_version = subprocess.Popen(['exonerate', '--version'], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
exo_version = exo_version.split('version ')[-1]
blast_version = subprocess.Popen(['tblastn', '-version'], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
blast_version = blast_version.split(': ')[-1]
lib.log.debug("BLAST v%s; Exonerate v%s" % (blast_version, exo_version))

def tblastnFilter(input, query, cpus, output):
    #start by formatting blast db/dustmasker filtered format
    cmd = ['dustmasker', '-in', input, '-infmt', 'fasta', '-parse_seqids', '-outfmt', 'maskinfo_asn1_bin', '-out', 'genome_dust.asnb']
    lib.runSubprocess(cmd, output, lib.log)
    cmd = ['makeblastdb', '-in', input, '-dbtype', 'nucl', '-parse_seqids', '-mask_data', 'genome_dust.asnb', '-out', 'genome']
    lib.runSubprocess(cmd, output, lib.log)
    cmd = ['tblastn', '-num_threads', str(cpus), '-db', 'genome', '-query', query, '-max_target_seqs', '1', '-db_soft_mask', '11', '-threshold', '999', '-max_intron_length', str(args.maxintron), '-evalue', '1e-10', '-outfmt', '6', '-out', 'filter.tblastn.tab']
    lib.runSubprocess(cmd, output, lib.log)

def parseBlast(blastresult):
    global HitList
    HitList = []
    #now parse through results, generating a list for exonerate function
    with open(blastresult, 'rU') as input:
        reader = csv.reader(input, delimiter='\t')
        for cols in reader:
            hit = cols[0] + '::' + cols[1]
            if hit not in HitList:
                HitList.append(hit)

def runExonerate(input):
    FNULL = open(os.devnull, 'w')
    s = input.split('::')
    if s[0].startswith('sp|'):
        name = s[0].split("|")[1] + '_' + s[1]
    else:
        name = s[0].split()[0] + '_' + s[1]
    query = os.path.join(tmpdir, name+'.fa')
    with open(query, 'w') as output:
        rec = record_dict[s[0]]
        output.write(">%s\n%s\n" % (rec.id, rec.seq))
    scaffold = s[1] + '.fa'
    scaffold = os.path.join(tmpdir, scaffold)
    exonerate_out = 'exonerate_' + name + '.out'
    exonerate_out = os.path.join(tmpdir, exonerate_out)
    ryo = "AveragePercentIdentity: %pi\n"
    with open(exonerate_out, 'w') as output:
        subprocess.call(['exonerate', '--model', 'p2g', '--showvulgar', 'no', '--showalignment', 'no', '--showquerygff', 'no', '--showtargetgff', 'yes', '--maxintron', str(args.maxintron), '--percent', '80', '--ryo', ryo , query, scaffold], stdout = output, stderr = FNULL)
    os.remove(query)
    #check filesize of exonerate output, no hits are 285 bytes, but lets just filter everything smaller than 310
    if lib.getSize(exonerate_out) < 310:
        os.remove(exonerate_out)

#make tmpdir
tmpdir = 'p2g_'+ str(os.getpid())
if not os.path.isdir(tmpdir):
    os.makedirs(tmpdir)

if args.tblastn:
    lib.log.info("Using pre-calculated tBLASTN result")
    BlastResult = args.tblastn
else:
    lib.log.info("Running pre-filter tBlastn step")
    BlastResult = os.path.join(tmpdir, 'filter.tblastn.tab')
    tblastnFilter(os.path.abspath(args.genome), os.path.abspath(args.proteins), args.cpus, tmpdir)

#parse the results
parseBlast(BlastResult)
lib.log.info("found %i preliminary alignments" % (len(HitList)))

#split genome fasta into individual scaffolds
if not os.path.exists(tmpdir):
    os.makedirs(tmpdir)
with open(os.path.abspath(args.genome), 'rU') as input:
    for record in SeqIO.parse(input, "fasta"):
        SeqIO.write(record, os.path.join(tmpdir, record.id + ".fa"), "fasta")

#Now run exonerate on hits
lib.log.info("Polishing alignments with Exonerate")
record_dict = SeqIO.index(os.path.abspath(args.proteins), 'fasta') #do index here in case memory problems?

#run multiprocessing exonerate
lib.runMultiProgress(runExonerate, HitList, args.cpus)
lib.log.info("Exonerate finished")

#now collect all exonerate results into one
with open(args.out, 'wb') as output:
    for root, dirs, files in os.walk(tmpdir):
        for file in files:
            if file.endswith('.out'):
                filename = os.path.join(root, file)
                with open(filename, 'rU') as readfile:
                    for line in itertools.islice(readfile, 3, None):
                        output.write(line)

#finally clean-up your mess
shutil.rmtree(tmpdir)
