#!/usr/bin/env python

import sys, os, subprocess, shutil, inspect, itertools, argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser
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
parser.add_argument('-t','--tblastn_out', help='Save tblastn output')
parser.add_argument('--maxintron', default = 3000, help='Maximum intron size')
parser.add_argument('--logfile', default ='funannotate-p2g.log', help='logfile')
parser.add_argument('--ploidy', default =1, type=int, help='Ploidy of assembly')
parser.add_argument('--debug', action='store_true', help='Keep intermediate folders if error detected')
parser.add_argument('-f','--filter', required=True, default='tblastn', choices=['diamond', 'tblastn'], help='Method to use for pre-filter for exonerate')
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
if args.filter == 'diamond':
    diamond_version = subprocess.Popen(['diamond', '--version'], stdout=subprocess.PIPE).communicate()[0].split('\n')[0]
    diamond_version = diamond_version.split('version ')[-1]

def runDiamond(input, query, cpus, output):
    #create DB of protein sequences
    cmd = ['diamond', 'makedb', '--threads', str(cpus), '--in', query, '--db', 'diamond']
    lib.runSubprocess(cmd, output, lib.log)
    #now run search
    cmd = ['diamond', 'blastx', '--threads', str(cpus), '-q', input, '--db', 'diamond', '-o', 'diamond.matches.tab', '-e', '1e-10', '-k', '0', '--more-sensitive', '-f', '6', 'sseqid', 'slen', 'sstart', 'send', 'qseqid', 'qlen', 'qstart', 'qend', 'pident', 'length', 'evalue', 'score', 'qcovhsp', 'qframe']
    lib.runSubprocess(cmd, output, lib.log)
    
def runtblastn(input, query, cpus, output, maxhits):
    #start by formatting blast db/dustmasker filtered format
    cmd = ['dustmasker', '-in', input, '-infmt', 'fasta', '-parse_seqids', '-outfmt', 'maskinfo_asn1_bin', '-out', 'genome_dust.asnb']
    lib.runSubprocess(cmd, output, lib.log)
    cmd = ['makeblastdb', '-in', input, '-dbtype', 'nucl', '-parse_seqids', '-mask_data', 'genome_dust.asnb', '-out', 'genome']
    lib.runSubprocess(cmd, output, lib.log)
    cmd = ['tblastn', '-num_threads', str(cpus), '-db', 'genome', '-query', query, '-max_target_seqs', str(maxhits), '-db_soft_mask', '11', '-threshold', '999', '-max_intron_length', str(args.maxintron), '-evalue', '1e-10', '-outfmt', '6', '-out', 'filter.tblastn.tab']
    lib.runSubprocess(cmd, output, lib.log)

def parseDiamond(blastresult):
    Results = {}
    with open(blastresult, 'rU') as input:
        for line in input:
            cols = line.split('\t')
            hit = cols[0] + ':::' + cols[4]
            if int(cols[13]) > 0:
                start = cols[6]
                end = cols[7]
            else:
                start = cols[7]
                end = cols[6]
            start_extend = int(cols[2]) - 1 *3
            end_extend = int(cols[1]) - int(cols[3]) *3
            start = int(start) - start_extend
            if start < 1:
                start = 1
            end = int(end) + end_extend
            if end > int(cols[5]):
                end = int(cols[5])
            Results[hit] = (start, end)
    #convert Dictionary to a list that has  hit:::scaffold:::start:::stop
    HitList = []
    for k,v in Results.items():
        finalhit = k+':::'+str(v[0])+':::'+str(v[1])
        HitList.append(finalhit)
    return HitList            
            
def parseBlast(blastresult):
    Results = {}
    with open(blastresult, 'rU') as input:
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
                #get old start stop
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
    #convert Dictionary to a list that has  hit:::scaffold:::start:::stop
    HitList = []
    for k,v in Results.items():
        finalhit = k+':::'+str(v[0])+':::'+str(v[1])
        HitList.append(finalhit)
    return HitList

def runExonerate(input):
    s = input.split(':::')
    ProtID = s[0]
    ScaffID = s[1]
    ScaffStart = int(s[2])
    ScaffEnd = int(s[3])
    #get the protein model
    query = os.path.join(tmpdir, ProtID+'.'+str(os.getpid())+'.fa')
    with open(query, 'w') as output:
        SeqIO.write(protein_dict[ProtID], output, 'fasta')
    #now get the genome region, use different variable names for SeqRecords to avoid collision
    scaffold = ScaffID+'.'+ProtID+'.'+str(ScaffStart)+'-'+str(ScaffEnd)+'.fa'
    scaffold = os.path.join(tmpdir, scaffold)
    with open(scaffold, 'w') as output2:
        with open(os.path.join(tmpdir, 'scaffolds', ScaffID+'.fa'), 'rU') as fullscaff:
            for header, Sequence in SimpleFastaParser(fullscaff):
                #grab a 3 kb cushion on either side of hit region, careful of scaffold ends      
                start = ScaffStart - 3000
                if start < 1:
                    start = 1
                end = ScaffEnd + 3000
                if end > len(Sequence):
                    end = len(Sequence)
                output2.write('>%s\n%s\n' % (header, Sequence[start:end]))
    exoname = ProtID+'.'+ScaffID+'__'+str(start)+'__'
    #check that input files are created and valid
    exonerate_out = 'exonerate.' + exoname + '.out'
    exonerate_out = os.path.join(tmpdir, exonerate_out)
    ryo = "AveragePercentIdentity: %pi\n"
    cmd = ['exonerate', '--model', 'p2g', '--showvulgar', 'no', '--showalignment', 'no', '--showquerygff', 'no', '--showtargetgff', 'yes', '--maxintron', str(args.maxintron), '--percent', '80', '--ryo', ryo , query, scaffold]
    #run exonerate, capture errors
    with open(exonerate_out, 'w') as output3:
        proc = subprocess.Popen(cmd, stdout = output3, stderr=subprocess.PIPE)
    stderr = proc.communicate()
    if 'WARNING' in stderr[1]:
    	lib.log.debug('%s, Len=%i, %i-%i; %i-%i' % (header, len(Sequence), ScaffStart, ScaffEnd, start, end))
        os.rename(query, os.path.join(tmpdir, 'failed', os.path.basename(query)))
        os.rename(scaffold, os.path.join(tmpdir, 'failed', os.path.basename(scaffold)))
    else:   
        for y in [query, scaffold]:
            try:
                os.remove(y)
            except OSError:
                lib.log.debug("Error removing %s" % (y))   
    #check filesize of exonerate output, no hits still have some output data in them, should be safe dropping anything smaller than 500 bytes
    if lib.getSize(exonerate_out) < 500:
        os.remove(exonerate_out)

#count number of proteins to look for
total = lib.countfasta(args.proteins)
lib.log.info('Using {0:,}'.format(total) + ' proteins as queries')

#make tmpdir
tmpdir = 'p2g_'+ str(os.getpid())
if not os.path.isdir(tmpdir):
    os.makedirs(tmpdir)
    os.makedirs(os.path.join(tmpdir, 'failed'))
    os.makedirs(os.path.join(tmpdir, 'scaffolds'))

if args.filter == 'tblastn':
    lib.log.debug("BLAST v%s; Exonerate v%s" % (blast_version, exo_version))
    #check for tblastn input
    if args.tblastn:
        lib.log.info("Using pre-calculated tBLASTN result")
        BlastResult = args.tblastn
    else:
        lib.log.info("Running pre-filter tBLASTN step")
        BlastResult = os.path.join(tmpdir, 'filter.tblastn.tab')
        runtblastn(os.path.abspath(args.genome), os.path.abspath(args.proteins), args.cpus, tmpdir, args.ploidy*5) #2X ploidy for tBLASTn filter
    #parse results
    Hits = parseBlast(BlastResult)
else:
    lib.log.debug("Diamond v%s; Exonerate v%s" % (diamond_version, exo_version))
    #run Diamond
    lib.log.info("Running Diamond pre-filter search")
    BlastResult = os.path.join(tmpdir, 'diamond.matches.tab')
    runDiamond(os.path.abspath(args.genome), os.path.abspath(args.proteins), args.cpus, tmpdir)
    Hits = parseDiamond(BlastResult)
    
lib.log.info('Found {0:,}'.format(len(Hits))+' preliminary alignments')

#index the genome and proteins
protein_dict = SeqIO.index(os.path.abspath(args.proteins), 'fasta') #do index here in case memory problems?

#split genome fasta into individual scaffolds
with open(os.path.abspath(args.genome), 'rU') as input:
    for record in SeqIO.parse(input, "fasta"):
        SeqIO.write(record, os.path.join(tmpdir, 'scaffolds', record.id + ".fa"), "fasta")

#run multiprocessing exonerate
lib.runMultiProgress(runExonerate, Hits, args.cpus)

#now need to loop through and offset exonerate predictions back to whole scaffolds
with open(args.out, 'w') as output:
    for file in os.listdir(tmpdir):
        if file.endswith('.out'):
            with open(os.path.join(tmpdir, file), 'rU') as exoresult:
                offset = int(file.split('__')[1])
                for line in itertools.islice(exoresult, 3, None):
                    if line.startswith('#') or line.startswith('Average') or line.startswith('-- completed'):
                        output.write(line)
                    else:
                        cols = line.split('\t')
                        cols[3] = str(int(cols[3])+offset)
                        cols[4] = str(int(cols[4])+offset)
                        output.write('\t'.join(cols))

#output some quick summary of exonerate alignments that you found
Found = lib.countGFFgenes(args.out)
lib.log.info('Exonerate finished: found {0:,}'.format(Found)+' alignments')

#check for saving output of tblastn
if args.tblastn_out:
    shutil.copyfile(BlastResult, args.tblastn_out)

#finally clean-up your mess if failed is empty
if args.debug:
	try:
		os.rmdir(os.path.join(tmpdir, 'failed'))
		empty = True
	except OSError:
		empty = False
	if empty:
		shutil.rmtree(tmpdir)
	else:
		lib.log.error("Failed exonerate alignments found, see files in %s" % os.path.join(tmpdir, 'failed'))
else:
	if os.path.isfile(tmpdir):
		shutil.rmtree(tmpdir)
sys.exit(1)
