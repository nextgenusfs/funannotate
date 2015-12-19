#!/usr/bin/env python

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir) 
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate.py', usage="%(prog)s [options] -i genome.fasta",
    description='''Script that does it all.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use')
args=parser.parse_args() 

#create log file
log_name = args.out + '.funannotate.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.log.info("Operating system: %s" % sys.platform)

programs = ['hmmscan','blastp','blastn','gag.py','tbl2asn','runiprscan','gmes_petap.pl', 'RepeatModeler', 'RepeatMasker']
missing = []
for p in programs:
    if lib.which(p) == False:
        missing.append(p)
if missing != []:
    error = ", ".join(missing)
    lib.log.error("Missing Dependencies: %s.  Please install missing dependencies and re-run script" % (error))
    sys.exit(1)

#create temp folder
if not os.path.exists(args.out):
    os.makedirs(args.out)
    
#Run genemark-ES on genome
gmes_log = os.path.join(args.out, 'genemark.log'
if os.path.isfile(gmes_log):
    os.remove(gmes_log)
contigCount = lib.countfasta(args.input)
lib.log.info('Loading genome assembly: ' + '{0:,}'.format(total) + ' contigs'
lib.log.info("Running GeneMark-ES on assembly")
lib.log.debug("gmes_petap.pl --ES --fungus --cores %i --sequence %s" % (args.cpus, args.input))
with open(gmes_log, 'w') as logfile:
    subprocess.call(['gmes_petap.pl', '--ES', '--fungus', '--cores', args.cpus, '--sequence', args.input], stdout = logfile, stderr = logfile)

#filter GeneMark models based on fCEGMA to get training set for augustus


#this is my attempt to collate all of the genome annotation scripts into a "pipeline"
'''Things needed to be installed:
runiprscan
hmmer3.1
blast
RepeatMasker
RepeatModeler
GeneMark-ES
Augustus
GAG
tbl2asn
SNAP?
Maker? - maybe skip this and just do GeneMark/Augustus or just GeneMark as it seems to work fairly well and is fast.
EVM

Downloads
Eggnog: http://eggnogdb.embl.de/download/eggnog_4.5/data/fuNOG/fuNOG.hmm.tar.gz
        http://eggnogdb.embl.de/download/eggnog_4.5/data/fuNOG/fuNOG.annotations.tsv.gz
dbCAN:

SwissProt:


'''

