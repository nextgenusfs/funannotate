#!/usr/bin/env python

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

#get script path for directory
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-functional.py', usage="%(prog)s [options] -i genome.fasta -g genome.gff -o test -e youremail@mail.edu",
    description='''Script that adds functional annotation to a genome.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-p','--proteins', help='Protein fasta file, will get from GFF if not provided')
parser.add_argument('-g','--gff', required=True, help='GFF3 annotation of genome (if available)')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('-e','--email', required=True, help='Email address for IPRSCAN server')
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use')
parser.add_argument('--iprscan', help='Folder of pre-computed InterProScan results (1 xml per protein)')
args=parser.parse_args()

#create log file
log_name = args.out + '.funnannotate-functional.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.log.info("Operating system: %s, %i cores, %i GB RAM" % (sys.platform, multiprocessing.cpu_count(), lib.MemoryCheck()))

#check dependencies
programs = ['hmmscan','blastp','blastn','gag.py','tbl2asn','runiprscan']
lib.CheckDependencies(programs)

#create temp folder to house intermediate files
if not os.path.exists(args.out):
    os.makedirs(args.out)

#need to do some checks here of the input
if



#run interpro scan, in background hopefully....
if not os.path.exists('iprscan'):
    os.makedirs('iprscan')
#keep track of number of times you launched RunIprScan
IPRcount = 0
lib.log.info("Starting RunIprScan and running in background")
p = subprocess.Popen(['runiprscan', '-i', args.input, '-m', args.email, '-o', 'iprscan'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
IPRcount += 1
while p.poll() is None:
    #run PFAM-A search
    lib.log.info("Running HMMer search of PFAM domains")
    pfam_results = args.out + '.pfam.txt'
    lib.PFAMsearch(args.input, args.cpus, 1e-50, args.out, pfam_results)
    num_annotations = lib.line_count(pfam_results)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    if p.poll() is None:
        lib.log.info("RunIprScan still running, moving onto next process")
    else:   #run it again to recover any that did not work
        lib.log.info("RunIprScan finished, but will try again to recover all results")
        IPRcount +=1
        p = subprocess.Popen(['runiprscan', '-i', args.input, '-m', args.email, '-o', 'iprscan'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #run SwissProt Blast search
    lib.log.info("Running Blastp search of UniProt DB")
    blast_out = args.out + '.swissprot.txt'
    lib.SwissProtBlast(args.input, args.cpus, 1e-5, args.out, blast_out)
    num_annotations = lib.line_count(blast_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    if p.poll() is None:
        lib.log.info("RunIprScan still running, moving onto next process")
    else:
        lib.log.info("RunIprScan finished, but will try again to recover all results")
        IPRcount +=1
        p = subprocess.Popen(['runiprscan', '-i', args.input, '-m', args.email, '-o', 'iprscan'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #run MEROPS Blast search
    lib.log.info("Running Blastp search of MEROPS protease DB")
    blast_out = args.out + '.merops.txt'
    lib.MEROPSBlast(args.input, args.cpus, 1e-5, args.out, blast_out)
    num_annotations = lib.line_count(blast_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    if p.poll() is None:
        lib.log.info("RunIprScan still running, moving onto next process")
    else:
        lib.log.info("RunIprScan finished, but will try again to recover all results")
        IPRcount +=1
        p = subprocess.Popen(['runiprscan', '-i', args.input, '-m', args.email, '-o', 'iprscan'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #run EggNog search
    eggnog_out = args.out + '.eggnog.txt'
    lib.log.info("Annotating proteins with EggNog 4.5 database")
    lib.runEggNog(args.input, args.cpus, 1e-10, args.out, eggnog_out)
    num_annotations = lib.line_count(eggnog_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    if p.poll() is None:
        lib.log.info("RunIprScan still running, moving onto next process")
    else:
        lib.log.info("RunIprScan finished, but will try again to recover all results")
        IPRcount +=1
        p = subprocess.Popen(['runiprscan', '-i', args.input, '-m', args.email, '-o', 'iprscan'], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    #run dbCAN search
    dbCAN_out = args.out + '.dbCAN.txt'
    lib.log.info("Annotating CAZYmes using dbCAN")
    lib.dbCANsearch(args.input, args.cpus, 1e-17, args.out, dbCAN_out)
    num_annotations = lib.line_count(dbCAN_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    if p.poll() is None:
        lib.log.info("RunIprScan still running, now waiting until it finishes")
    p.wait()

#if RunIprScan has not been run at least 3 times, run again, this time just wait for it to finish
if IPRcount < 3:
    lib.log.info("RunIprScan has been called less than 3 times, running again")
    subprocess.call(['runiprscan', '-i', args.input, '-m', args.email, '-o', 'iprscan'], stdout = FNULL, stderr = FNULL)
lib.log.info("RunIprScan has finished, now pulling out annotations from results")

#now collect the results from InterProscan
