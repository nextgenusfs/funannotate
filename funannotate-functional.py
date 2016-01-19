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
parser.add_argument('-i','--input', required=True, help='Annotated genome in GenBank format')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('-e','--email', required=True, help='Email address for IPRSCAN server')
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use')
parser.add_argument('--iprscan', help='Folder of pre-computed InterProScan results (1 xml per protein)')
parser.add_argument('--force', action='store_true', help='Over-write output folder')
args=parser.parse_args()

#create log file
log_name = 'funnannotate-functional.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.log.info("Operating system: %s, %i cores, ~ %i GB RAM" % (sys.platform, multiprocessing.cpu_count(), lib.MemoryCheck()))

#check dependencies
programs = ['hmmscan','blastp','gag.py','runiprscan']
lib.CheckDependencies(programs)

#create temp folder to house intermediate files
if not os.path.exists(args.out):
    os.makedirs(args.out)
else:
    if not args.force:
        lib.log.error("Output directory %s already exists, provide a unique name for output folder" % (args.out))
        os._exit(1)
    else:
        shutil.rmtree(args.out)
        os.makedirs(args.out)

#need to do some checks here of the input
if not args.input.endswith('.gbk' or '.gb'):
    lib.log.error("Input does not appear to be a Genbank file (it does not end in .gbk or .gb) can't run functional annotation.")
    shutil.rmtree(args.out)
    os._exit(1)
else:
    Scaffolds = os.path.join(args.out, 'genome.scaffolds.fasta')
    Proteins = os.path.join(args.out, 'genome.proteins.fasta')
    Transcripts = os.path.join(args.out, 'genome.transcripts.fasta')
    lib.gb2output(args.input, Proteins, Transcripts, Scaffolds)
    #get absolute path for all so no confusion
    for i in Scaffolds, Proteins, Transcripts:
        i = os.path.abspath(i)


#temp exit to test code up to here
#os._exit(1)
   

#run interpro scan, in background hopefully....
if not os.path.exists(os.path.join(args.out, 'iprscan')):
    os.makedirs(os.path.join(args.out, 'iprscan'))

if not args.iprscan: #here run the routine of IPRscan in the background
    #keep track of number of times you launched RunIprScan, want to be at least 3 times to make sure you capture all output.
    IPRcount = 0
    lib.log.info("Starting RunIprScan and running in background")
    p = subprocess.Popen(['runiprscan', '-i', Proteins, '-m', args.email, '-o', 'iprscan'], cwd = args.out, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    IPRcount += 1
    #while RunIprScan is running in background, run more functional annotation methods
    while p.poll() is None:
        #run PFAM-A search
        lib.log.info("Running HMMer search of PFAM domains")
        pfam_results = os.path.join(args.out, 'annotations.pfam.txt')
        lib.PFAMsearch(Proteins, args.cpus, 1e-50, args.out, pfam_results)
        num_annotations = lib.line_count(pfam_results)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:   #run it again to recover any that did not work
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            IPRcount +=1
            p = subprocess.Popen(['runiprscan', '-i', Proteins, '-m', args.email, '-o', 'iprscan'], cwd = args.out, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run SwissProt Blast search
        lib.log.info("Running Blastp search of UniProt DB")
        blast_out = os.path.join(args.out, 'annotations.swissprot.txt')
        lib.SwissProtBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
        num_annotations = lib.line_count(blast_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            IPRcount +=1
            p = subprocess.Popen(['runiprscan', '-i', Proteins, '-m', args.email, '-o', 'iprscan'], cwd = args.out, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run MEROPS Blast search
        lib.log.info("Running Blastp search of MEROPS protease DB")
        blast_out = os.path.join(args.out, 'annotations.merops.txt')
        lib.MEROPSBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
        num_annotations = lib.line_count(blast_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            IPRcount +=1
            p = subprocess.Popen(['runiprscan', '-i', Proteins, '-m', args.email, '-o', 'iprscan'], cwd = args.out, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run EggNog search
        eggnog_out = os.path.join(args.out, 'annotations.eggnog.txt')
        lib.log.info("Annotating proteins with EggNog 4.5 database")
        lib.runEggNog(Proteins, args.cpus, 1e-10, args.out, eggnog_out)
        num_annotations = lib.line_count(eggnog_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            IPRcount +=1
            p = subprocess.Popen(['runiprscan', '-i', Proteins, '-m', args.email, '-o', 'iprscan'], cwd = args.out, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run dbCAN search
        dbCAN_out = os.path.join(args.out, 'annotations.dbCAN.txt')
        lib.log.info("Annotating CAZYmes using dbCAN")
        lib.dbCANsearch(Proteins, args.cpus, 1e-17, args.out, dbCAN_out)
        num_annotations = lib.line_count(dbCAN_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, now waiting until it finishes")
        p.wait()

    #if RunIprScan has not been run at least 3 times, run again, this time just wait for it to finish
    if IPRcount < 3:
        lib.log.info("RunIprScan has been called less than 3 times, running again")
        subprocess.call(['runiprscan', '-i', Proteins, '-m', args.email, '-o', 'iprscan'], cwd = args.out, stdout = FNULL, stderr = FNULL)
    lib.log.info("RunIprScan has finished, now pulling out annotations from results")
else:
    #now collect the results from InterProscan, then run remaining searches
    os._exit(1)
    

