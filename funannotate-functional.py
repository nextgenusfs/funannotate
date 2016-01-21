#!/usr/bin/env python

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib


#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-functional.py', usage="%(prog)s [options] -i genome.fasta -g genome.gff -o test -e youremail@mail.edu",
    description='''Script that adds functional annotation to a genome.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Annotated genome in GenBank format')
parser.add_argument('-g','--gff', required=True, help='GFF file from funannotate-predict')
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

#get executable for runiprscan
PATH2JAR = os.path.join(currentdir, 'util', 'RunIprScan-1.1.0', 'RunIprScan.jar')
RUNIPRSCAN_PATH = os.path.join(currentdir, 'util', 'RunIprScan-1.1.0')
IPROUT = os.path.join(args.out, 'iprscan')

#check dependencies
programs = ['hmmscan','blastp','gag.py', 'java']
lib.CheckDependencies(programs)

#temp exit to test code up to here
#os._exit(1)


#create temp folder to house intermediate files
if not os.path.exists(args.out):
    os.makedirs(args.out)
else:
    lib.log.error("Output directory %s already exists, will use any existing data.  If this is not what you want, exit, and provide a unique name for output folder" % (args.out))


#need to do some checks here of the input
if not args.input.endswith('.gbk' or '.gb'):
    lib.log.error("Input does not appear to be a Genbank file (it does not end in .gbk or .gb) can't run functional annotation.")
    shutil.rmtree(args.out)
    os._exit(1)
else:
    Scaffolds = os.path.join(args.out, 'genome.scaffolds.fasta')
    Proteins = os.path.join(args.out, 'genome.proteins.fasta')
    Transcripts = os.path.join(args.out, 'genome.transcripts.fasta')
    GFF = args.gff
    if not os.path.isfile(Proteins) or not os.path.isfile(Transcripts) or not os.path.isfile(Scaffolds):
        lib.gb2output(args.input, Proteins, Transcripts, Scaffolds)
        #lib.convert_to_GFF3(args.input, GFF)
    #get absolute path for all so no confusion
    for i in Scaffolds, Proteins, Transcripts, GFF:
        i = os.path.abspath(i)

ProtCount = lib.countfasta(Proteins)
lib.log.info("Loading %i protein records" % ProtCount)   

#run interpro scan, in background hopefully....
if not os.path.exists(os.path.join(args.out, 'iprscan')):
    os.makedirs(os.path.join(args.out, 'iprscan'))

if not args.iprscan: #here run the routine of IPRscan in the background
    #keep track of number of times you launched RunIprScan, want to be at least 3 times to make sure you capture all output.
    IPRcount = 0
    lib.log.info("Starting RunIprScan and running in background")
    p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    IPRcount += 1
    #while RunIprScan is running in background, run more functional annotation methods
    while p.poll() is None:
        #run PFAM-A search
        lib.log.info("Running HMMer search of PFAM domains")
        pfam_results = os.path.join(args.out, 'annotations.pfam.txt')
        if not os.path.isfile(pfam_results):
            lib.PFAMsearch(Proteins, args.cpus, 1e-50, args.out, pfam_results)
        num_annotations = lib.line_count(pfam_results)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:   #run it again to recover any that did not work
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            IPRcount +=1
            p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run SwissProt Blast search
        lib.log.info("Running Blastp search of UniProt DB")
        blast_out = os.path.join(args.out, 'annotations.swissprot.txt')
        if not os.path.isfile(blast_out):
            lib.SwissProtBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
        num_annotations = lib.line_count(blast_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            IPRcount +=1
            p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run MEROPS Blast search
        lib.log.info("Running Blastp search of MEROPS protease DB")
        blast_out = os.path.join(args.out, 'annotations.merops.txt')
        if not os.path.isfile(blast_out):
            lib.MEROPSBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
        num_annotations = lib.line_count(blast_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            IPRcount +=1
            p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run EggNog search
        eggnog_out = os.path.join(args.out, 'annotations.eggnog.txt')
        lib.log.info("Annotating proteins with EggNog 4.5 database")
        if not os.path.isfile(eggnog_out):
            lib.runEggNog(Proteins, args.cpus, 1e-10, args.out, eggnog_out)
        num_annotations = lib.line_count(eggnog_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, moving onto next process")
        else:
            lib.log.info("RunIprScan finished, but will try again to recover all results")
            IPRcount +=1
            p = subprocess.Popen(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        #run dbCAN search
        dbCAN_out = os.path.join(args.out, 'annotations.dbCAN.txt')
        lib.log.info("Annotating CAZYmes using dbCAN")
        if not os.path.isfile(dbCAN_out):
            lib.dbCANsearch(Proteins, args.cpus, 1e-17, args.out, dbCAN_out)
        num_annotations = lib.line_count(dbCAN_out)
        lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
        if p.poll() is None:
            lib.log.info("RunIprScan still running, now waiting until it finishes")
        p.wait()

    #if RunIprScan has not been run at least 3 times, run again, this time just wait for it to finish
    if IPRcount < 3:
        lib.log.info("RunIprScan has been called less than 3 times, running again")
        subprocess.call(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout = FNULL, stderr = FNULL)
    
    #lets do one final check of IPRresults, i.e. the number of proteins should equal number of files in iprscan folder
    num_prots = lib.countfasta(Proteins)
    num_ipr = len([name for name in os.listdir(IPROUT) if os.path.isfile(os.path.join(IPROUT, name))])
    if num_prots != num_ipr:
        lib.log.info("Number of IPR xml files (%i) does not equal number of proteins (%i), will try RunIprScan once more" % (num_ipr, num_prots))
        subprocess.call(['java', '-jar', PATH2JAR, '$@', '-i', Proteins, '-m', args.email, '-o', IPROUT], stdout = FNULL, stderr = FNULL)
else:   
    #check that remaining searches have been done, if not, do them.
    #run PFAM-A search
    lib.log.info("Running HMMer search of PFAM domains")
    pfam_results = os.path.join(args.out, 'annotations.pfam.txt')
    if not os.path.isfile(pfam_results):
        lib.PFAMsearch(Proteins, args.cpus, 1e-50, args.out, pfam_results)
    num_annotations = lib.line_count(pfam_results)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    #run SwissProt Blast search
    lib.log.info("Running Blastp search of UniProt DB")
    blast_out = os.path.join(args.out, 'annotations.swissprot.txt')
    if not os.path.isfile(blast_out):
        lib.SwissProtBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
    num_annotations = lib.line_count(blast_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    #run MEROPS Blast search
    lib.log.info("Running Blastp search of MEROPS protease DB")
    blast_out = os.path.join(args.out, 'annotations.merops.txt')
    if not os.path.isfile(blast_out):
        lib.MEROPSBlast(Proteins, args.cpus, 1e-5, args.out, blast_out)
    num_annotations = lib.line_count(blast_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    #run EggNog search
    eggnog_out = os.path.join(args.out, 'annotations.eggnog.txt')
    lib.log.info("Annotating proteins with EggNog 4.5 database")
    if not os.path.isfile(eggnog_out):
        lib.runEggNog(Proteins, args.cpus, 1e-10, args.out, eggnog_out)
    num_annotations = lib.line_count(eggnog_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    #run dbCAN search
    dbCAN_out = os.path.join(args.out, 'annotations.dbCAN.txt')
    lib.log.info("Annotating CAZYmes using dbCAN")
    if not os.path.isfile(dbCAN_out):
        lib.dbCANsearch(Proteins, args.cpus, 1e-17, args.out, dbCAN_out)
    num_annotations = lib.line_count(dbCAN_out)
    lib.log.info('{0:,}'.format(num_annotations) + ' annotations added')
    
#now collect the results from InterProscan, then start to reformat results
lib.log.info("RunIprScan has finished, now pulling out annotations from results")
IPR_terms = os.path.join(args.out, 'annotations.iprscan.txt')
if not os.path.isfile(IPR_terms):
    IPR2TSV = os.path.join(RUNIPRSCAN_PATH, 'ipr2tsv.py')
    with open(IPR_terms, 'w') as output:
        subprocess.call([sys.executable, IPR2TSV, IPROUT], stdout = output, stderr = FNULL)
GO_terms = os.path.join(args.out, 'annotations.GO.txt')
if not os.path.isfile(GO_terms):
    IPR2GO = os.path.join(RUNIPRSCAN_PATH, 'ipr2go.py')
    OBO = os.path.join(currentdir, 'DB', 'go.obo')
    with open(IPR_terms, 'w') as output:
        subprocess.call([sys.executable, IPR2GO, OBO, IPROUT], stdout = output, stderr = FNULL)

#now bring all annotations together and annotated genome using gag
ANNOTS = os.path.join(args.out, 'all.annotations.txt')
with open(ANNOTS, 'w') as output:
    for file in os.listdir(args.out):
        if file.startswith('annotations'):
            file = os.path.join(args.out, file)
            with open(file) as input:
                output.write(input.read())
ANNOTS = os.path.abspath(ANNOTS)
#launch gag
GAG = os.path.join(args.out, 'gag')
lib.log.info("Adding annotations to GFF using GAG")
subprocess.call(['gag.py', '-f', Scaffolds, '-g', GFF, '-a', ANNOTS, '-o', GAG], stdout = FNULL, stderr = FNULL)
os._exit(1)
    

