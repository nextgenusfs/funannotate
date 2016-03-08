#!/usr/bin/env python
from __future__ import division

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time, csv, glob
from Bio import SeqIO
from natsort import natsorted
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-fix.py',
    description='''Script that parses NCBI contamination report and fixes your genome..''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Folder from funannotate.')
parser.add_argument('-c','--contamination', help='Contamination report from NCBI.')
parser.add_argument('--sbt', default='SBT', help='Basename of output files')
parser.add_argument('-s','--species', help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('--isolate', help='Isolate/strain name (e.g. Af293)')
args=parser.parse_args()

#create log file
log_name = os.path.join(args.input, 'logfiles', 'funnannotate-fix.log')
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
programs = ['gag.py']
lib.CheckDependencies(programs)

#take care of some preliminary checks
if args.sbt == 'SBT':
    SBT = os.path.join(parentdir, 'lib', 'test.sbt')
    lib.log.info("No NCBI SBT file given, will use default, however if you plan to submit to NCBI, create one and pass it here '--sbt'")
else:
    SBT = args.sbt

#need to do some checks here of the input
#should be a folder, with funannotate files, thus store results there, no need to create output folder
if not os.path.isdir(args.input):
    lib.log.error("%i directory does not exist" % args.input)
    os._exit(1)
if os.path.isdir(os.path.join(args.input, 'predict_results')): #funannotate results should be here
    inputdir = os.path.join(args.input, 'predict_results')
    outputdir = args.input
else:
    inputdir = os.path.join(args.input) #here user specified the predict_results folder, or it is all together wrong, find out in next few lines
    outputdir = lib.get_parent_dir(args.input)
#get files that you need
for file in os.listdir(inputdir):
    if file.endswith('.scaffolds.fa'):
        Scaffolds = os.path.join(inputdir, file)
    if file.endswith('.gff3'):
        GFF = os.path.join(inputdir, file)

#now get the AGP file - this will be in annotate_results folder
if os.path.isdir(os.path.join(outputdir, 'annotate_results')):
    for file in os.listdir(os.path.join(outputdir, 'annotate_results')):
        if file.endswith('.agp'):
            AGP = os.path.join(outputdir, 'annotate_results', file)
        if file.endswith('.gbk'):
            GBK = os.path.join(outputdir, 'annotate_results', file)

#get absolute path for all input so there are no problems later
for i in Scaffolds, GFF, GBK, AGP:
    i = os.path.abspath(i)
        
#get organism and isolate from GBK file
if not args.species:
    with open(GBK, 'rU') as gbk:
        SeqRecords = SeqIO.parse(gbk, 'genbank')
        for record in SeqRecords:
            for f in record.features:
                if f.type == "source":
                    organism = f.qualifiers.get("organism", ["???"])[0]
                    if not args.isolate:
                        isolate = f.qualifiers.get("isolate", ["???"])[0]
                    else:
                        isolate = args.isolate
                    break
else:
    organism = args.species
    if not args.isolate:
        isolate = '???'
    else:
        isolate = args.isolate

############################################################################

#load the AGP file into a dictionary using contig name as key
agpDict = {}
scaffDict = {}
with open(AGP, 'rU') as input:
    for line in input:
        col = line.split('\t')
        if col[4] == 'W':
            contig = col[5]
            c_start = col[6]
            c_end = col[7]
            scaffold = col[0]
            s_start = col[1]
            s_end = col[2]
            if contig not in agpDict:
                agpDict[contig] = [c_start, c_end, scaffold, s_start, s_end]
            if scaffold not in scaffDict:
                scaffDict[scaffold] = s_end
            else:
                if int(scaffDict.get(scaffold)) < int(s_end):
                    scaffDict[scaffold] = s_end
            
#parse the error report
TRIM = os.path.join(outputdir, 'annotate_misc', 'trim.bed')
with open(TRIM, 'w') as output:
    with open(args.contamination, 'rU') as input:
        for line in input:
            if line.startswith('contig_'):
                line = line.replace('..', '\t')
                col = line.split('\t')
                scaffhits = agpDict.get(col[0])
                if scaffhits[3] == '1':
                    start = int(col[2])
                    stop = int(col[3])
                else:
                    start = int(col[2]) + int(scaffhits[3]) - 1
                    stop = int(col[3]) + int(scaffhits[3]) - 1
                output.write("%s\t%i\t%i\n" % (scaffhits[2], start, stop))

ANNOTS = os.path.join(outputdir, 'annotate_misc', 'all.annotations.txt')

#launch gag
GAG = os.path.join(outputdir, 'annotate_misc', 'gag2')
lib.log.info("Adding annotations to GFF using GAG")
subprocess.call(['gag.py', '-f', Scaffolds, '-g', GFF, '-t', TRIM, '-a', ANNOTS, '-o', GAG], stdout = FNULL, stderr = FNULL)

#fix the tbl file for tRNA genes
lib.log.info("Fixing tRNA annotations in GenBank tbl file")
original = os.path.join(outputdir, 'annotate_misc','gag2', 'genome.tbl')
tmp_tbl = os.path.join(outputdir, 'annotate_misc','gag2', 'genome.tbl.original')
os.rename(original, tmp_tbl)
lib.CleantRNAtbl(GFF, tmp_tbl, original)

#write to GBK file
if not isolate == '???':
    ORGANISM = "[organism=" + organism + "] " + "[isolate=" + isolate + "]"
    baseOUTPUT = organism + '_' + isolate
else:
    ORGANISM = "[organism=" + organism + "]"
    baseOUTPUT = organism
#remove any spaces from baseoutput 
baseOUTPUT = baseOUTPUT.replace(' ', '_')

#launch tbl2asn to create genbank submission files
shutil.copyfile(os.path.join(GAG, 'genome.fasta'), os.path.join(GAG, 'genome.fsa'))
discrep = 'discrepency.report.txt'
lib.log.info("Converting to final Genbank format, good luck!.....")
subprocess.call(['tbl2asn', '-p', GAG, '-t', SBT, '-M', 'n', '-Z', discrep, '-a', 'r10u', '-l', 'paired-ends', '-j', ORGANISM, '-V', 'b', '-c', 'fx'], stdout = FNULL, stderr = FNULL)

#collected output files and rename accordingly
ResultsFolder = os.path.join(outputdir, 'annotate_results')
os.rename(discrep, os.path.join(ResultsFolder, baseOUTPUT+'.discrepency.report.txt'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag2', 'genome.gbf'), os.path.join(ResultsFolder, baseOUTPUT+'.gbk'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag2', 'genome.gff'), os.path.join(ResultsFolder, baseOUTPUT+'.gff3'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag2', 'genome.tbl'), os.path.join(ResultsFolder, baseOUTPUT+'.tbl'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag2', 'genome.sqn'), os.path.join(ResultsFolder, baseOUTPUT+'.sqn'))
os.rename(os.path.join(outputdir, 'annotate_misc', 'gag2', 'genome.fasta'), os.path.join(ResultsFolder, baseOUTPUT+'.scaffolds.fa'))
shutil.rmtree(PROTS)

#write AGP output so all files in correct directory
lib.log.info("Creating AGP file and corresponding contigs file")
agp2fasta = os.path.join(parentdir, 'util', 'fasta2agp.pl')
AGP = os.path.join(ResultsFolder, baseOUTPUT+'.agp')
with open(AGP, 'w') as output:
    subprocess.call(['perl', agp2fasta, baseOUTPUT+'.scaffolds.fa'], cwd = ResultsFolder, stdout = output, stderr = FNULL)