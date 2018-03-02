#!/usr/bin/env python

import sys, os, inspect, shutil, argparse
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib
from natsort import natsorted

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='updateGBK.py', usage="%(prog)s [options] -f genome.GBK -t genome.tbl",
    description = '''Script will update annotation of a Genbank file with new tbl.''',
    epilog = """Written by Jon Palmer (2017) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--input', required=True, help='Genome in GBK format')
parser.add_argument('-t', '--tbl', required=True, help='Genome annotation in NCBI tbl format')
parser.add_argument('-d', '--drop', help='List of locus_tag to remove/drop from annotation')
parser.add_argument('-o', '--out', help='Basename of output files')
parser.add_argument('--tbl2asn', default='-l paired-ends', help='Parameters for tbl2asn, linkage and gap info')
args=parser.parse_args()

#create log file
log_name = 'funannotate-fix.log'
if os.path.isfile(log_name):
    os.remove(log_name)
    
#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.SystemInfo()

#get version of funannotate
version = lib.get_version()
lib.log.info("Running %s" % version)

#create output and temporary directory
if args.out:
    basedir = args.out
else:
    #get location from tbl file
    basedir = os.path.dirname(args.tbl)
    
if not os.path.isdir(basedir):
    os.makedirs(basedir)
if not os.path.isdir(os.path.join(basedir, 'tbl2asn')):
    os.makedirs(os.path.join(basedir, 'tbl2asn'))

#copy over the annotation file to tbl2asn folder, or process if args.drop passed
if args.drop:
    lib.tblfilter(args.tbl, args.drop, os.path.join(basedir, 'tbl2asn', 'genome.tbl'))
else:
    shutil.copyfile(args.tbl, os.path.join(basedir, 'tbl2asn', 'genome.tbl'))

#get information info from GBK file
organism, strain, isolate, accession, WGS_accession, gb_gi, version = lib.getGBKinfo(args.input)
locustag, genenum, justify = lib.getGBKLocusTag(args.input)
if strain:
    organism_name = organism+'_'+strain
elif isolate:
    organism_name = organism+'_'+isolate
else:
    organism_name = organism
organism_name = organism_name.replace(' ', '_')

#extract fasta file from genbank file,
lib.log.info('Extracting genome sequence and parsing meta information')
contigs, genes, trnas = lib.countGenBank(args.input)
lib.log.info('{:,} contigs containing {:,} protein coding genes and {:,} tRNA genes'.format(contigs,genes,trnas))
lib.gb2dna(args.input, os.path.join(basedir, 'tbl2asn', 'genome.fsa'))

#assuming that this is the predict_results dir or update_results dir, but check first and then archive
if '_results' in basedir:
    archivedir = os.path.join(basedir, 'archive_'+str(os.getpid()))
    lib.log.info('Found pre-existing funannotate files, archiving to %s' % archivedir)
    os.makedirs(archivedir)
    #move files in results to archive dir
    for file in os.listdir(basedir):
        if 'pasa-reannotation' in file or 'WGS_accession' in file or 'ncbi.p2g' in file:
            continue
        if os.path.isfile(os.path.join(basedir, file)):
            os.rename(os.path.join(basedir, file), os.path.join(archivedir, file))

#now we can run tbl2asn
SBT = os.path.join(parentdir, 'lib', 'test.sbt')
discrep = os.path.join(basedir, organism_name+'.discrepency.txt')
if not version:
    version = 1
lib.log.info('Converting to GenBank format')
tbl2asn_cmd = lib.runtbl2asn(os.path.join(basedir, 'tbl2asn'), SBT, discrep, organism, isolate, strain, args.tbl2asn, version)

#now get GBK files from folder
lib.log.info('Generating output files.')
#setup final output files
final_fasta = os.path.join(basedir, organism_name + '.scaffolds.fa')
final_gff = os.path.join(basedir, organism_name + '.gff3')
final_gbk = os.path.join(basedir, organism_name + '.gbk')
final_tbl = os.path.join(basedir, organism_name + '.tbl')
final_proteins = os.path.join(basedir, organism_name + '.proteins.fa')
final_transcripts = os.path.join(basedir, organism_name + '.transcripts.fa')
final_validation = os.path.join(basedir, organism_name+'.validation.txt')
final_error = os.path.join(basedir, organism_name+'.error.summary.txt')
final_fixes = os.path.join(basedir, organism_name+'.models-need-fixing.txt')

#retrieve files/reorganize
shutil.copyfile(os.path.join(basedir, 'tbl2asn', 'genome.gbf'), final_gbk)
shutil.copyfile(os.path.join(basedir, 'tbl2asn', 'genome.tbl'), final_tbl)
shutil.copyfile(os.path.join(basedir, 'tbl2asn', 'genome.val'), final_validation)
shutil.copyfile(os.path.join(basedir, 'tbl2asn', 'errorsummary.val'), final_error)
lib.gb2allout(final_gbk, final_gff, final_proteins, final_transcripts, final_fasta)
errors = lib.ncbiCheckErrors(final_error, final_validation, locustag, final_fixes)
if errors > 0:
    lib.log.info("Manually edit the tbl file %s, then run:\n\nfunannotate fix -i %s -t %s\n" % (final_tbl, final_gbk, final_tbl))
else:
    contigs, genes, trnas = lib.countGenBank(final_gbk)
    lib.log.info('Output genome consists of: {:,} contigs containing {:,} protein coding genes and {:,} tRNA genes'.format(contigs,genes,trnas))

#clean up
shutil.rmtree(os.path.join(basedir, 'tbl2asn'))
