#!/usr/bin/env python

import sys, os, subprocess, inspect, shutil, argparse, re
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
import lib.library as lib
from natsort import natsorted

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self, prog):
        super(MyFormatter, self).__init__(prog, max_help_position=48)
parser = argparse.ArgumentParser(prog='funannotate-predict.py', usage="%(prog)s [options] -i genome.fasta",
    description = '''Script that does it all.''',
    epilog = """Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i', '--input', help='Genome in FASTA format')
parser.add_argument('-o', '--out', required=True, help='Basename of output files')
parser.add_argument('-s', '--species', required=True, help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('--isolate', default=False, help='Isolate name (e.g. Af293)')
parser.add_argument('--strain', default=False, help='Strain name (e.g. CEA10)')
parser.add_argument('--masked_genome', help='Soft-masked Genome in FASTA format ')
parser.add_argument('--repeatmasker_gff3', help='RepeatMasker GFF3 file')
parser.add_argument('--repeatmasker_species', help='RepeatMasker species, will skip repeatmodeler')
parser.add_argument('--header_length', default=16, type=int, help='Max length for fasta headers')
parser.add_argument('--name', default="FUN_", help='Shortname for genes, perhaps assigned by NCBI, eg. VC83')
parser.add_argument('--augustus_species', help='Specify species for Augustus')
parser.add_argument('--genemark_mod', help='Use pre-existing Genemark training file (e.g. gmhmm.mod)')
parser.add_argument('--protein_evidence', nargs='+', help='Specify protein evidence (multiple files can be separaed by a space)')
parser.add_argument('--use_diamond', action='store_true', help='Use Diamond instead of tblastn for mapping protein evidence')
parser.add_argument('--exonerate_proteins', help='Pre-computed Exonerate protein alignments (see README for how to run exonerate)')
parser.add_argument('--transcript_evidence', nargs='+', help='Transcript evidence (map to genome with GMAP)')
parser.add_argument('--gmap_gff', help='Pre-computed GMAP transcript alignments (GFF3)')
parser.add_argument('--pasa_gff', help='Pre-computed PASA/TransDecoder high quality models')
parser.add_argument('--other_gff', help='GFF gene prediction pass-through to EVM')
parser.add_argument('--augustus_gff', help='Pre-computed Augustus gene models (GFF3)')
parser.add_argument('--genemark_gtf', help='Pre-computed GeneMark gene models (GTF)')
parser.add_argument('--maker_gff', help='MAKER2 GFF output')
parser.add_argument('--repeatmodeler_lib', help='Pre-computed RepeatModeler (or other) repetitive elements')
parser.add_argument('--rna_bam', help='BAM (sorted) of RNAseq aligned to reference for BRAKER1')
parser.add_argument('--min_intronlen', default=10, help='Minimum intron length for gene models')
parser.add_argument('--max_intronlen', default=3000, help='Maximum intron length for gene models')
parser.add_argument('--min_protlen', default=50, type=int, help='Minimum amino acid length for valid gene model')
parser.add_argument('--keep_no_stops', action='store_true', help='Keep gene models without valid stop codons')
parser.add_argument('--ploidy', default=1, type=int, help='Ploidy of assembly')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to use')
parser.add_argument('--busco_seed_species', default='anidulans', help='Augustus species to use as initial training point for BUSCO')
parser.add_argument('--optimize_augustus', action='store_true', help='Run "long" training of Augustus')
parser.add_argument('--busco_db', default='dikarya', help='BUSCO model database')
parser.add_argument('-t','--tbl2asn', default='-l paired-ends', help='Parameters for tbl2asn, linkage and gap info')
parser.add_argument('--organism', default='fungus', choices=['fungus', 'other'], help='Fungal specific settings')
parser.add_argument('-d','--database', help='Path to funannotate database, $FUNANNOTATE_DB')
parser.add_argument('--EVM_HOME', help='Path to Evidence Modeler home directory, $EVM_HOME')
parser.add_argument('--AUGUSTUS_CONFIG_PATH', help='Path to Augustus config directory, $AUGUSTUS_CONFIG_PATH')
parser.add_argument('--GENEMARK_PATH', help='Path to GeneMark exe (gmes_petap.pl) directory, $GENEMARK_PATH')
parser.add_argument('--BAMTOOLS_PATH', help='Path to BamTools exe directory, $BAMTOOLS_PATH')
args=parser.parse_args()

#check for conflicting folder names to avoid problems
conflict = ['busco', 'busco_proteins', 'RepeatMasker', 'RepeatModeler', 'genemark', 'EVM_tmp', 'braker']
if args.out in conflict:
    lib.log.error("%s output folder conflicts with a hard coded tmp folder, please change -o parameter" % args.out)
    sys.exit(1)

#create folder structure
if not os.path.isdir(args.out):
    os.makedirs(args.out)
    os.makedirs(os.path.join(args.out, 'predict_misc'))
    os.makedirs(os.path.join(args.out, 'predict_results'))
    os.makedirs(os.path.join(args.out, 'logfiles'))
else:
    if os.path.isdir(os.path.join(args.out, 'predict_results')):
        shutil.rmtree(os.path.join(args.out, 'predict_results'))
        os.makedirs(os.path.join(args.out, 'predict_results'))
    #make sure subdirectories exist
    dirs = [os.path.join(args.out, 'predict_misc'), os.path.join(args.out, 'logfiles'), os.path.join(args.out, 'predict_results')]
    for d in dirs:
        if not os.path.isdir(d):
            os.makedirs(d)
    
#create log file
log_name = os.path.join(args.out, 'logfiles', 'funannotate-predict.log')
if os.path.isfile(log_name):
    os.remove(log_name)

#create debug log file for Repeats(capture stderr)
debug = os.path.join(args.out, 'logfiles', 'funannotate-repeats.log')
if os.path.isfile(debug):
    os.remove(debug)

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

#setup funannotate DB path
if args.database:
    FUNDB = args.database
else:
    try:
        FUNDB = os.environ["FUNANNOTATE_DB"]
    except KeyError:
        FUNDB = os.path.join(parentdir,'DB')
        
#check if database setup        
blastdb = os.path.join(FUNDB,'REPEATS.psq')
if not os.path.isfile(blastdb):
    lib.log.error("Can't find Repeat Database at {:}, you may need to re-run funannotate setup or funannotate database".format(os.path.join(FUNDB, 'REPEATS')))
    sys.exit(1)
#check buscos, download if necessary
if not os.path.isdir(os.path.join(FUNDB, args.busco_db)):
    lib.download_buscos(args.busco_db)
    
#do some checks and balances
try:
    EVM = os.environ["EVM_HOME"]
except KeyError:
    if not args.EVM_HOME:
        lib.log.error("$EVM_HOME environmental variable not found, Evidence Modeler is not properly configured.  You can use the --EVM_HOME argument to specifiy a path at runtime")
        sys.exit(1)
    else:
        EVM = args.EVM_HOME

try:
    AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
except KeyError:
    if not args.AUGUSTUS_CONFIG_PATH:
        lib.log.error("$AUGUSTUS_CONFIG_PATH environmental variable not found, Augustus is not properly configured. You can use the --AUGUSTUS_CONFIG_PATH argument to specify a path at runtime.")
        sys.exit(1)
    else:
        AUGUSTUS = args.AUGUSTUS_CONFIG_PATH
        
#if you want to use BRAKER1, you also need some additional config paths
try:
    GENEMARK_PATH = os.environ["GENEMARK_PATH"]
except KeyError:
    if not lib.which('gmes_petap.pl'):
        if not args.GENEMARK_PATH:
            lib.log.error("GeneMark not found and $GENEMARK_PATH environmental variable missing, BRAKER1 is not properly configured. You can use the --GENEMARK_PATH argument to specify a path at runtime.")
            sys.exit(1)
        else:
            GENEMARK_PATH = args.GENEMARK_PATH

try:
    BAMTOOLS_PATH = os.environ["BAMTOOLS_PATH"]
except KeyError:
    #check if it is in PATH, if it is, no problem, else through warning
    if not lib.which('bamtools'):
        if not args.BAMTOOLS_PATH:
            lib.log.error("Bamtools not found and $BAMTOOLS_PATH environmental variable missing, BRAKER1 is not properly configured. You can use the --BAMTOOLS_PATH argument to specify a path at runtime.")
            sys.exit(1)
        else:
            BAMTOOLS_PATH = args.BAMTOOLS_PATH

if AUGUSTUS.endswith('config'):
    AUGUSTUS_BASE = AUGUSTUS.replace('config', '')
elif AUGUSTUS.endswith('config'+os.sep):
    AUGUSTUS_BASE = AUGUSTUS.replace('config'+os.sep, '')
AutoAug = os.path.join(AUGUSTUS_BASE, 'scripts', 'autoAug.pl')
GeneMark2GFF = os.path.join(parentdir, 'util', 'genemark_gtf2gff3.pl')

programs = ['tblastn', 'exonerate', 'makeblastdb', 'dustmasker', 'gag.py', 'tbl2asn', 'gmes_petap.pl', 'BuildDatabase', 'RepeatModeler', 'RepeatMasker', GeneMark2GFF, AutoAug, 'bedtools', 'gmap', 'gmap_build', 'blat', 'pslCDnaFilter', 'augustus', 'etraining', 'rmOutToGFF3.pl']
lib.CheckDependencies(programs)
if args.use_diamond:
    lib.CheckDependencies(['diamond'])

#check augustus species now, so that you don't get through script and then find out it is already in DB
if not args.augustus_species:
    aug_species = args.species.replace(' ', '_').lower()
else:
    aug_species = args.augustus_species
augspeciescheck = lib.CheckAugustusSpecies(aug_species)
if augspeciescheck and not args.maker_gff:
    lib.log.error("Augustus training set for %s already exists, thus funannotate will use those parameters.\n\t\tIf you want to re-train, provide a unique name for the --augustus_species argument" % (aug_species))

#check augustus functionality
augustuscheck = lib.checkAugustusFunc(AUGUSTUS_BASE)
system_os = lib.systemOS()
if args.rna_bam:
    if augustuscheck[1] == 0:
        lib.log.error("ERROR: %s is not installed properly for BRAKER1 (check bam2hints compilation)" % augustuscheck[0])
        sys.exit(1)
if not augspeciescheck: #means training needs to be done
    if augustuscheck[2] == 0:
        if 'MacOSX' in system_os:
            lib.log.error("ERROR: %s is not installed properly and this version not work with BUSCO, on %s you should try manual compilation with gcc-6 of v3.2.1." % (augustuscheck[0], system_os))
        elif 'Ubuntu' in system_os:
            lib.log.error("ERROR: %s is not installed properly and this version not work with BUSCO, on %s you should install like this: `brew install augustus`." % (augustuscheck[0], system_os))
        elif 'centos' in system_os:
            lib.log.error("ERROR: %s is not installed properly and this version not work with BUSCO, on %s you may need a much older version ~ v3.03." % (augustuscheck[0], system_os))
        else:
            lib.log.error("ERROR: %s is not installed properly and this version not work with BUSCO, this is a problem with Augustus compliatation, you may need to compile manually on %s." % (augustuscheck[0], system_os))
        if not args.pasa_gff: #first training will use pasa, otherwise BUSCO
            sys.exit(1)
        else:
            lib.log.info("Will proceed with PASA models to train Augustus")
            
#if made it here output Augustus version
lib.log.info("%s detected, version seems to be compatible with BRAKER1 and BUSCO" % augustuscheck[0])

#check input files to make sure they are not empty, first check if multiple files passed to transcript/protein evidence
input_checks = [args.input, args.masked_genome, args.repeatmasker_gff3, args.genemark_mod, args.exonerate_proteins, args.gmap_gff, args.repeatmodeler_lib, args.pasa_gff, args.other_gff, args.rna_bam]
if not args.protein_evidence:
    args.protein_evidence = os.path.join(FUNDB, 'uniprot_sprot.fasta')
input_checks = input_checks + args.protein_evidence
if args.transcript_evidence:  #if transcripts passed, otherwise ignore
    input_checks = input_checks + args.transcript_evidence
#now check the inputs
for i in input_checks:
    if i:
        lib.checkinputs(i)

#convert PASA GFF and/or GFF pass-through
#convert PASA to have 'pasa_pred' in second column to make sure weights work with EVM
PASA_GFF = os.path.join(args.out, 'predict_misc', 'pasa_predictions.gff3')
PASA_weight = '10'
if args.pasa_gff:
    if ':' in args.pasa_gff:
        pasacol = args.pasa_gff.split(':')
        PASA_weight = pasacol[1]
        args.pasa_gff = pasacol[0]
    lib.renameGFF(os.path.abspath(args.pasa_gff), 'pasa_pred', PASA_GFF)
    #validate it will work with EVM
    if not lib.evmGFFvalidate(PASA_GFF, EVM, lib.log):
        lib.log.error("ERROR: %s is not a properly formatted PASA GFF file, please consult EvidenceModeler docs" % args.pasa_gff)
        sys.exit(1)
OTHER_GFF = os.path.join(args.out, 'predict_misc', 'other_predictions.gff3')
OTHER_weight = '1'
if args.other_gff:
    if ':' in args.other_gff:
        othercol = args.other_gff.split(':')
        OTHER_weight = othercol[1]
        args.other_gff = othercol[0]
    lib.renameGFF(os.path.abspath(args.other_gff), 'other_pred', OTHER_GFF)
    #validate it will work with EVM
    if not lib.evmGFFvalidate(OTHER_GFF, EVM, lib.log):
        lib.log.error("ERROR: %s is not a properly formatted GFF file, please consult EvidenceModeler docs" % args.other_gff)
        sys.exit(1)

#setup the genome fasta file, need either args.input or need to have args.masked_genome + args.repeatmasker_gff3
#declare output location
genome_input = os.path.join(args.out, 'predict_misc', 'genome.fasta')
MaskGenome = os.path.join(args.out, 'predict_misc', 'genome.softmasked.fa')
RepeatMasker = os.path.join(args.out, 'predict_misc', 'repeatmasker.gff3')
Scaffoldsort = os.path.join(args.out, 'predict_misc', 'scaffold.sort.order.txt')
Renamingsort = os.path.join(args.out, 'predict_misc', 'scaffold.sort.rename.txt')
#check inputs
if args.input:
    #check fasta header length
    header_test = lib.checkFastaHeaders(args.input, args.header_length)
    if not header_test[0]:
        lib.log.error("Fasta headers on your input have more characters than the max (%i), reformat headers to continue." % args.header_length)
        lib.log.error("First 5 header names:\n%s" % '\n'.join(header_test[1][:5]))
        sys.exit(1)
    else:
        with open(Scaffoldsort, 'w') as contigsout:
            sortedHeaders = natsorted(header_test[1])
            contigsout.write('%s' % '\n'.join(sortedHeaders))
        with open(Renamingsort, 'w') as renameout:
            counter = 0
            with open(Scaffoldsort, 'rU') as contigsin:
                for line in contigsin:
                    counter +=1
                    line = line.replace('\n', '') 
                    renameout.write('%s\t%i\n' % (line, counter))
                    
    #if BAM file passed, check if headers are same as input
    if args.rna_bam:
        if not lib.BamHeaderTest(args.input, args.rna_bam):
            lib.log.error("Fasta headers in BAM file do not match genome, exiting.")
            sys.exit(1)
    #just copy the input fasta to the misc folder and move on.
    shutil.copyfile(args.input, genome_input)
    Genome = os.path.abspath(genome_input)
else:
    if not args.masked_genome or not args.repeatmasker_gff3:
        lib.log.error("Input Error: either need -i,--input or need both --masked_genome and --repeatmasker_gff3")
        sys.exit(1)
    for x in [args.masked_genome, args.repeatmasker_gff3]:
        lib.checkinputs(x)
    #if BAM file passed, check if headers are same as input
    if args.rna_bam:
        if not lib.BamHeaderTest(args.masked_genome, args.rna_bam):
            lib.log.error("Fasta headers in BAM file do not match genome, exiting.")
            sys.exit(1)
    #now copy over masked genome and repeatmasker
    shutil.copyfile(args.masked_genome, genome_input)
    shutil.copyfile(args.masked_genome, MaskGenome)
    shutil.copyfile(args.repeatmasker_gff3, RepeatMasker)
    
#setup augustus parallel command
AUGUSTUS_PARALELL = os.path.join(parentdir, 'bin', 'augustus_parallel.py')

#EVM command line scripts
Converter = os.path.join(EVM, 'EvmUtils', 'misc', 'augustus_GFF3_to_EVM_GFF3.pl')
ExoConverter = os.path.join(EVM, 'EvmUtils', 'misc', 'exonerate_gff_to_alignment_gff3.pl')
Converter2 = os.path.join(EVM, 'EvmUtils', 'misc', 'augustus_GTF_to_EVM_GFF3.pl')
EVM2proteins = os.path.join(EVM, 'EvmUtils', 'gff3_file_to_proteins.pl')
    
#repeatmasker, run if not passed from command line
if not os.path.isfile(MaskGenome):
    if not args.repeatmodeler_lib:
        if not args.repeatmasker_species:
            lib.RepeatModelMask(Genome, args.cpus, os.path.join(args.out, 'predict_misc'), MaskGenome, debug)
        else:
            lib.RepeatMaskSpecies(Genome, args.repeatmasker_species, args.cpus, os.path.join(args.out, 'predict_misc'), MaskGenome, debug)
    else:
        lib.RepeatMask(Genome, args.repeatmodeler_lib, args.cpus, os.path.join(args.out, 'predict_misc'), MaskGenome, debug)
#make sure absolute path
RepeatMasker = os.path.abspath(RepeatMasker)
MaskGenome = os.path.abspath(MaskGenome)
#final output for augustus hints, declare ahead of time for checking portion of script
hintsE = os.path.join(args.out, 'predict_misc', 'hints.E.gff')
hintsP = os.path.join(args.out, 'predict_misc', 'hints.P.gff')
hints_all = os.path.join(args.out, 'predict_misc', 'hints.PE.gff')

#check for masked genome here
if not os.path.isfile(MaskGenome) or lib.getSize(MaskGenome) < 10:
    lib.log.error("RepeatMasking failed, check log files.")
    sys.exit(1)

#load contig names and sizes into dictionary, get masked repeat stats
ContigSizes = {}
GenomeLength = 0
maskedSize = 0
with open(MaskGenome, 'rU') as input:
    for rec in SeqIO.parse(input, 'fasta'):
        if not rec.id in ContigSizes:
            ContigSizes[rec.id] = len(rec.seq)
            GenomeLength += len(rec.seq)
            maskedSize += lib.n_lower_chars(str(rec.seq))
        else:
            lib.log.error("Error, duplicate contig names, exiting")
            sys.exit(1)
percentMask = maskedSize / float(GenomeLength)
MaskedStats = '{0:.2f}%'.format(percentMask*10)
lib.log.info('Masked genome: {0:,}'.format(len(ContigSizes))+' scaffolds; {0:,}'.format(GenomeLength)+ ' bp; '+MaskedStats+' repeats masked')
            
#check for previous files and setup output files
Predictions = os.path.join(args.out, 'predict_misc', 'gene_predictions.gff3')
Exonerate = os.path.join(args.out, 'predict_misc', 'protein_alignments.gff3')
Transcripts = os.path.join(args.out, 'predict_misc', 'transcript_alignments.gff3')
Weights = os.path.join(args.out, 'predict_misc', 'weights.evm.txt')
EVM_out = os.path.join(args.out, 'predict_misc', 'evm.round1.gff3')
evminput = [Predictions, Exonerate, Transcripts]
for i in evminput:
    if os.path.isfile(i+'.old'):
        os.remove(i+'.old')
    if os.path.isfile(i):
        shutil.copyfile(i, i+'.old')

#if maker_gff passed, use that info and move on, if pasa present than run EVM.
if args.maker_gff:
    lib.log.info("Parsing Maker2 GFF for use in EVidence Modeler")
    maker2evm = os.path.join(parentdir, 'util', 'maker2evm.py')
    cmd = [sys.executable, maker2evm, os.path.abspath(args.maker_gff)]
    lib.runSubprocess(cmd, os.path.join(args.out, 'predict_misc'), lib.log)
    #append PASA data if exists
    if args.pasa_gff:
        with open(Predictions, 'a') as output:
            with open(PASA_GFF) as input:
                output.write(input.read())       
    #setup weights file for EVM
    with open(Weights, 'w') as output:
        genesources = []
        with open(Predictions, 'rU') as preds:
            for line in preds:
                if line.startswith('\n'):
                    continue
                source = line.split('\t')[1]
                if not source in genesources:
                    genesources.append(source)
        if not genesources:
            lib.log.error("Maker2 GFF not parsed correctly, no gene models found, exiting.")
            sys.exit(1)
        for i in genesources:
            if i == 'maker':
                output.write("ABINITIO_PREDICTION\t%s\t1\n" % i)
            else:
                output.write("OTHER_PREDICTION\t%s\t10\n" % i)  #set PASA to higher weight
        tr_sources = []
        with open(Transcripts, 'rU') as trns:
            for line in trns:
                source = line.split('\t')[1]
                if source not in tr_sources:
                    tr_sources.append(source)
        for i in tr_sources:
            output.write("TRANSCRIPT\t%s\t1\n" % i)
        output.write("PROTEIN\tprotein2genome\t1\n")

    Exonerate = os.path.abspath(Exonerate)
    Transcripts = os.path.abspath(Transcripts)
    
else:
    #no maker_gff, so let funannotate handle gene prediction
    #check for transcript evidence/format as needed
    trans_out = os.path.join(args.out, 'predict_misc', 'transcript_alignments.gff3')
    trans_temp = os.path.join(args.out, 'predict_misc', 'transcripts.combined.fa')
    blat_out = os.path.join(args.out, 'predict_misc', 'blat.psl')
    blat_filt = os.path.join(args.out, 'predict_misc', 'blat.filt.psl')
    blat_sort1 = os.path.join(args.out, 'predict_misc', 'blat.sort.tmp.psl')
    blat_sort2 = os.path.join(args.out, 'predict_misc', 'blat.sort.psl')
    maxINT = '-maxIntron='+str(args.max_intronlen)
    b2h_input = '--in='+blat_sort2
    b2h_output = '--out='+hintsE
    #combine transcript evidence into a single file
    if args.transcript_evidence:
        if os.path.isfile(trans_temp):
            shutil.copyfile(trans_temp, trans_temp+'.old')  
        with open(trans_temp, 'w') as output:
            for f in args.transcript_evidence:
                with open(f) as input:
                    output.write(input.read())
        else:
            shutil.copyfile(args.transcript_evidence, trans_temp)
        #check if old transcripts same as new ones, if different re-run GMAP/BLAT, otherwise use old if exists
        if os.path.isfile(trans_temp+'.old'):
            if not lib.sha256_check(trans_temp, trans_temp+'.old'): #they are not the same, re-run GMAP
                lib.log.info("Aligning transcript evidence to genome with GMAP")
                lib.runGMAP(trans_temp, MaskGenome, args.cpus, args.max_intronlen, os.path.join(args.out, 'predict_misc'), trans_out)
            else:
                os.remove(trans_temp+'.old')
                lib.log.info("Using existing transcript evidence alignments")
        #run Gmap of transcripts to genome
        if not os.path.isfile(trans_out):
            lib.log.info("Aligning transcript evidence to genome with GMAP")
            lib.runGMAP(trans_temp, MaskGenome, args.cpus, args.max_intronlen, os.path.join(args.out, 'predict_misc'), trans_out)
        Transcripts = os.path.abspath(trans_out)
    else:
        Transcripts = False
        if args.gmap_gff:
            shutil.copyfile(args.gmap_gff, trans_out)
            Transcripts = os.path.abspath(trans_out)
    if Transcripts:
        total = lib.countGMAPtranscripts(Transcripts)
        lib.log.info('{0:,}'.format(total) + ' transcripts aligned with GMAP')
    if not os.path.isfile(hintsE): #use previous hints file if exists
        if os.path.isfile(trans_temp): #if transcripts are available to algin, run BLAT
            #now run BLAT for Augustus hints
            lib.log.info("Aligning transcript evidence to genome with BLAT")
            cmd = ['blat', '-noHead', '-minIdentity=80', maxINT, MaskGenome, trans_temp, blat_out]
            lib.runSubprocess(cmd, '.', lib.log)
            cmd = ['pslCDnaFilter', '-minId=0.9', '-localNearBest=0.005', '-ignoreNs', '-bestOverlap', blat_out, blat_filt]
            lib.runSubprocess(cmd, '.', lib.log)
            cmd = ['sort', '-n', '-k', '16,16', blat_filt]
            lib.runSubprocess2(cmd, '.', lib.log, blat_sort1)
            cmd = ['sort', '-s', '-k', '14,14', blat_sort1]
            lib.runSubprocess2(cmd, '.', lib.log, blat_sort2)
            #run blat2hints
            blat2hints = os.path.join(AUGUSTUS_BASE, 'scripts', 'blat2hints.pl')
            cmd = [blat2hints, b2h_input, b2h_output, '--minintronlen=20', '--trunkSS']
            lib.runSubprocess(cmd, '.', lib.log)
            total = lib.line_count(blat_sort2)
            lib.log.info('{0:,}'.format(total) + ' filtered BLAT alignments')
        else:
            lib.log.error("No transcripts available to generate Augustus hints, provide --transcript_evidence")
            
    #check for protein evidence/format as needed
    p2g_out = os.path.join(args.out, 'predict_misc', 'exonerate.out')
    prot_temp = os.path.join(args.out, 'predict_misc', 'proteins.combined.fa')
    P2G = os.path.join(parentdir, 'bin', 'funannotate-p2g.py')
    if not args.exonerate_proteins:
        if args.protein_evidence:
            if os.path.isfile(prot_temp):
                shutil.copyfile(prot_temp, prot_temp+'.old')     
            #clean up headers, etc
            lib.cleanProteins(args.protein_evidence, prot_temp)
            #run funannotate-p2g to map to genome
            if args.use_diamond:
                p2g_cmd = [sys.executable, P2G, '-p', prot_temp, '-g', MaskGenome, '-o', p2g_out, '--maxintron', str(args.max_intronlen), '--cpus', str(args.cpus), '--ploidy', str(args.ploidy), '-f', 'diamond', '--tblastn_out', os.path.join(args.out, 'predict_misc', 'p2g.diamond.out'), '--logfile', os.path.join(args.out, 'logfiles', 'funannotate-p2g.log')]
            else:
                p2g_cmd = [sys.executable, P2G, '-p', prot_temp, '-g', MaskGenome, '-o', p2g_out, '--maxintron', str(args.max_intronlen), '--cpus', str(args.cpus), '--ploidy', str(args.ploidy), '-f', 'tblastn', '--tblastn_out', os.path.join(args.out, 'predict_misc', 'p2g.tblastn.out'), '--logfile', os.path.join(args.out, 'logfiles', 'funannotate-p2g.log')]               
            #check if protein evidence is same as old evidence
            if os.path.isfile(prot_temp+'.old'):
                if not lib.sha256_check(prot_temp, prot_temp+'.old'):
                    lib.log.info("Mapping proteins to genome using tBlastn/Exonerate")
                    subprocess.call(p2g_cmd)
                else:
                    lib.log.info("Using existing protein evidence alignments")
                    os.remove(prot_temp+'.old')
            if not os.path.isfile(p2g_out):
                lib.log.info("Mapping proteins to genome using tBlastn/Exonerate")
                subprocess.call(p2g_cmd)
            exonerate_out = os.path.abspath(p2g_out)
        else:
            exonerate_out = False
    else:
        shutil.copyfile(args.exonerate_proteins, p2g_out)
        exonerate_out = os.path.abspath(p2g_out)

    if exonerate_out:
        Exonerate = os.path.join(args.out, 'predict_misc', 'protein_alignments.gff3')
        with open(Exonerate, 'w') as output:
            try:
                subprocess.call([ExoConverter, exonerate_out], stdout = output, stderr = FNULL)
            except OSError:
                lib.log.error("$EVM_HOME variable is incorrect, please double-check: %s" % EVM)
                sys.exit(1)
        Exonerate = os.path.abspath(Exonerate)
        #now run exonerate2 hints for Augustus
        exonerate2hints = os.path.join(AUGUSTUS_BASE, 'scripts', 'exonerate2hints.pl')
        e2h_in = '--in='+p2g_out
        e2h_out = '--out='+hintsP
        e2h_minINT = '--minintronlen='+str(args.min_intronlen)
        e2h_maxINT = '--maxintronlen='+str(args.max_intronlen)
        cmd = [exonerate2hints, e2h_in, e2h_out, e2h_minINT, e2h_maxINT]
        lib.runSubprocess(cmd, '.', lib.log)

    #combine hints for Augustus
    if os.path.isfile(hintsP) or os.path.isfile(hintsE):
        if os.path.isfile(hints_all):
            os.remove(hints_all)
        with open(hints_all, 'a') as out:
            if os.path.isfile(hintsP):
                with open(hintsP) as input:
                    out.write(input.read())
            if os.path.isfile(hintsE):
                with open(hintsE) as input2:
                    out.write(input2.read())
    
    Augustus = ''
    GeneMark = ''

    #Walk thru data available and determine best approach. 
    if args.genemark_gtf:
        #convert the predictors to EVM format and merge
        #convert GeneMark
        GeneMarkGFF3 = os.path.join(args.out, 'predict_misc', 'genemark.gff')
        cmd = [GeneMark2GFF, args.genemark_gtf]
        lib.runSubprocess2(cmd, '.', lib.log, GeneMarkGFF3)
        GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
        cmd = ['perl', Converter, GeneMarkGFF3]
        lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp)
        GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
        with open(GeneMark, 'w') as output:
            with open(GeneMarkTemp, 'rU') as input:
                lines = input.read().replace("Augustus", "GeneMark")
                output.write(lines)

    if args.augustus_gff:
        #convert Augustus
        aug_out = args.augustus_gff
        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        cmd = ['perl', Converter, aug_out]
        lib.runSubprocess2(cmd, '.', lib.log, Augustus)

    if args.rna_bam and not any([GeneMark, Augustus]):
        if not args.augustus_species:
            aug_species = args.species.replace(' ', '_').lower()
        else:
            aug_species = args.augustus_species
        if lib.CheckAugustusSpecies(aug_species):
            lib.log.error("%s as already been trained, using existing parameters" % (aug_species))
        #now need to run BRAKER1
        braker_log = os.path.join(args.out, 'logfiles', 'braker.log')
        lib.log.info("Now launching BRAKER to train GeneMark and Augustus")
        species = '--species=' + aug_species
        genome = '--genome=' + MaskGenome
        bam = '--bam=' + os.path.abspath(args.rna_bam)
        Option1 = '--AUGUSTUS_CONFIG_PATH=' + AUGUSTUS
        Option2 = '--BAMTOOLS_PATH=' + BAMTOOLS_PATH
        Option3 = '--GENEMARK_PATH=' + GENEMARK_PATH
        #check if output is already there
        if not os.path.isdir(os.path.join(args.out, 'predict_misc', 'braker')):
            with open(braker_log, 'w') as logfile:
                if args.organism == 'fungus':
                    subprocess.call(['braker.pl', '--fungus', '--cores', str(args.cpus), Option1, Option2, Option3, '--gff3', '--softmasking', '1', genome, species, bam], stdout = logfile, stderr = logfile)
                else:
                    subprocess.call(['braker.pl', '--cores', str(args.cpus), Option1, Option2, Option3, '--gff3', '--softmasking', '1', genome, species, bam], stdout = logfile, stderr = logfile)
        #move braker output folder
        if os.path.isdir('braker'):
            if os.path.isdir(os.path.join(args.out, 'predict_misc', 'braker')):
                shutil.rmtree(os.path.join(args.out, 'predict_misc', 'braker'))
            os.rename('braker', os.path.join(args.out, 'predict_misc', 'braker'))
        #okay, now need to fetch the Augustus GFF and Genemark GTF files
        aug_out = os.path.join(args.out, 'predict_misc', 'braker', aug_species, 'augustus.gff')
        gene_out = os.path.join(args.out, 'predict_misc', 'braker', aug_species, 'GeneMark-ET', 'genemark.gtf')
        #now convert to EVM format
        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        cmd = ['perl', Converter2, aug_out]
        lib.runSubprocess2(cmd, '.', lib.log, Augustus)
        GeneMarkGFF3 = os.path.join(args.out, 'predict_misc', 'genemark.gff')
        cmd = [GeneMark2GFF, gene_out]
        lib.runSubprocess2(cmd, '.', lib.log, GeneMarkGFF3)
        GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
        cmd = ['perl', Converter, GeneMarkGFF3]
        lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp)
        GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
        with open(GeneMark, 'w') as output:
            with open(GeneMarkTemp, 'rU') as input:
                lines = input.read().replace("Augustus","GeneMark")
                output.write(lines)
    
    if args.pasa_gff and not Augustus:
        #use pasa models to train Augustus if not already trained
        if not args.augustus_species:
            aug_species = args.species.replace(' ', '_').lower()
        else:
            aug_species = args.augustus_species
        #setup final output
        aug_out = os.path.join(args.out, 'predict_misc', 'augustus.gff3')
        #check for training data, if no training data, then train using PASA
        if not lib.CheckAugustusSpecies(aug_species):
            lib.log.info("Training Augustus using PASA data, this may take awhile")
            GFF2GB = os.path.join(AUGUSTUS_BASE, 'scripts', 'gff2gbSmallDNA.pl')
            trainingset = os.path.join(args.out, 'predict_misc', 'augustus.pasa.gb')
            cmd = [GFF2GB, PASA_GFF, MaskGenome, '500', trainingset]
            lib.runSubprocess(cmd, '.', lib.log)
            if args.optimize_augustus:
                lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset, MaskGenome, args.out, args.cpus, True)   
            else:
                lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset, MaskGenome, args.out, args.cpus, False)
                
        #now run whole genome Augustus using trained parameters.
        lib.log.info("Running Augustus gene prediction")
        if not os.path.isfile(aug_out):     
            if os.path.isfile(hints_all):
                cmd = [AUGUSTUS_PARALELL, '--species', aug_species, '--hints', hints_all, '-i', MaskGenome, '-o', aug_out, '--cpus', str(args.cpus), '--logfile', os.path.join(args.out, 'logfiles', 'augustus-parallel.log')]
            else:
                cmd = [AUGUSTUS_PARALELL, '--species', aug_species, '-i', MaskGenome, '-o', aug_out, '--cpus', str(args.cpus), '--logfile', os.path.join(args.out, 'logfiles', 'augustus-parallel.log')]
            subprocess.call(cmd)
        #convert for EVM
        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        cmd = ['perl', Converter, aug_out]
        lib.runSubprocess2(cmd, '.', lib.log, Augustus)
   
    if not GeneMark:
        GeneMarkGFF3 = os.path.join(args.out, 'predict_misc', 'genemark.gff')
        #count contigs
        num_contigs = lib.countfasta(MaskGenome)
        #now run GeneMark-ES, first check for gmhmm mod file, use if available otherwise run ES
        if not args.genemark_mod:
            #if there are less than 2 data points (contigs, self-training fails), count contigs
            if num_contigs < 2:
                lib.log.error("GeneMark-ES cannot run with only a single contig, you must provide --ini_mod file to run GeneMark")
            else:
                if not os.path.isfile(GeneMarkGFF3):
                    if args.organism == 'fungus':
                        lib.RunGeneMarkES(MaskGenome, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, True)
                    else:
                        lib.RunGeneMarkES(MaskGenome, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, False)
                GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
                cmd = ['perl', Converter, GeneMarkGFF3]
                lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp)
                GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
                with open(GeneMark, 'w') as output:
                    with open(GeneMarkTemp, 'rU') as input:
                        lines = input.read().replace("Augustus", "GeneMark")
                        output.write(lines)
        else:   #have training parameters file, so just run genemark with
            if num_contigs < 2: #now can run modified prediction on single contig
                with open(MaskGenome, 'rU') as genome:
                    for line in genome:
                        if line.startswith('>'):
                            header = line.replace('>', '')
                            header = header.replace('\n', '')
                GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
                GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
                if not os.path.isfile(GeneMarkGFF3):
                    lib.log.info("Running GeneMark on single-contig assembly")
                    cmd = ['gmhmme3', '-m', args.genemark_mod, '-o', GeneMarkGFF3, '-f', 'gff3', MaskGenome]
                    lib.runSubprocess(cmd, '.', lib.log)
                #now open output and reformat
                lib.log.info("Converting GeneMark GTF file to GFF3")
                with open(GeneMarkTemp, 'w') as geneout:
                    with open(GeneMarkGFF3, 'rU') as genein:
                        for line in genein:
                            if not line.startswith('#'):
                                if not '\tIntron\t' in line:
                                    newline = line.replace('seq', header)
                                    newline = newline.replace('.hmm3', '')
                                    geneout.write(newline)
                GeneMarkTemp2 = os.path.join(args.out, 'predict_misc', 'genemark.temp2.gff')
                cmd = ['perl', Converter, GeneMarkTemp]
                lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp2)
                with open(GeneMark, 'w') as output:
                    with open(GeneMarkTemp2, 'rU') as input:
                        lines = input.read().replace("Augustus", "GeneMark")
                        output.write(lines)           
            else:
                if not os.path.isfile(GeneMarkGFF3):
                    if args.organism == 'fungus':
                        lib.RunGeneMark(MaskGenome, args.genemark_mod, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, True)
                    else:
                        lib.RunGeneMark(MaskGenome, args.genemark_mod, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, False)                  
                GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
                cmd = ['perl', Converter, GeneMarkGFF3]
                lib.runSubprocess2(cmd, '.', lib.log, GeneMarkTemp) 
                GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
                with open(GeneMark, 'w') as output:
                    with open(GeneMarkTemp, 'rU') as input:
                        lines = input.read().replace("Augustus", "GeneMark")
                        output.write(lines)

    if not Augustus: 
        if not args.augustus_species:
            aug_species = args.species.replace(' ', '_').lower()
        else:
            aug_species = args.augustus_species
        aug_out = os.path.join(args.out, 'predict_misc', 'augustus.gff3')
        busco_log = os.path.join(args.out, 'logfiles', 'busco.log')
        if not lib.CheckAugustusSpecies(aug_species):   
            #run BUSCO
            #define BUSCO and FUNGI models
            BUSCO = os.path.join(parentdir, 'util', 'funannotate-BUSCO2.py')
            BUSCO_FUNGI = os.path.join(FUNDB, args.busco_db)
            runbusco = True
            if os.path.isdir(os.path.join(args.out, 'predict_misc', 'busco')):
                #check if complete run
                if lib.checkannotations(os.path.join(args.out, 'predict_misc', 'busco', 'run_'+aug_species, 'training_set_'+aug_species)):
                    lib.log.info("BUSCO has already been run, using existing data")
                    runbusco = False
                else:
                    shutil.rmtree(os.path.join(args.out, 'predict_misc', 'busco'))                
            if runbusco:
                lib.log.info("Running BUSCO to find conserved gene models for training Augustus")
                if not os.path.isdir('busco'):
                    os.makedirs('busco')
                else:
                    shutil.rmtree('busco') #delete if it is there
                    os.makedirs('busco') #create fresh folder
                if lib.CheckAugustusSpecies(args.busco_seed_species):
                    busco_seed = args.busco_seed_species
                else:
                    busco_seed = 'generic'
                with open(busco_log, 'w') as logfile:
                    subprocess.call([sys.executable, BUSCO, '-i', MaskGenome, '-m', 'genome', '--lineage', BUSCO_FUNGI, '-o', aug_species, '-c', str(args.cpus), '--species', busco_seed, '-f'], cwd = 'busco', stdout = logfile, stderr = logfile)
                #check if BUSCO found models for training, if not error out and exit.
                busco_training = os.path.join('busco', 'run_'+aug_species, 'augustus_output', 'training_set_'+aug_species+'.txt')
                if not lib.checkannotations(busco_training):
                    lib.log.error("BUSCO training of Augusus failed, check busco logs, exiting")
                    sys.exit(1)
            #move the busco folder now where it should reside
            if os.path.isdir('busco'):
                if os.path.isdir(os.path.join(args.out, 'predict_misc', 'busco')):
                    shutil.rmtree(os.path.join(args.out, 'predict_misc', 'busco'))
                os.rename('busco', os.path.join(args.out, 'predict_misc', 'busco'))
            
            #open output and pull locations to make bed file
            busco_bed = os.path.join(args.out, 'predict_misc', 'buscos.bed')
            busco_fulltable = os.path.join(args.out, 'predict_misc', 'busco', 'run_'+aug_species, 'full_table_'+aug_species+'.tsv')
            busco_complete = lib.parseBUSCO2genome(busco_fulltable, args.ploidy, ContigSizes, busco_bed)
            
            #proper training files exist, now run EVM on busco models to get high quality predictions.
            lib.log.info('{0:,}'.format(len(busco_complete)) +' valid BUSCO predictions found, now formatting for EVM')

            #now get BUSCO GFF models
            busco_augustus_tmp = os.path.join(args.out, 'predict_misc', 'busco_augustus.tmp')
            with open(busco_augustus_tmp, 'w') as output:
                for i in busco_complete:
                    file = os.path.join(args.out, 'predict_misc', 'busco', 'run_'+aug_species, 'augustus_output', 'gffs', i+'.gff')
                    subprocess.call(['perl', Converter2, file], stderr = FNULL, stdout = output)
            #finally rename models so they are not redundant
            busco_augustus = os.path.join(args.out, 'predict_misc', 'busco_augustus.gff3')
            cmd = [os.path.join(parentdir, 'util', 'fix_busco_naming.py'), busco_augustus_tmp, busco_fulltable, busco_augustus]
            lib.runSubprocess(cmd, '.', lib.log)
            #now get genemark-es models in this region
            busco_genemark = os.path.join(args.out, 'predict_misc', 'busco_genemark.gff3')
            cmd = ['bedtools', 'intersect', '-a', GeneMark, '-b', busco_bed]
            lib.runSubprocess2(cmd, '.', lib.log, busco_genemark)
            #combine predictions
            busco_predictions = os.path.join(args.out, 'predict_misc', 'busco_predictions.gff3')
            with open(busco_predictions, 'w') as output:
                with open(busco_augustus) as input:
                    output.write(input.read())
                with open(busco_genemark) as input:
                    output.write(input.read())
            #get evidence if exists
            if Transcripts:
                #get transcript alignments in this region
                busco_transcripts = os.path.join(args.out, 'predict_misc', 'busco_transcripts.gff3')
                cmd = ['bedtools', 'intersect', '-a', Transcripts, '-b', busco_bed]
                lib.runSubprocess2(cmd, '.', lib.log, busco_transcripts)
            if Exonerate:
                #get protein alignments in this region
                busco_proteins = os.path.join(args.out, 'predict_misc', 'busco_proteins.gff3')
                cmd = ['bedtools', 'intersect', '-a', Exonerate, '-b', busco_bed]
                lib.runSubprocess2(cmd, '.', lib.log, busco_proteins)
            #set Weights file dependent on which data is present.
            busco_weights = os.path.join(args.out, 'predict_misc', 'busco_weights.txt')
            with open(busco_weights, 'w') as output:
                output.write("OTHER_PREDICTION\tAugustus\t2\n")
                output.write("ABINITIO_PREDICTION\tGeneMark\t1\n")
                if Exonerate:
                    output.write("PROTEIN\texonerate\t1\n")
                if Transcripts:
                    output.write("TRANSCRIPT\tgenome\t1\n")
            #setup EVM run
            EVM_busco = os.path.join(args.out, 'predict_misc', 'busco.evm.gff3')
            EVM_script = os.path.join(parentdir, 'bin', 'funannotate-runEVM.py')
            #get absolute paths for everything
            Busco_Weights = os.path.abspath(busco_weights)
            EVM_busco = os.path.abspath(EVM_busco)
            Busco_Predictions = os.path.abspath(busco_predictions)
            #parse entire EVM command to script, must be absolute paths for everything
            if Exonerate and Transcripts:
                evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM_busco.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Busco_Predictions, '--protein_alignments', os.path.abspath(busco_proteins), '--transcript_alignments', os.path.abspath(busco_transcripts), '--weights', Busco_Weights, '--min_intron_length', str(args.min_intronlen), EVM_busco]
            elif not Exonerate and Transcripts:
                evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM_busco.log'),str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Busco_Predictions, '--transcript_alignments', os.path.abspath(busco_transcripts), '--weights', Busco_Weights, '--min_intron_length', str(args.min_intronlen), EVM_busco]
            elif not Transcripts and Exonerate:
                evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM_busco.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Busco_Predictions, '--protein_alignments', os.path.abspath(busco_proteins), '--weights', Busco_Weights, '--min_intron_length', str(args.min_intronlen), EVM_busco]
            elif not any([Transcripts,Exonerate]):
                evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM_busco.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Busco_Predictions, '--weights', Busco_Weights, '--min_intron_length', str(args.min_intronlen), EVM_busco]
            #run EVM
            if not os.path.isfile(EVM_busco):
                subprocess.call(evm_cmd)
            try:
                total = lib.countGFFgenes(EVM_busco)
            except IOError:
                lib.log.error("EVM did not run correctly, output file missing")
                sys.exit(1)
            #check number of gene models, if 0 then failed, delete output file for re-running
            if total < 1:
                lib.log.error("Evidence modeler has failed, exiting")
                os.remove(EVM_busco)
                sys.exit(1)
            else:
                lib.log.info('{0:,}'.format(total) + ' total gene models from EVM')
            #move EVM folder to predict folder
            if os.path.isdir('EVM_tmp'):
                if os.path.isdir(os.path.join(args.out, 'predict_misc', 'EVM_busco')):
                    shutil.rmtree(os.path.join(args.out, 'predict_misc', 'EVM_busco'))
                os.rename('EVM_tmp', os.path.join(args.out, 'predict_misc', 'EVM_busco'))
            #convert to proteins and screen with busco
            lib.log.info("Checking BUSCO protein models for accuracy")
            evm_proteins = os.path.join(args.out, 'predict_misc', 'busco.evm.proteins.fa')
            busco_final = os.path.join(args.out, 'predict_misc', 'busco.final.gff3')
            cmd = [EVM2proteins, EVM_busco, MaskGenome]
            lib.runSubprocess2(cmd, '.', lib.log, evm_proteins)
            if not os.path.isdir(os.path.join(args.out, 'predict_misc', 'busco_proteins')):
                os.makedirs(os.path.join(args.out, 'predict_misc', 'busco_proteins'))
            with open(busco_log, 'a') as logfile:
                subprocess.call([sys.executable, BUSCO, '-i', os.path.abspath(evm_proteins), '-m', 'proteins', '--lineage', BUSCO_FUNGI, '-o', aug_species, '--cpu', str(args.cpus), '--species', busco_seed, '-f' ], cwd = os.path.join(args.out, 'predict_misc', 'busco_proteins'), stdout = logfile, stderr = logfile)
            subprocess.call([os.path.join(parentdir, 'util', 'filter_buscos.py'), EVM_busco, os.path.join(args.out, 'predict_misc', 'busco_proteins', 'run_'+aug_species, 'full_table_'+aug_species+'.tsv'), busco_final], stdout = FNULL, stderr = FNULL)
            total = lib.countGFFgenes(busco_final)
            lib.log.info('{0:,}'.format(total) + ' gene models validated, using for training Augustus')
        
            ###Run Augustus training
            GFF2GB = os.path.join(AUGUSTUS_BASE, 'scripts', 'gff2gbSmallDNA.pl')
            trainingset = os.path.join(args.out, 'predict_misc', 'busco.training.gb')
            cmd = [GFF2GB, busco_final, MaskGenome, '500', trainingset]
            lib.runSubprocess(cmd, '.', lib.log)
            if args.optimize_augustus:
                lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset, MaskGenome, args.out, args.cpus, True)
            else:
                lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset, MaskGenome, args.out, args.cpus, False)

        lib.log.info("Running Augustus gene prediction")
        #now run Augustus multithreaded...
        if not os.path.isfile(aug_out):
            if os.path.isfile(hints_all):
                cmd = [AUGUSTUS_PARALELL, '--species', aug_species, '--hints', hints_all, '-i', MaskGenome, '-o', aug_out, '--cpus', str(args.cpus), '--logfile', os.path.join(args.out, 'logfiles', 'augustus-parallel.log')]
            else:
                cmd = [AUGUSTUS_PARALELL, '--species', aug_species, '-i', MaskGenome, '-o', aug_out, '--cpus', str(args.cpus), '--logfile', os.path.join(args.out, 'logfiles', 'augustus-parallel.log')]
            subprocess.call(cmd)
        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        cmd = ['perl', Converter, aug_out]
        lib.runSubprocess2(cmd, '.', lib.log, Augustus)
        
    #just double-check that you've gotten here and both Augustus/GeneMark are finished
    if not any([Augustus, GeneMark]):
        lib.log.error("Augustus or GeneMark prediction is missing, check log files for errors")
        sys.exit(1)
    
    #GeneMark can fail if you try to pass a single contig, check file length
    GM_check = lib.line_count(GeneMark)
    gmc = 1
    if GM_check < 3:
        gmc = 0
        lib.log.error("GeneMark predictions failed, proceeding with only Augustus")
    
    #make sure Augustus finished successfully
    if not lib.checkannotations(Augustus):
        lib.log.error("Augustus prediction failed, check `logfiles/augustus-parallel.log`")
        sys.exit(1)
    
    #if hints used for Augustus, get high quality models > 90% coverage to pass to EVM
    if os.path.isfile(hints_all) or args.rna_bam:
        lib.log.info("Pulling out high quality Augustus predictions")
        hiQ_models = []
        with open(aug_out, 'rU') as augustus:
            for pred in lib.readBlocks(augustus, '# start gene'):
                values = []
                geneID = ''
                support = ''
                if pred[0].startswith('# This output'):
                    continue
                if pred[0].startswith('##gff-version 3'):
                    continue
                for line in pred:
                    line = line.replace('\n', '')
                    if line.startswith('# start gene'):
                        geneID = line.split(' ')[-1]
                        values.append(geneID)
                    if not args.rna_bam:
                        if line.startswith('# % of transcript supported by hints'):
                            support = line.split(' ')[-1]
                            values.append(support)
                    else: #if BRAKER is run then only intron CDS evidence is passed, so get models that fullfill that check
                        if line.startswith('# CDS introns:'):
                            intronMatch = line.split(' ')[-1]
                            try:
                                support = int(intronMatch.split('/')[0]) / float(intronMatch.split('/')[1]) * 100
                            except ZeroDivisionError:
                                support = 0
                            values.append(support)
                if float(values[1]) > 89: #greater than ~90% of exons supported, this is really stringent which is what we want here, as we are going to weight these models 5 to 1 over genemark
                    hiQ_models.append(values[0])

        #now open evm augustus and pull out models
        HiQ = set(hiQ_models)
        lib.log.info("Found %i high quality predictions from Augustus (>90%% exon evidence)"  % len(HiQ))
        HiQ_match = re.compile(r'\b(?:%s)[\.t1;]+\b' % '|'.join(HiQ))
        AugustusHiQ = os.path.join(args.out, 'predict_misc', 'augustus-HiQ.evm.gff3')
        with open(AugustusHiQ, 'w') as HiQ_out:
            with open(Augustus, 'rU') as evm_aug:
                for line in evm_aug:
                    if HiQ_match.search(line):
                        newline = line.replace('\tAugustus\t', '\tHiQ\t')
                        HiQ_out.write(newline)


    #EVM related input tasks, find all predictions and concatenate together
    pred_in = [Augustus, GeneMark]
    if args.pasa_gff:
        pred_in.append(PASA_GFF)
    if args.other_gff:
        pred_in.append(OTHER_GFF)
    if os.path.isfile(hints_all) or args.rna_bam:
        pred_in.append(AugustusHiQ)
    
    #write gene predictions file
    with open(Predictions+'.tmp', 'w') as output:
        for f in sorted(pred_in):
            with open(f) as input:
                output.write(input.read())
    #sort the predictions file
    lib.sortGFF(Predictions+'.tmp', Predictions, Scaffoldsort)

    #set Weights file dependent on which data is present.
    Weights = os.path.join(args.out, 'predict_misc', 'weights.evm.txt')
    with open(Weights, 'w') as output:
        output.write("ABINITIO_PREDICTION\tAugustus\t1\n")
        output.write("ABINITIO_PREDICTION\tGeneMark\t1\n")
        if os.path.isfile(hints_all) and not args.rna_bam:
            output.write("OTHER_PREDICTION\tHiQ\t5\n")
        if args.pasa_gff:
            output.write("OTHER_PREDICTION\tpasa_pred\t%s\n" % PASA_weight)
        if exonerate_out:
            output.write("PROTEIN\texonerate\t1\n")
        if Transcripts:
            output.write("TRANSCRIPT\tgenome\t1\n")
        if args.other_gff:
            output.write("OTHER_PREDICTION\tother_pred\t1\n" % OTHER_weight)

#total up Predictions
total = lib.countGFFgenes(Predictions)
lib.log.info('{0:,}'.format(total) + ' total gene models from all sources')
    
#setup EVM run
EVM_script = os.path.join(parentdir, 'bin', 'funannotate-runEVM.py')

#check if EVM input is identical as before
if os.path.isfile(Predictions+'.old'):
    if not lib.sha256_check(Predictions, Predictions+'.old'):
        #need to run EVM again, so delete output
        if os.path.isfile(EVM_out):
            os.remove(EVM_out)
    else:
        lib.log.info("Using existing EVM run data")

#get absolute paths for everything
Weights = os.path.abspath(Weights)
EVM_out = os.path.abspath(EVM_out)
Predictions = os.path.abspath(Predictions)

#parse entire EVM command to script
if Exonerate and Transcripts:
    Transcripts = os.path.abspath(Transcripts)
    Exonerate = os.path.abspath(Exonerate)
    evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--protein_alignments', Exonerate, '--transcript_alignments', Transcripts, '--weights', Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
elif not Exonerate and Transcripts:
    Transcripts = os.path.abspath(Transcripts)
    evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM.log'),str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--transcript_alignments', Transcripts, '--weights', Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
elif not Transcripts and Exonerate:
    Exonerate = os.path.abspath(Exonerate)
    evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--protein_alignments', Exonerate, '--weights', Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
elif not any([Transcripts,Exonerate]):
    evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--weights', Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]

#run EVM
if not os.path.isfile(EVM_out):
    subprocess.call(evm_cmd)
try:
    total = lib.countGFFgenes(EVM_out)
except IOError:
    lib.log.error("EVM did not run correctly, output file missing")
    sys.exit(1)
#check number of gene models, if 0 then failed, delete output file for re-running
if total < 1:
    lib.log.error("Evidence modeler has failed, exiting")
    os.remove(EVM_out)
    sys.exit(1)
else:
    lib.log.info('{0:,}'.format(total) + ' total gene models from EVM')

#move EVM folder to predict folder
if os.path.isdir('EVM_tmp'):
    if os.path.isdir(os.path.join(args.out, 'predict_misc', 'EVM')):
        shutil.rmtree(os.path.join(args.out, 'predict_misc', 'EVM'))
    os.rename('EVM_tmp', os.path.join(args.out, 'predict_misc', 'EVM'))

#run tRNAscan
lib.log.info("Predicting tRNAs")
tRNAscan = os.path.join(args.out, 'predict_misc', 'trnascan.gff3')
if not os.path.isfile(tRNAscan):
    lib.runtRNAscan(MaskGenome, os.path.join(args.out,'predict_misc'), tRNAscan)

#combine tRNAscan with EVM gff, dropping tRNA models if they overlap with EVM models
cleanTRNA = os.path.join(args.out, 'predict_misc', 'trnascan.no-overlaps.gff3')
cmd = ['bedtools', 'intersect', '-v', '-a', tRNAscan, '-b', EVM_out]
lib.runSubprocess2(cmd, '.', lib.log, cleanTRNA)
lib.log.debug("{:,} tRNAscan models are valid (do not overlap EVM models)".format(lib.countGFFgenes(cleanTRNA)))
lib.log.info("Merging EVM output with tRNAscan output")
gffs = [cleanTRNA, EVM_out]
GFF = os.path.join(args.out, 'predict_misc', 'evm.trnascan.gff')
if os.path.isfile(GFF):
    os.remove(GFF)
with open(GFF, 'w') as output:
    for f in gffs:
        with open(f, 'rU') as input:
            for line in input:
                if not line.startswith('\n'):
                    output.write(line)

#run GAG to get gff and proteins file for screening
lib.log.info("Reformatting GFF file using GAG")
gag1dir = os.path.join(args.out, 'predict_misc', 'gag1')
if os.path.isdir(gag1dir):
    shutil.rmtree(gag1dir)
cmd = ['gag.py', '-f', MaskGenome, '-g', GFF, '-o', gag1dir,'--fix_start_stop']
lib.runSubprocess(cmd, '.', lib.log)
GAG_gff = os.path.join(gag1dir, 'genome.gff')
GAG_proteins_original = os.path.join(gag1dir, 'genome.proteins.fasta')
GAG_proteins = os.path.join(args.out, 'predict_misc', 'gag1.proteins.fasta')
#clean up GAG proteins so names are consistent between versions
lib.GAGprotClean(GAG_proteins_original, GAG_proteins)
total = lib.countGFFgenes(GAG_gff)
lib.log.info('{0:,}'.format(total) + ' total gene models')

#filter bad models
lib.log.info("Filtering out bad gene models (< %i aa in length, transposable elements, etc)." % (args.min_protlen))
Blast_rep_remove = os.path.join(args.out, 'predict_misc', 'repeat.gene.models.txt')
if os.path.isfile(Blast_rep_remove): #need to run this every time if gene models have changed from a re-run
    os.remove(Blast_rep_remove)
lib.RepeatBlast(GAG_proteins, args.cpus, 1e-10, FUNDB, os.path.join(args.out, 'predict_misc'), Blast_rep_remove)
CleanGFF = os.path.join(args.out, 'predict_misc', 'cleaned.gff3')
if os.path.isfile(CleanGFF):
    os.remove(CleanGFF)
lib.RemoveBadModels(GAG_proteins, GAG_gff, args.min_protlen, RepeatMasker, Blast_rep_remove, os.path.join(args.out, 'predict_misc'), CleanGFF) 
total = lib.countGFFgenes(CleanGFF)
lib.log.info('{0:,}'.format(total) + ' gene models remaining')

#need to write to tbl2asn twice to fix errors, run first time and then parse error report
lib.log.info("Converting to preliminary Genbank format")
gag2dir = os.path.join(args.out, 'predict_misc', 'gag2')
if os.path.isdir(gag2dir):
    shutil.rmtree(gag2dir)
cmd = ['gag.py', '-f', MaskGenome, '-g', CleanGFF, '-o', gag2dir,'--fix_start_stop']
lib.runSubprocess(cmd, '.', lib.log)
shutil.copyfile(os.path.join(gag2dir, 'genome.fasta'), os.path.join(gag2dir, 'genome.fsa'))
SBT = os.path.join(parentdir, 'lib', 'test.sbt')
discrep = 'discrepency.report.txt'
lib.runtbl2asn(gag2dir, SBT, discrep, args.species, args.isolate, args.strain, args.tbl2asn)

#now parse error reports and remove bad models
lib.log.info("Cleaning models flagged by tbl2asn")
NCBIcleanGFF = os.path.join(args.out, 'predict_misc', 'ncbi.cleaned.gff3')
ErrSum = os.path.join(gag2dir, 'errorsummary.val')
Val = os.path.join(gag2dir, 'genome.val')
DirtyGFF = os.path.join(gag2dir, 'genome.gff')
if os.path.isfile(NCBIcleanGFF):
    os.remove(NCBIcleanGFF)
if args.keep_no_stops:
    lib.ParseErrorReport(DirtyGFF, ErrSum, Val, discrep, NCBIcleanGFF, keep_stops=True)
else:
    lib.ParseErrorReport(DirtyGFF, ErrSum, Val, discrep, NCBIcleanGFF, keep_stops=False)
total = lib.countGFFgenes(NCBIcleanGFF)
lib.log.info('{0:,}'.format(total) + ' gene models remaining')
shutil.copyfile(discrep, os.path.join(gag2dir, discrep))

#now we can rename gene models
lib.log.info("Re-naming gene models")
if os.path.isfile(os.path.join(args.out, 'predict_misc', 'ncbi.cleaned.gff3.bak')):
    os.remove(os.path.join(args.out, 'predict_misc', 'ncbi.cleaned.gff3.bak'))
shutil.copyfile(NCBIcleanGFF, os.path.join(args.out, 'predict_misc', 'ncbi.cleaned.gff3.bak'))
MAP = os.path.join(parentdir, 'util', 'maker_map_ids.pl')
MAPGFF = os.path.join(parentdir, 'util', 'map_gff_ids.pl')
mapping = os.path.join(args.out, 'predict_misc', 'mapping.ids')
if os.path.isfile(mapping):
    os.remove(mapping)
if not args.name.endswith('_'):
    args.name = args.name + '_'
cmd = ['perl', MAP, '--prefix', args.name, '--sort_order', Renamingsort, '--justify', '6', '--suffix', '-T', '--iterate', '1', NCBIcleanGFF]
lib.runSubprocess2(cmd, '.', lib.log, mapping)
cmd = ['perl', MAPGFF, mapping, NCBIcleanGFF]
lib.runSubprocess4(cmd, '.', lib.log)

#run GAG again with clean dataset, fix start/stops
gag3dir = os.path.join(args.out, 'predict_misc', 'tbl2asn')
if os.path.isdir(gag3dir):
    shutil.rmtree(gag3dir)
cmd = ['gag.py', '-f', MaskGenome, '-g', NCBIcleanGFF, '-o', gag3dir, '--fix_start_stop']
lib.runSubprocess(cmd, '.', lib.log)

#setup final output files
base = args.species.replace(' ', '_').lower()
final_fasta = os.path.join(args.out, 'predict_results', base + '.scaffolds.fa')
final_gff = os.path.join(args.out, 'predict_results', base + '.gff3')
final_gbk = os.path.join(args.out, 'predict_results', base + '.gbk')
final_tbl = os.path.join(args.out, 'predict_results', base + '.tbl')
final_proteins = os.path.join(args.out, 'predict_results', base + '.proteins.fa')
final_transcripts = os.path.join(args.out, 'predict_results', base + '.transcripts.fa')
final_validation = os.path.join(args.out, 'predict_results', base+'.validation.txt')
final_error = os.path.join(args.out, 'predict_results', base+'.error.summary.txt')

#run tbl2asn in new directory directory
shutil.copyfile(os.path.join(gag3dir, 'genome.fasta'), os.path.join(gag3dir, 'genome.fsa'))
discrep = os.path.join(args.out, 'predict_results', base + '.discrepency.report.txt')
lib.log.info("Converting to final Genbank format")
lib.runtbl2asn(gag3dir, SBT, discrep, args.species, args.isolate, args.strain, args.tbl2asn)

#retrieve files/reorganize
shutil.copyfile(os.path.join(gag3dir, 'genome.gff'), final_gff)
shutil.copyfile(os.path.join(gag3dir, 'genome.gbf'), final_gbk)
shutil.copyfile(os.path.join(gag3dir, 'genome.tbl'), final_tbl)
shutil.copyfile(os.path.join(gag3dir, 'genome.val'), final_validation)
shutil.copyfile(os.path.join(gag3dir, 'errorsummary.val'), final_error)
lib.log.info("Collecting final annotation files")
total = lib.countGFFgenes(final_gff)
lib.log.info('{0:,}'.format(total) + ' gene models')
lib.gb2output(final_gbk, final_proteins, final_transcripts, final_fasta)

lib.log.info("Funannotate predict is finished, output files are in the %s/predict_results folder" % (args.out))
lib.log.info("Note, you should fix any tbl2asn errors now before running functional annotation.")
if args.rna_bam and args.pasa_gff and os.path.isdir(args.out, 'training'): #give a suggested command
    lib.log.info("Your next step to capture UTRs and update annotation using PASA:\n\n\
funannotate update -i {:} --cpus {:}".format(args.out, args.cpus))
elif args.rna_bam: #means you have RNA-seq, but did not use funannotate train
    lib.log.info("Your next step to capture UTRs and update annotation using PASA:\n\n\
funannotate update -i {:} --cpus {:} \\\
        --left illumina_forward_RNAseq_R1.fastq.gz \\\
        --right illumina_forward_RNAseq_R2.fastq.gz \\\
        --jaccard clip\n".format(args.out, args.cpus))
else:
    lib.log.info("Your next step might be functional annotation, suggested commands:\
                \n\tRun EggNog-mapper: \n\t\temapper.py -i {:} -d fuNOG -o {:} --cpu {:}\
                \n\tRun InterProScan (Docker required): \n\t\t{:} -i={:} -c={:}\
                \n\tRun antiSMASH: \n\t\tfunannotate remote -i {:} -m antismash -e youremail@server.edu\
                \n\tAnnotate Genome: \n\t\tfunannotate annotate -i {:} --eggnog {:} \\\n\t\t--iprscan {:} --cpus {:} --sbt yourSBTfile.txt\
                ".format(os.path.join(args.out, 'predict_results', base+'.proteins.fa'), \
                base, \
                args.cpus, \
                os.path.join(parentdir, 'util', 'interproscan_docker.sh'), \
                os.path.join(args.out, 'predict_results', base+'.proteins.fa'), \
                args.cpus, \
                args.out, \
                args.out, \
                organism_name+'.emapper.annotations', \
                os.path.join(args.out, 'predict_results', base+'.proteins.fa.xml'), \
                args.cpus))
print("-------------------------------------------------------")

#clean up intermediate folders
if os.path.isfile('discrepency.report.txt'):
    os.rename('discrepency.report.txt', os.path.join(gag3dir, 'discrepency.report.txt'))
if os.path.isfile('funannotate-EVM.log'):
    os.rename('funannotate-EVM.log', os.path.join(args.out, 'logfiles', 'funannotate-EVM.log'))
if os.path.isfile('funannotate-p2g.log'):
    os.rename('funannotate-p2g.log', os.path.join(args.out, 'logfiles', 'funannotate-p2g.log'))
sys.exit(1)
