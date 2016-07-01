#!/usr/bin/env python

import sys, os, subprocess, inspect, multiprocessing, shutil, argparse, time, re, platform
from Bio import SeqIO
currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0,parentdir)
import lib.library as lib

#setup menu with argparse
class MyFormatter(argparse.ArgumentDefaultsHelpFormatter):
    def __init__(self,prog):
        super(MyFormatter,self).__init__(prog,max_help_position=48)
parser=argparse.ArgumentParser(prog='funannotate-predict.py', usage="%(prog)s [options] -i genome.fasta",
    description='''Script that does it all...''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('-s','--species', required=True, help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('--isolate', help='Isolate/strain name (e.g. Af293)')
parser.add_argument('--header_length', default=16, type=int, help='Max length for fasta headers')
parser.add_argument('--name', default="FUN_", help='Shortname for genes, perhaps assigned by NCBI, eg. VC83')
parser.add_argument('--augustus_species', help='Specify species for Augustus')
parser.add_argument('--genemark_mod', help='Use pre-existing Genemark training file (e.g. gmhmm.mod)')
parser.add_argument('--protein_evidence', default='uniprot.fa', help='Specify protein evidence (multiple files can be separaed by a comma)')
parser.add_argument('--exonerate_proteins', help='Pre-computed Exonerate protein alignments (see README for how to run exonerate)')
parser.add_argument('--transcript_evidence', help='Transcript evidence (map to genome with GMAP)')
parser.add_argument('--gmap_gff', help='Pre-computed GMAP transcript alignments (GFF3)')
parser.add_argument('--pasa_gff', help='Pre-computed PASA/TransDecoder high quality models')
parser.add_argument('--augustus_gff', help='Pre-computed Augustus gene models (GFF3)')
parser.add_argument('--genemark_gtf', help='Pre-computed GeneMark gene models (GTF)')
parser.add_argument('--maker_gff', help='MAKER2 GFF output')
parser.add_argument('--repeatmodeler_lib', help='Pre-computed RepeatModeler (or other) repetitive elements')
parser.add_argument('--rna_bam', help='BAM (sorted) of RNAseq aligned to reference for BRAKER1')
parser.add_argument('--min_intronlen', default=10, help='Minimum intron length for gene models')
parser.add_argument('--max_intronlen', default=3000, help='Maximum intron length for gene models')
parser.add_argument('--min_protlen', default=50, type=int, help='Minimum amino acid length for valid gene model')
parser.add_argument('--keep_no_stops', action='store_true', help='Keep gene models without valid stop codons')
parser.add_argument('--cpus', default=2, type=int, help='Number of CPUs to use')
parser.add_argument('--busco_seed_species', default='anidulans', help='Augustus species to use as initial training point for BUSCO')
parser.add_argument('--optimize_augustus', action='store_true', help='Run "long" training of Augustus')
parser.add_argument('--busco_db', default='fungi', choices=['fungi', 'metazoa', 'eukaryota', 'arthropoda', 'vertebrata'], help='BUSCO model database')
parser.add_argument('--organism', default='fungus', choices=['fungus', 'other'], help='Fungal specific settings')
parser.add_argument('--EVM_HOME', help='Path to Evidence Modeler home directory, $EVM_HOME')
parser.add_argument('--AUGUSTUS_CONFIG_PATH', help='Path to Augustus config directory, $AUGUSTUS_CONFIG_PATH')
parser.add_argument('--GENEMARK_PATH', help='Path to GeneMark exe (gmes_petap.pl) directory, $GENEMARK_PATH')
parser.add_argument('--BAMTOOLS_PATH', help='Path to BamTools exe directory, $BAMTOOLS_PATH')
args=parser.parse_args()

#check for conflicting folder names to avoid problems
conflict = ['busco', 'busco_proteins', 'RepeatMasker', 'RepeatModeler', 'genemark', 'EVM_tmp', 'braker']
if args.out in conflict:
    lib.log.error("%s output folder conflicts with a hard coded tmp folder, please change -o parameter" % args.out)
    os._exit(1)

#create folder structure
if not os.path.exists(args.out):
    os.makedirs(args.out)
    os.makedirs(os.path.join(args.out, 'predict_misc'))
    os.makedirs(os.path.join(args.out, 'predict_results'))
    os.makedirs(os.path.join(args.out, 'logfiles'))
    
#create log file
log_name = os.path.join(args.out, 'logfiles', 'funannotate-predict.log')
if os.path.isfile(log_name):
    os.remove(log_name)

#create debug log file (capture stderr)
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

#check for DB files needed for funanntoate predict, should only need REPEAT DB
blastdb = os.path.join(parentdir,'DB','REPEATS.psq')
if not os.path.isfile(blastdb):
    lib.log.error("funannotate database is not properly configured, please run `./setup.sh` in the %s directory" % parentdir)
    os._exit(1)
#check buscos, download if necessary
if not os.path.isdir(os.path.join(parentdir, 'DB', args.busco_db)):
    lib.download_buscos(args.busco_db)
    

#do some checks and balances
try:
    EVM = os.environ["EVM_HOME"]
except KeyError:
    if not args.EVM_HOME:
        lib.log.error("$EVM_HOME environmental variable not found, Evidence Modeler is not properly configured.  You can use the --EVM_HOME argument to specifiy a path at runtime")
        os._exit(1)
    else:
        EVM = args.EVM_HOME

try:
    AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
except KeyError:
    if not args.AUGUSTUS_CONFIG_PATH:
        lib.log.error("$AUGUSTUS_CONFIG_PATH environmental variable not found, Augustus is not properly configured. You can use the --AUGUSTUS_CONFIG_PATH argument to specify a path at runtime.")
        os._exit(1)
    else:
        AUGUSTUS = args.AUGUSTUS_CONFIG_PATH
        
#if you want to use BRAKER1, you also need some additional config paths
try:
    GENEMARK_PATH = os.environ["GENEMARK_PATH"]
except KeyError:
    if not lib.which('gmes_petap.pl'):
        if not args.GENEMARK_PATH:
            lib.log.error("GeneMark not found and $GENEMARK_PATH environmental variable missing, BRAKER1 is not properly configured. You can use the --GENEMARK_PATH argument to specify a path at runtime.")
            os._exit(1)
        else:
            GENEMARK_PATH = args.GENEMARK_PATH

try:
    BAMTOOLS_PATH = os.environ["BAMTOOLS_PATH"]
except KeyError:
    #check if it is in PATH, if it is, no problem, else through warning
    if not lib.which('bamtools'):
        if not args.BAMTOOLS_PATH:
            lib.log.error("Bamtools not found and $BAMTOOLS_PATH environmental variable missing, BRAKER1 is not properly configured. You can use the --BAMTOOLS_PATH argument to specify a path at runtime.")
            os._exit(1)
        else:
            BAMTOOLS_PATH = args.BAMTOOLS_PATH

if AUGUSTUS.endswith('config'):
    AUGUSTUS_BASE = AUGUSTUS.replace('config', '')
elif AUGUSTUS.endswith('config'+os.sep):
    AUGUSTUS_BASE = AUGUSTUS.replace('config'+os.sep, '')
AutoAug = os.path.join(AUGUSTUS_BASE, 'scripts', 'autoAug.pl')
GeneMark2GFF = os.path.join(parentdir, 'util', 'genemark_gtf2gff3.pl')

programs = ['tblastn', 'exonerate', 'makeblastdb','dustmasker','gag.py','tbl2asn','gmes_petap.pl', 'BuildDatabase', 'RepeatModeler', 'RepeatMasker', GeneMark2GFF, AutoAug, 'bedtools', 'gmap', 'gmap_build', 'blat', 'pslCDnaFilter', 'augustus', 'etraining', 'rmOutToGFF3.pl']
lib.CheckDependencies(programs)

#check augustus species now, so that you don't get through script and then find out it is already in DB
if not args.augustus_species:
    aug_species = args.species.replace(' ', '_').lower()
else:
    aug_species = args.augustus_species
augspeciescheck = lib.CheckAugustusSpecies(aug_species)
if augspeciescheck:
    lib.log.error("Augustus training set for %s already exists, thus funannotate will use those parameters. If you want to re-train, provide a unique name for the --augustus_species argument" % (aug_species))

#check augustus functionality
augustuscheck = lib.checkAugustusFunc(AUGUSTUS_BASE)
if args.rna_bam:
    if augustuscheck[1] == 0:
        lib.log.error("ERROR: %s is not installed properly for BRAKER1 (check bam2hints compilation)" % augustuscheck[0])
        os._exit(1)
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
            os._exit(1)
        else:
            lib.log.info("Will proceed with PASA models to train Augustus")
            
#if made it here output Augustus version
lib.log.info("%s detected, version seems to be compatible with BRAKER1 and BUSCO" % augustuscheck[0])

if args.protein_evidence == 'uniprot.fa':
    args.protein_evidence = os.path.join(parentdir, 'DB', 'uniprot_sprot.fasta')
    
#check input files to make sure they are not empty, first check if multiple files passed to transcript/protein evidence
input_checks = [args.input, args.genemark_mod, args.exonerate_proteins, args.gmap_gff, args.pasa_gff, args.repeatmodeler_lib, args.rna_bam]
if ',' in args.protein_evidence: #there will always be something here, since defaults to uniprot
    prot_evid = args.protein_evidence.split(',')
    for x in prot_evid:
        input_checks.append(x)
else:
    input_checks.append(args.protein_evidence)
if args.transcript_evidence:  #if transcripts passed, otherwise ignore
    if ',' in args.transcript_evidence:
        trans_evid = args.transcript_evidence.split(',')
        for y in trans_evid:
            input_checks.append(y)
    else:
        input_checks.append(args.transcript_evidence)
#now check the inputs
for i in input_checks:
    if i:
        lib.checkinputs(i)

#check fasta header length
header_test = lib.checkFastaHeaders(args.input, args.header_length)
if not header_test:
    lib.log.error("Fasta headers on your input have more characters than the max (16), reformat headers to continue.")
    os._exit(1)
    
#setup augustus parallel command
AUGUSTUS_PARALELL = os.path.join(parentdir, 'bin', 'augustus_parallel.py')

#EVM command line scripts
Converter = os.path.join(EVM, 'EvmUtils', 'misc', 'augustus_GFF3_to_EVM_GFF3.pl')
ExoConverter = os.path.join(EVM, 'EvmUtils', 'misc', 'exonerate_gff_to_alignment_gff3.pl')
Validator = os.path.join(EVM, 'EvmUtils', 'gff3_gene_prediction_file_validator.pl')
Converter2 = os.path.join(EVM, 'EvmUtils', 'misc', 'augustus_GTF_to_EVM_GFF3.pl')
EVM2proteins = os.path.join(EVM, 'EvmUtils', 'gff3_file_to_proteins.pl')

#so first thing is to reformat genome fasta only if there is no aligned evidence already
sort_out = os.path.join(args.out, 'predict_misc', 'genome.fasta')
#just copy the input fasta to the misc folder and move on.
shutil.copyfile(args.input, sort_out)
Genome = os.path.abspath(sort_out)

#repeatmasker, run if not passed from command line
if not args.repeatmodeler_lib:
    MaskGenome = os.path.join(args.out, 'predict_misc', 'genome.softmasked.fa')
    if not os.path.isfile(MaskGenome):
        lib.RepeatModelMask(Genome, args.cpus, os.path.join(args.out, 'predict_misc'), MaskGenome, debug)
else:
    MaskGenome = os.path.join(args.out, 'predict_misc', 'genome.softmasked.fa')
    if not os.path.isfile(MaskGenome):
        lib.RepeatMask(Genome, args.repeatmodeler_lib, args.cpus, os.path.join(args.out, 'predict_misc'), MaskGenome, debug)
RepeatMasker = os.path.join(args.out, 'predict_misc', 'repeatmasker.gff3')
RepeatMasker = os.path.abspath(RepeatMasker)
MaskGenome = os.path.abspath(MaskGenome)
#clean up folders
if os.path.isdir('RepeatModeler'):
    os.rename('RepeatModeler', os.path.join(args.out, 'predict_misc', 'RepeatModeler'))
if os.path.isdir('RepeatMasker'):
    os.rename('RepeatMasker', os.path.join(args.out, 'predict_misc', 'RepeatMasker'))

#final output for augustus hints, declare ahead of time for checking portion of script
hintsE = os.path.join(args.out, 'predict_misc', 'hints.E.gff')
hintsP = os.path.join(args.out, 'predict_misc', 'hints.P.gff')
hints_all = os.path.join(args.out, 'predict_misc', 'hints.PE.gff')

#check for masked genome here
if not os.path.isfile(MaskGenome) or lib.getSize(MaskGenome) < 10:
    lib.log.error("RepeatMasking failed, check log files.")
    os._exit(1)

#if maker_gff passed, use that info and move directly to EVM
if args.maker_gff:
    lib.log.info("Maker2 GFF passed, parsing results and proceeding directly to EVidence Modeler")
    maker2evm = os.path.join(parentdir, 'util', 'maker2evm.pl')
    subprocess.call(['perl', maker2evm, os.path.abspath(args.maker_gff)], cwd = os.path.join(args.out, 'predict_misc'))
    Predictions = os.path.join(args.out, 'predict_misc', 'gene_predictions.gff3')
    Exonerate = os.path.join(args.out, 'predict_misc', 'protein_alignments.gff3')
    Transcripts = os.path.join(args.out, 'predict_misc', 'transcript_alignments.gff3')
    
    #append PASA data if passed
    if args.pasa_gff:
        with open(Predictions, 'a') as output:
            with open(args.pasa_gff) as input:
                output.write(input.read())       
    #setup weights file for EVM
    Weights = os.path.join(args.out, 'predict_misc', 'weights.evm.txt')
    with open(Weights, 'w') as output:
        sources = []
        with open(Predictions, 'rU') as preds:
            for line in preds:
                source = line.split('\t')[1]
                if source not in sources:
                    sources.append(source)
        if args.pasa_gff:
            output.write("OTHER_PREDICTION\ttransdecoder\t10\n")
        for i in sources:
            output.write("ABINITIO_PREDICTION\t%s\t1\n" % i)
        output.write("PROTEIN\tprotein2genome\t1\n")
        output.write("TRANSCRIPT\test2genome\t1\n")
    Exonerate = os.path.abspath(Exonerate)
    Transcripts = os.path.abspath(Transcripts)
    
else:
    #no maker_gff, so let funannotate handle gene prediction
    #check for transcript evidence/format as needed
    trans_out = os.path.join(args.out, 'predict_misc', 'transcript_alignments.gff3')
    if not args.gmap_gff:
        if args.transcript_evidence:
            trans_temp = os.path.join(args.out, 'predict_misc', 'transcripts.combined.fa')
            if ',' in args.transcript_evidence:
                files = args.transcript_evidence.split(",")
                with open(trans_temp, 'w') as output:
                    for f in files:
                        with open(f) as input:
                            output.write(input.read())
            else:
                shutil.copyfile(args.transcript_evidence, trans_temp)     
            #run Gmap of transcripts to genome
            lib.log.info("Aligning transcript evidence to genome with GMAP")
            if not os.path.isfile(trans_out):
                lib.runGMAP(trans_temp, MaskGenome, args.cpus, args.max_intronlen, os.path.join(args.out, 'predict_misc'), trans_out)
            Transcripts = os.path.abspath(trans_out)
            #now run BLAT for Augustus hints
            blat_out = os.path.join(args.out, 'predict_misc', 'blat.psl')
            blat_filt = os.path.join(args.out, 'predict_misc', 'blat.filt.psl')
            blat_sort1 = os.path.join(args.out, 'predict_misc', 'blat.sort.tmp.psl')
            blat_sort2 = os.path.join(args.out, 'predict_misc', 'blat.sort.psl')
            maxINT = '-maxIntron='+str(args.max_intronlen)
            lib.log.info("Aligning transcript evidence to genome with BLAT")
            if not os.path.isfile(hints_all):
                subprocess.call(['blat', '-noHead', '-minIdentity=80', maxINT, MaskGenome, trans_temp, blat_out], stdout=FNULL, stderr=FNULL)
                subprocess.call(['pslCDnaFilter', '-minId=0.9', '-localNearBest=0.005', '-ignoreNs', '-bestOverlap', blat_out, blat_filt], stdout=FNULL, stderr=FNULL)
                with open(blat_sort1, 'w') as output:
                    subprocess.call(['sort', '-n', '-k', '16,16', blat_filt], stdout=output, stderr=FNULL)
                with open(blat_sort2, 'w') as output:
                    subprocess.call(['sort', '-s', '-k', '14,14', blat_sort1], stdout=output, stderr=FNULL)
                #run blat2hints
                blat2hints = os.path.join(AUGUSTUS_BASE, 'scripts', 'blat2hints.pl')
                b2h_input = '--in='+blat_sort2
                b2h_output = '--out='+hintsE
                subprocess.call([blat2hints, b2h_input, b2h_output, '--minintronlen=20', '--trunkSS'], stdout=FNULL, stderr=FNULL)
        else:
            Transcripts = False
    else:
        shutil.copyfile(args.gmap_gff, trans_out)
        Transcripts = os.path.abspath(trans_out)

    #check for protein evidence/format as needed
    p2g_out = os.path.join(args.out, 'predict_misc', 'exonerate.out')
    if not args.exonerate_proteins:
        if args.protein_evidence:
            prot_temp = os.path.join(args.out, 'predict_misc', 'proteins.combined.fa')
            if ',' in args.protein_evidence:
                files = args.protein_evidence.split(",")
                with open(prot_temp, 'w') as output:
                    for f in files:
                        with open(f) as input:
                            output.write(input.read())
            else:
                shutil.copyfile(args.protein_evidence, prot_temp)
            #run funannotate-p2g to map to genome
            lib.log.info("Mapping proteins to genome using tBlastn/Exonerate")
            P2G = os.path.join(parentdir, 'bin','funannotate-p2g.py')
            p2g_cmd = [sys.executable, P2G, prot_temp, MaskGenome, p2g_out, str(args.max_intronlen), str(args.cpus), os.path.join(args.out, 'logfiles', 'funannotate-p2g.log')]
            if not os.path.isfile(p2g_out):
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
            subprocess.call([ExoConverter, exonerate_out], stdout = output, stderr = FNULL)
        Exonerate = os.path.abspath(Exonerate)
        #now run exonerate2 hints for Augustus
        exonerate2hints = os.path.join(AUGUSTUS_BASE, 'scripts', 'exonerate2hints.pl')
        e2h_in = '--in='+p2g_out
        e2h_out = '--out='+hintsP
        e2h_minINT = '--minintronlen='+str(args.min_intronlen)
        e2h_maxINT = '--maxintronlen='+str(args.max_intronlen)
        subprocess.call([exonerate2hints, e2h_in, e2h_out, e2h_minINT, e2h_maxINT], stdout=FNULL, stderr=FNULL)

    #combine hints for Augustus
    if os.path.isfile(hintsP) or os.path.isfile(hintsE):
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
        with open(GeneMarkGFF3, 'w') as output:
            subprocess.call([GeneMark2GFF, args.genemark_gtf], stdout = output, stderr = FNULL)
        GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
        with open(GeneMarkTemp, 'w') as output:
            subprocess.call(['perl', Converter, GeneMarkGFF3], stdout = output, stderr = FNULL)
        GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
        with open(GeneMark, 'w') as output:
            with open(GeneMarkTemp, 'rU') as input:
                lines = input.read().replace("Augustus","GeneMark")
                output.write(lines)

    if args.augustus_gff:
        #convert Augustus
        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        with open(Augustus, 'w') as output:
            subprocess.call(['perl', Converter, args.augustus_gff], stdout = output, stderr = FNULL)

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
        with open(braker_log, 'w') as logfile:
            if args.organism == 'fungus':
                subprocess.call(['braker.pl', '--fungus', '--cores', str(args.cpus), Option1, Option2, Option3, '--gff3', '--softmasking', '1', genome, species, bam], stdout = logfile, stderr = logfile)
            else:
                subprocess.call(['braker.pl', '--cores', str(args.cpus), Option1, Option2, Option3, '--gff3', '--softmasking', '1', genome, species, bam], stdout = logfile, stderr = logfile)
        #okay, now need to fetch the Augustus GFF and Genemark GTF files
        aug_out = os.path.join('braker', aug_species, 'augustus.gff3')
        gene_out = os.path.join('braker', aug_species, 'GeneMark-ET', 'genemark.gtf')
        #now convert to EVM format
        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        with open(Augustus, 'w') as output:
            subprocess.call(['perl', Converter, aug_out], stdout = output, stderr = FNULL)
        GeneMarkGFF3 = os.path.join(args.out, 'predict_misc', 'genemark.gff')
        with open(GeneMarkGFF3, 'w') as output:
            subprocess.call([GeneMark2GFF, gene_out], stdout = output, stderr = FNULL)
        GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
        with open(GeneMarkTemp, 'w') as output:
            subprocess.call(['perl', Converter, GeneMarkGFF3], stdout = output, stderr = FNULL)
        GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
        with open(GeneMark, 'w') as output:
            with open(GeneMarkTemp, 'rU') as input:
                lines = input.read().replace("Augustus","GeneMark")
                output.write(lines
        if os.path.isdir('braker'):
            if os.path.isdir(os.path.join(args.out, 'predict_misc', 'braker')):
                shutil.rmtree(os.path.join(args.out, 'predict_misc', 'braker'))
            os.rename('braker', os.path.join(args.out, 'predict_misc', 'braker'))

    
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
            subprocess.call([GFF2GB, args.pasa_gff, MaskGenome, '250', trainingset], stderr = FNULL)
            if args.optimize_augustus:
                lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset, MaskGenome, args.out, args.cpus, True)   
            else:
                lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset, MaskGenome, args.out, args.cpus, False)
                
        #now run whole genome Augustus using trained parameters.
        lib.log.info("Running Augustus gene prediction")
        if not os.path.isfile(aug_out):     
            if os.path.isfile(hints_all):
                subprocess.call([AUGUSTUS_PARALELL, '--species', aug_species, '--hints', hints_all, '-i', MaskGenome, '-o', aug_out, '--cpus', str(args.cpus)], stderr = FNULL, stdout = FNULL)
            else:
                subprocess.call([AUGUSTUS_PARALELL, '--species', aug_species, '-i', MaskGenome, '-o', aug_out, '--cpus', str(args.cpus)], stderr = FNULL, stdout = FNULL)
        #convert for EVM
        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        with open(Augustus, 'w') as output:
            subprocess.call(['perl', Converter, aug_out], stdout = output, stderr = FNULL)
   
    if not GeneMark:
        #count contigs
        num_contigs = lib.countfasta(MaskGenome)
        #now run GeneMark-ES, first check for gmhmm mod file, use if available otherwise run ES
        if not args.genemark_mod:
            #if there are less than 2 data points (contigs, self-training fails), count contigs
            if num_contigs < 2:
                lib.log.error("GeneMark-ES cannot run with only a single contig, you must provide --ini_mod file to run GeneMark")
            else:
                GeneMarkGFF3 = os.path.join(args.out, 'predict_misc', 'genemark.gff')
                if not os.path.isfile(GeneMarkGFF3):
                    if args.organism == 'fungus':
                        lib.RunGeneMarkES(MaskGenome, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, True)
                    else:
                        lib.RunGeneMarkES(MaskGenome, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, False)
                GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
                with open(GeneMarkTemp, 'w') as output:
                    subprocess.call(['perl', Converter, GeneMarkGFF3], stdout = output, stderr = FNULL)
                GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
                with open(GeneMark, 'w') as output:
                    with open(GeneMarkTemp, 'rU') as input:
                        lines = input.read().replace("Augustus","GeneMark")
                        output.write(lines)
        else:   #have training parameters file, so just run genemark with
            GeneMarkGFF3 = os.path.join(args.out, 'predict_misc', 'genemark.gff')
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
                    subprocess.call(['gmhmme3', '-m', args.genemark_mod, '-o', GeneMarkGFF3, '-f', 'gff3', MaskGenome])
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
                with open(GeneMarkTemp2, 'w') as output:
                    subprocess.call(['perl', Converter, GeneMarkTemp], stdout = output, stderr = FNULL)                    
                with open(GeneMark, 'w') as output:
                    with open(GeneMarkTemp2, 'rU') as input:
                        lines = input.read().replace("Augustus","GeneMark")
                        output.write(lines)
             
            else:
                GeneMarkGFF3 = os.path.join(args.out, 'predict_misc', 'genemark.gff')
                if not os.path.isfile(GeneMarkGFF3):
                    if args.organism == 'fungus':
                        lib.RunGeneMark(MaskGenome, args.genemark_mod, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, True)
                    else:
                        lib.RunGeneMark(MaskGenome, args.genemark_mod, args.cpus, os.path.join(args.out, 'predict_misc'), GeneMarkGFF3, False)                  
                GeneMarkTemp = os.path.join(args.out, 'predict_misc', 'genemark.temp.gff')
                with open(GeneMarkTemp, 'w') as output:
                    subprocess.call(['perl', Converter, GeneMarkGFF3], stdout = output, stderr = FNULL)
                GeneMark = os.path.join(args.out, 'predict_misc', 'genemark.evm.gff3')
                with open(GeneMark, 'w') as output:
                    with open(GeneMarkTemp, 'rU') as input:
                        lines = input.read().replace("Augustus","GeneMark")
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
            BUSCO = os.path.join(parentdir, 'util', 'funannotate-BUSCO.py')
            BUSCO_FUNGI = os.path.join(parentdir, 'DB', args.busco_db)
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
                    subprocess.call([sys.executable, BUSCO, '--genome', MaskGenome, '--lineage', BUSCO_FUNGI, '-o', aug_species, '--cpu', str(args.cpus), '--species', busco_seed], cwd = 'busco', stdout = logfile, stderr = logfile)
                #check if BUSCO found models for training, if not error out and exit.
                busco_training = os.path.join('busco', 'run_'+aug_species, 'training_set_'+aug_species)
                if not lib.checkannotations(busco_training):
                    lib.log.error("BUSCO training of Augusus failed, check busco logs, exiting")
                    #remove the augustus training config folder
                    shutil.rmtree(os.path.join(AUGUSTUS, 'species', aug_species))
                    os._exit(1)
            #proper training files exist, now run EVM on busco models to get high quality predictions.
            lib.log.info("BUSCO predictions complete, now formatting for EVM")
            #move the busco folder now where it should reside
            if os.path.isdir('busco'):
                if os.path.isdir(os.path.join(args.out, 'predict_misc', 'busco')):
                    shutil.rmtree(os.path.join(args.out, 'predict_misc', 'busco'))
                os.rename('busco', os.path.join(args.out, 'predict_misc', 'busco'))
            
            #open output and pull locations to make bed file
            busco_bed = os.path.join(args.out, 'predict_misc', 'buscos.bed')
            busco_fulltable = os.path.join(args.out, 'predict_misc', 'busco', 'run_'+aug_species, 'full_table_'+aug_species)
            busco_complete = {}
            with open(busco_bed, 'w') as bedfile:
                with open(busco_fulltable, 'rU') as buscoinput:
                    for line in buscoinput:
                        line = line.replace(']', '')
                        cols = line.split('\t')
                        if cols[1] == 'Complete':
                            if not cols[0] in busco_complete:
                                busco_complete[cols[0]] = cols[2]+':'+cols[3]+'-'+cols[4]
                            start = int(cols[3]) - 100
                            end = int(cols[4]) + 100
                            bedfile.write('%s\t%i\t%i\t%s\n' % (cols[2],start,end,cols[0]))
            #now get BUSCO GFF models
            busco_augustus_tmp = os.path.join(args.out, 'predict_misc', 'busco_augustus.tmp')
            with open(busco_augustus_tmp, 'w') as output:
                for i in busco_complete:
                    file = os.path.join(args.out, 'predict_misc', 'busco', 'run_'+aug_species, 'gffs', i+'.gff')
                    subprocess.call(['perl', Converter2, file], stdout = output)
            #finally rename models so they are not redundant
            busco_augustus = os.path.join(args.out, 'predict_misc', 'busco_augustus.gff3')
            subprocess.call([os.path.join(parentdir, 'util', 'fix_busco_naming.py'), busco_augustus_tmp, busco_fulltable, busco_augustus])
            #now get genemark-es models in this region
            busco_genemark = os.path.join(args.out, 'predict_misc', 'busco_genemark.gff3')
            with open(busco_genemark, 'w') as output:
                subprocess.call(['bedtools', 'intersect', '-a', GeneMark, '-b', busco_bed], stdout = output)
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
                with open(busco_transcripts, 'w') as output:
                    subprocess.call(['bedtools', 'intersect', '-a', Transcripts, '-b', busco_bed], stdout = output, stderr = FNULL)
            if Exonerate:
                #get protein alignments in this region
                busco_proteins = os.path.join(args.out, 'predict_misc', 'busco_proteins.gff3')
                with open(busco_proteins, 'w') as output:
                    subprocess.call(['bedtools', 'intersect', '-a', Exonerate, '-b', busco_bed], stdout = output, stderr = FNULL)           
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
            EVM_out = os.path.join(args.out, 'predict_misc', 'busco.evm.gff3')
            EVM_script = os.path.join(parentdir, 'bin', 'funannotate-runEVM.py')
            #get absolute paths for everything
            Busco_Weights = os.path.abspath(busco_weights)
            EVM_out = os.path.abspath(EVM_out)
            Busco_Predictions = os.path.abspath(busco_predictions)
            #parse entire EVM command to script
            if Exonerate and Transcripts:
                evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM_busco.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Busco_Predictions, '--protein_alignments', Exonerate, '--transcript_alignments', Transcripts, '--weights', Busco_Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
            elif not Exonerate and Transcripts:
                evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM_busco.log'),str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Busco_Predictions, '--transcript_alignments', Transcripts, '--weights', Busco_Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
            elif not Transcripts and Exonerate:
                evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM_busco.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Busco_Predictions, '--protein_alignments', Exonerate, '--weights', Busco_Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
            elif not any([Transcripts,Exonerate]):
                evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM_busco.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Busco_Predictions, '--weights', Busco_Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
            #run EVM
            if not os.path.isfile(EVM_out):
                subprocess.call(evm_cmd)
            try:
                total = lib.countGFFgenes(EVM_out)
            except IOError:
                lib.log.error("EVM did not run correctly, output file missing")
                os._exit(1)
            #check number of gene models, if 0 then failed, delete output file for re-running
            if total < 1:
                lib.log.error("Evidence modeler has failed, exiting")
                os.remove(EVM_out)
                os._exit(1)
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
            with open(evm_proteins, 'w') as output:
                subprocess.call([EVM2proteins, EVM_out, MaskGenome], stdout = output)
            if not os.path.isdir(os.path.join(args.out, 'predict_misc', 'busco_proteins')):
                os.makedirs(os.path.join(args.out, 'predict_misc', 'busco_proteins'))
            with open(busco_log, 'a') as logfile:
                subprocess.call([sys.executable, BUSCO, '-in', os.path.abspath(evm_proteins), '--lineage', BUSCO_FUNGI, '-a', aug_species, '--cpu', str(args.cpus), '-m', 'OGS', '-f'], cwd = os.path.join(args.out, 'predict_misc', 'busco_proteins'), stdout = logfile, stderr = logfile)
            subprocess.call([os.path.join(parentdir, 'util', 'filter_buscos.py'), EVM_out, os.path.join(args.out, 'predict_misc', 'busco_proteins', 'run_'+aug_species, 'full_table_'+aug_species), busco_final], stdout = FNULL, stderr = FNULL)
            total = lib.countGFFgenes(busco_final)
            lib.log.info('{0:,}'.format(total) + ' gene models validated, using for training Augustus')
        
            ###Run Augustus training
            GFF2GB = os.path.join(AUGUSTUS_BASE, 'scripts', 'gff2gbSmallDNA.pl')
            trainingset = os.path.join(args.out, 'predict_misc', 'busco.training.gb')
            subprocess.call([GFF2GB, busco_final, MaskGenome, '250', trainingset], stderr = FNULL)
            if args.optimize_augustus:
                lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset, MaskGenome, args.out, args.cpus, True)
            else:
                lib.trainAugustus(AUGUSTUS_BASE, aug_species, trainingset, MaskGenome, args.out, args.cpus, False)

        lib.log.info("Running Augustus gene prediction")
        #now run Augustus multithreaded...
        if not os.path.isfile(aug_out):
            if os.path.isfile(hints_all):
                subprocess.call([AUGUSTUS_PARALELL, '--species', aug_species, '--hints', hints_all, '-i', MaskGenome, '-o', aug_out, '--cpus', str(args.cpus)], stderr = FNULL, stdout = FNULL)
            else:
                subprocess.call([AUGUSTUS_PARALELL, '--species', aug_species, '-i', MaskGenome, '-o', aug_out, '--cpus', str(args.cpus)], stderr = FNULL, stdout = FNULL)

        Augustus = os.path.join(args.out, 'predict_misc', 'augustus.evm.gff3')
        with open(Augustus, 'w') as output:
            subprocess.call(['perl', Converter, aug_out], stdout = output, stderr = FNULL)
        
        
    #just double-check that you've gotten here and both Augustus/GeneMark are finished
    if not any([Augustus, GeneMark]):
        lib.log.error("Augustus or GeneMark prediction is missing, check log files for errors")
        os._exit(1)
    
    #GeneMark can fail if you try to pass a single contig, check file length
    GM_check = lib.line_count(GeneMark)
    gmc = 1
    if GM_check < 3:
        gmc = 0
        lib.log.error("GeneMark predictions failed, proceeding with only Augustus")
    
    #if hints used for Augustus, get high quality models > 80% coverage to pass to EVM
    if os.path.isfile(hints_all):
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
                    if line.startswith('# % of transcript supported by hints'):
                        support = line.split(' ')[-1]
                        values.append(support)
                if float(values[1]) > 65: #greater than ~66% of exons supported, i.e. 2/3 or 3/4, but not less
                    hiQ_models.append(values[0])

        #now open evm augustus and pull out models
        HiQ = set(hiQ_models)
        lib.log.info("Found %i high quality predictions from Augustus (>66%% exon evidence)"  % len(HiQ))
        HiQ_match = re.compile(r'\b(?:%s)[\.t1;]+\b' % '|'.join(HiQ))
        AugustusHiQ = os.path.join(args.out, 'predict_misc', 'augustus-HiQ.evm.gff3')
        with open(AugustusHiQ, 'w') as HiQ_out:
            with open(Augustus, 'rU') as evm_aug:
                for line in evm_aug:
                    if HiQ_match.search(line):
                        newline = line.replace('\tAugustus\t', '\tHiQ\t')
                        HiQ_out.write(newline)


    #EVM related input tasks, find all predictions and concatenate together
    if args.pasa_gff:
        if os.path.isfile(hints_all):
            pred_in = [Augustus, GeneMark, args.pasa_gff, AugustusHiQ]
        else:
            pred_in = [Augustus, GeneMark, args.pasa_gff]
    else:
        if os.path.isfile(hints_all):
            pred_in = [Augustus, GeneMark, AugustusHiQ]
        else:
            pred_in = [Augustus, GeneMark]
    Predictions = os.path.join(args.out, 'predict_misc', 'predictions.gff3')
    with open(Predictions, 'w') as output:
        for f in pred_in:
            with open(f) as input:
                output.write(input.read())

    #set Weights file dependent on which data is present.
    Weights = os.path.join(args.out, 'predict_misc', 'weights.evm.txt')
    with open(Weights, 'w') as output:
        output.write("ABINITIO_PREDICTION\tAugustus\t1\n")
        output.write("ABINITIO_PREDICTION\tGeneMark\t1\n")
        if os.path.isfile(hints_all):
            output.write("OTHER_PREDICTION\tHiQ\t10\n")
        if args.pasa_gff:
            output.write("OTHER_PREDICTION\ttransdecoder\t10\n")
        if exonerate_out:
            output.write("PROTEIN\texonerate\t1\n")
        if Transcripts:
            output.write("TRANSCRIPT\tgenome\t1\n")

#total up Predictions
total = lib.countGFFgenes(Predictions)
lib.log.info('{0:,}'.format(total) + ' total gene models from all sources')

#setup EVM run
EVM_out = os.path.join(args.out, 'predict_misc', 'evm.round1.gff3')
EVM_script = os.path.join(parentdir, 'bin', 'funannotate-runEVM.py')
#get absolute paths for everything
Weights = os.path.abspath(Weights)
EVM_out = os.path.abspath(EVM_out)
Predictions = os.path.abspath(Predictions)

#parse entire EVM command to script
if Exonerate and Transcripts:
    evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM.log'), str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--protein_alignments', Exonerate, '--transcript_alignments', Transcripts, '--weights', Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
elif not Exonerate and Transcripts:
    evm_cmd = [sys.executable, EVM_script, os.path.join(args.out, 'logfiles', 'funannotate-EVM.log'),str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--transcript_alignments', Transcripts, '--weights', Weights, '--min_intron_length', str(args.min_intronlen), EVM_out]
elif not Transcripts and Exonerate:
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
    os._exit(1)
#check number of gene models, if 0 then failed, delete output file for re-running
if total < 1:
    lib.log.error("Evidence modeler has failed, exiting")
    os.remove(EVM_out)
    os._exit(1)
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

#combine tRNAscan with EVM gff
lib.log.info("Merging EVM output with tRNAscan output")
gffs = [tRNAscan, EVM_out]
GFF = os.path.join(args.out, 'predict_misc', 'evm.trnascan.gff')
with open(GFF, 'w') as output:
    for f in gffs:
        with open(f) as input:
            output.write(input.read())

#run GAG to get gff and proteins file for screening
lib.log.info("Reformatting GFF file using GAG")
subprocess.call(['gag.py', '-f', MaskGenome, '-g', GFF, '-o', 'gag1','--fix_start_stop', '-ril', str(args.max_intronlen)], stdout = FNULL, stderr = FNULL)
GAG_gff = os.path.join('gag1', 'genome.gff')
GAG_proteins = os.path.join('gag1', 'genome.proteins.fasta')
total = lib.countGFFgenes(GAG_gff)
lib.log.info('{0:,}'.format(total) + ' total gene models')

#filter bad models
lib.log.info("Filtering out bad gene models (< %i aa in length, transposable elements, etc)." % (args.min_protlen))
Blast_rep_remove = os.path.join(args.out, 'predict_misc', 'repeat.gene.models.txt')
if not os.path.isfile(Blast_rep_remove):
    lib.RepeatBlast(GAG_proteins, args.cpus, 1e-10, os.path.join(args.out, 'predict_misc'), Blast_rep_remove)
CleanGFF = os.path.join(args.out, 'predict_misc', 'cleaned.gff3')
lib.RemoveBadModels(GAG_proteins, GAG_gff, args.min_protlen, RepeatMasker, Blast_rep_remove, os.path.join(args.out, 'predict_misc'), CleanGFF) 
total = lib.countGFFgenes(CleanGFF)
lib.log.info('{0:,}'.format(total) + ' gene models remaining')

#need to write to tbl2asn twice to fix errors, run first time and then parse error report
lib.log.info("Converting to preliminary Genbank format")
subprocess.call(['gag.py', '-f', MaskGenome, '-g', CleanGFF, '-o', 'gag2','--fix_start_stop'], stdout = FNULL, stderr = FNULL)
shutil.copyfile(os.path.join('gag2', 'genome.fasta'), os.path.join('gag2', 'genome.fsa'))
SBT = os.path.join(parentdir, 'lib', 'test.sbt')
discrep = 'discrepency.report.txt'
if args.isolate:
    ORGANISM = "[organism=" + args.species + "] " + "[isolate=" + args.isolate + "]"
else:
    ORGANISM = "[organism=" + args.species + "]"
subprocess.call(['tbl2asn', '-p', 'gag2', '-t', SBT, '-M', 'n', '-Z', discrep, '-a', 'r10u', '-l', 'paired-ends', '-j', ORGANISM, '-V', 'b', '-c', 'fx'], stdout = FNULL, stderr = FNULL)

#now parse error reports and remove bad models
lib.log.info("Cleaning models flagged by tbl2asn")
NCBIcleanGFF = os.path.join(args.out, 'predict_misc', 'ncbi.cleaned.gff3')
ErrSum = os.path.join('gag2', 'errorsummary.val')
Val = os.path.join('gag2', 'genome.val')
DirtyGFF = os.path.join('gag2', 'genome.gff')
if args.keep_no_stops:
    lib.ParseErrorReport(DirtyGFF, ErrSum, Val, discrep, NCBIcleanGFF, keep_stops=True)
else:
    lib.ParseErrorReport(DirtyGFF, ErrSum, Val, discrep, NCBIcleanGFF, keep_stops=False)
total = lib.countGFFgenes(NCBIcleanGFF)
lib.log.info('{0:,}'.format(total) + ' gene models remaining')
shutil.copyfile(discrep, os.path.join('gag2', discrep))

#now we can rename gene models
lib.log.info("Re-naming gene models")
MAP = os.path.join(parentdir, 'util', 'maker_map_ids.pl')
MAPGFF = os.path.join(parentdir, 'util', 'map_gff_ids.pl')
mapping = os.path.join(args.out, 'predict_misc', 'mapping.ids')
if not args.name.endswith('_'):
    args.name = args.name + '_'
with open(mapping, 'w') as output:
    subprocess.call(['perl', MAP, '--prefix', args.name, '--justify', '5', '--suffix', '-T', '--iterate', '1', NCBIcleanGFF], stdout = output, stderr = FNULL)
subprocess.call(['perl', MAPGFF, mapping, NCBIcleanGFF], stdout = FNULL, stderr = FNULL)

#run GAG again with clean dataset, fix start/stops
subprocess.call(['gag.py', '-f', MaskGenome, '-g', NCBIcleanGFF, '-o', 'tbl2asn', '--fix_start_stop'], stdout = FNULL, stderr = FNULL)

#setup final output files
base = args.species.replace(' ', '_').lower()
final_fasta = os.path.join(args.out, 'predict_results', base + '.scaffolds.fa')
final_gff = os.path.join(args.out, 'predict_results', base + '.gff3')
final_gbk = os.path.join(args.out, 'predict_results', base + '.gbk')
final_tbl = os.path.join(args.out, 'predict_results', base + '.tbl')
final_proteins = os.path.join(args.out, 'predict_results', base + '.proteins.fa')

#run tbl2asn in new directory directory
shutil.copyfile(os.path.join('tbl2asn', 'genome.fasta'), os.path.join('tbl2asn', 'genome.fsa'))
discrep = os.path.join(args.out, 'predict_results', base + '.discrepency.report.txt')
lib.log.info("Converting to final Genbank format")
subprocess.call(['tbl2asn', '-p', 'tbl2asn', '-t', SBT, '-M', 'n', '-Z', discrep, '-a', 'r10u', '-l', 'paired-ends', '-j', ORGANISM, '-V', 'b', '-c', 'fx'], stdout = FNULL, stderr = FNULL)
shutil.copyfile(os.path.join('tbl2asn', 'genome.fasta'), final_fasta)
shutil.copyfile(os.path.join('tbl2asn', 'genome.gff'), final_gff)
shutil.copyfile(os.path.join('tbl2asn', 'genome.gbf'), final_gbk)
shutil.copyfile(os.path.join('tbl2asn', 'genome.tbl'), final_tbl)
lib.log.info("Collecting final annotation files")

lib.log.info("Funannotate predict is finished, output files are in the %s/predict_results folder" % (args.out))
lib.log.info("Note, you should pay attention to any tbl2asn errors now before running functional annotation, although many automatic steps were taken to ensure NCBI submission compatibility, it is likely that some manual editing will be required.")

#clean up intermediate folders
if os.path.isdir('genemark_gag'):
    shutil.rmtree('genemark_gag')
if os.path.isdir('genemark'):
    if os.path.isdir(os.path.join(args.out, 'predict_misc', 'genemark')):
        shutil.rmtree(os.path.join(args.out, 'predict_misc', 'genemark'))
    os.rename('genemark', os.path.join(args.out, 'predict_misc', 'genemark'))
if os.path.isdir('gag1'):
    if os.path.isdir(os.path.join(args.out, 'predict_misc', 'gag1')):
        shutil.rmtree(os.path.join(args.out, 'predict_misc', 'gag1'))
    os.rename('gag1', os.path.join(args.out, 'predict_misc', 'gag1'))
if os.path.isdir('gag2'):
    if os.path.isdir(os.path.join(args.out, 'predict_misc', 'gag2')):
        shutil.rmtree(os.path.join(args.out, 'predict_misc', 'gag2'))
    os.rename('gag2', os.path.join(args.out, 'predict_misc', 'gag2'))
if os.path.isfile('discrepency.report.txt'):
    os.rename('discrepency.report.txt', os.path.join('tbl2asn', 'discrepency.report.txt'))
if os.path.isdir('tbl2asn'):
    if os.path.isdir(os.path.join(args.out, 'predict_misc', 'tbl2asn')):
        shutil.rmtree(os.path.join(args.out, 'predict_misc', 'tbl2asn'))
    os.rename('tbl2asn', os.path.join(args.out, 'predict_misc', 'tbl2asn'))
if os.path.isfile('funannotate-EVM.log'):
    os.rename('funannotate-EVM.log', os.path.join(args.out, 'logfiles', 'funannotate-EVM.log'))
if os.path.isfile('funannotate-p2g.log'):
    os.rename('funannotate-p2g.log', os.path.join(args.out, 'logfiles', 'funannotate-p2g.log'))
os._exit(1)
