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
parser=argparse.ArgumentParser(prog='funannotate-predict.py', usage="%(prog)s [options] -i genome.fasta",
    description='''Script that does it all...''',
    epilog="""Written by Jon Palmer (2016) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('-s','--species', required=True, help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('--isolate', help='Isolate/strain name (e.g. Af293)')
parser.add_argument('-n','--name', default="FUN_", help='Shortname for genes, perhaps assigned by NCBI, eg. VC83')
parser.add_argument('--augustus_species', help='Specify species for Augustus')
parser.add_argument('--genemark_mod', help='Use pre-existing Genemark training file (e.g. gmhmm.mod)')
parser.add_argument('--protein_evidence', default='uniprot.fa', help='Specify protein evidence (multiple files can be separaed by a comma)')
parser.add_argument('--exonerate_proteins', help='Pre-computed Exonerate protein alignments')
parser.add_argument('--transcript_evidence', help='Transcript evidence (map to genome with GMAP)')
parser.add_argument('--gmap_gff', help='Pre-computed GMAP transcript alignments (GFF3)')
parser.add_argument('--pasa_gff', help='Pre-computed PASA/TransDecoder high quality models')
parser.add_argument('--augustus_gff', help='Pre-computed Augustus gene models (GFF3)')
parser.add_argument('--genemark_gtf', help='Pre-computed GeneMark gene models (GTF)')
parser.add_argument('--repeatmodeler_lib', help='Pre-computed RepeatModeler (or other) repetitive elements')
parser.add_argument('--rna_bam', help='BAM (sorted) of RNAseq aligned to reference for BRAKER1')
parser.add_argument('--min_intron_length', default=10, help='Minimum intron length for gene models')
parser.add_argument('--max_intron_length', default=3000, help='Maximum intron length for gene models')
parser.add_argument('--min_protein_length', default=51, type=int, help='Minimum amino acid length for valid gene model')
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use')
args=parser.parse_args()

#create log file
log_name = 'funannotate-predict.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.log.info("Operating system: %s, %i cores, ~ %i GB RAM" % (sys.platform, multiprocessing.cpu_count(), lib.MemoryCheck()))

programs = ['hmmscan','blastp','blastn','gag.py','tbl2asn','runiprscan','gmes_petap.pl', 'BuildDatabase', 'RepeatModeler', 'RepeatMasker', 'genemark_gtf2gff3','autoAug.pl', 'maker', 'bedtools', 'gmap', 'gmap_build']
lib.CheckDependencies(programs)

#create temp folder
if not os.path.exists(args.out):
    os.makedirs(args.out)

#do some checks and balances
try:
    EVM = os.environ["EVM_HOME"]
except KeyError:
    lib.log.error("$EVM_HOME enironmental variable not found, either Evidence Modeler is not installed or variable not in $PATH")
    os._exit(1)
    
#check augustus species now, so that you don't get through script and then find out it is already in DB
if not args.augustus_species:
    aug_species = args.species.replace(' ', '_').lower()
else:
    aug_species = args.augustus_species
if lib.CheckAugustusSpecies(aug_species):
    lib.log.error("Augustus training set for %s already exists, thus funannotate will use those parameters. If this is not what you want, exit script and provide a unique name for the --augustus_species argument" % (aug_species))
    
if args.protein_evidence == 'uniprot.fa':
    args.protein_evidence = os.path.join(script_path, 'DB', 'uniprot_sprot.fasta')

#EVM command line scripts
Converter = os.path.join(EVM, 'EvmUtils', 'misc', 'augustus_GFF3_to_EVM_GFF3.pl')
ExoConverter = os.path.join(EVM, 'EvmUtils', 'misc', 'exonerate_gff_to_alignment_gff3.pl')
Validator = os.path.join(EVM, 'EvmUtils', 'gff3_gene_prediction_file_validator.pl')

#so first thing is to reformat genome fasta only if there is no aligned evidence already
if not any([args.gmap_gff, args.pasa_gff, args.augustus_gff, args.genemark_gtf, args.rna_bam, args.exonerate_proteins]):
    #reformat fasta headers to avoid problems with Augustus
    lib.log.info("Re-formatting genome FASTA headers")
    sort_out = os.path.join(args.out, 'genome.fasta')
    lib.SortRenameHeaders(args.input, sort_out)
    Genome = os.path.abspath(sort_out)
else:
    Genome = os.path.abspath(args.input)

#repeatmasker, run if not passed from command line
if not args.repeatmodeler_lib:
    MaskGenome = os.path.join(args.out, 'genome.softmasked.fa')
    if not os.path.isfile(MaskGenome):
        lib.RepeatModelMask(Genome, args.cpus, args.out, MaskGenome)
else:
    MaskGenome = os.path.join(args.out, 'genome.softmasked.fa')
    if not os.path.isfile(MaskGenome):
        lib.RepeatMask(Genome, args.repeatmodeler_lib, args.cpus, args.out, MaskGenome)
RepeatMasker = os.path.join(args.out, 'repeatmasker.gff3')
RepeatMasker = os.path.abspath(RepeatMasker)
MaskGenome = os.path.abspath(MaskGenome)

#check for transcript evidence/format as needed
if not args.gmap_gff:
    if args.transcript_evidence:
        if ',' in args.transcript_evidence:
            trans_temp = os.path.join(args.out, 'transcripts.combined.fa')
            files = args.transcript_evidence.split(",")
            with open(trans_temp, 'w') as output:
                for f in files:
                    with open(f) as input:
                        output.write(input.read())
        else:
            trans_temp = args.transcript_evidence        
        #run Gmap of transcripts to genome
        trans_out = os.path.join(args.out, 'transcript_alignments.gff3')
        lib.log.info("Aligning transcript evidence to genome with GMAP")
        if not os.path.isfile(trans_out):
            lib.runGMAP(trans_temp, MaskGenome, args.cpus, args.max_intron_length, args.out, trans_out)
        Transcripts = os.path.abspath(trans_out)
    else:
        Transcripts = False
else:
    Transcripts = os.path.abspath(args.gmap_gff)

#check for protein evidence/format as needed
if not args.exonerate_proteins:
    if args.protein_evidence:
        if ',' in args.protein_evidence:
            prot_temp = os.path.join(args.out, 'proteins.combined.fa')
            files = args.protein_evidence.split(",")
            with open(prot_temp, 'w') as output:
                for f in files:
                    with open(f) as input:
                        output.write(input.read())
        else:
            prot_temp = args.protein_evidence
        #run funannotate-p2g to map to genome
        lib.log.info("Mapping proteins to genome using tBlastn/Exonerate")
        P2G = os.path.join(script_path,'funannotate-p2g.py')
        p2g_out = os.path.join(args.out, 'exonerate.out')
        p2g_cmd = [sys.executable, P2G, prot_temp, MaskGenome, p2g_out, str(args.max_intron_length), str(args.cpus)]
        if not os.path.isfile(p2g_out):
            subprocess.call(p2g_cmd)
        exonerate_out = os.path.abspath(p2g_out)
    else:
        exonerate_out = False
else:
    exonerate_out = os.path.abspath(args.exonerate_proteins)

if exonerate_out:
    Exonerate = os.path.join(args.out, 'protein_alignments.gff3')
    with open(Exonerate, 'w') as output:
        subprocess.call([ExoConverter, exonerate_out], stdout = output, stderr = FNULL)
    Exonerate = os.path.abspath(Exonerate)

Augustus = ''
GeneMark = ''

#Walk thru data available and determine best approach. 
if args.genemark_gtf:
    #convert the predictors to EVM format and merge
    #convert GeneMark
    GeneMarkGFF3 = os.path.join(args.out, 'genemark.gff')
    with open(GeneMarkGFF3, 'w') as output:
        subprocess.call(['genemark_gtf2gff3', args.genemark_gtf], stdout = output, stderr = FNULL)
    GeneMarkTemp = os.path.join(args.out, 'genemark.temp.gff')
    with open(GeneMarkTemp, 'w') as output:
        subprocess.call(['perl', Converter, GeneMarkGFF3], stdout = output, stderr = FNULL)
    GeneMark = os.path.join(args.out, 'genemark.evm.gff3')
    with open(GeneMark, 'w') as output:
        with open(GeneMarkTemp, 'rU') as input:
            lines = input.read().replace("Augustus","GeneMark")
            output.write(lines)

if args.augustus_gff:
    #convert Augustus
    Augustus = os.path.join(args.out, 'augustus.evm.gff3')
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
    lib.log.info("Now launching BRAKER to train GeneMark and Augustus")
    species = '--species=' + aug_species
    genome = '--genome=' + MaskGenome
    bam = '--bam=' + os.path.abspath(args.rna_bam)
    with open('braker.log') as logfile:
        subprocess.call(['braker.pl', '--fungus', '--cores', str(args.cpus), '--gff3', '--softmasking', '1', genome, species, bam], stdout = logfile, stderr = logfile)
    #okay, now need to fetch the Augustus GFF and Genemark GTF files
    aug_out = os.path.join('braker', aug_species, 'augustus.gff3')
    gene_out = os.path.join('braker', aug_species, 'GeneMark-ET', 'genemark.gtf')
    #now convert to EVM format
    Augustus = os.path.join(args.out, 'augustus.evm.gff3')
    with open(Augustus, 'w') as output:
        subprocess.call(['perl', Converter, aug_out], stdout = output, stderr = FNULL)
    GeneMarkGFF3 = os.path.join(args.out, 'genemark.gff')
    with open(GeneMarkGFF3, 'w') as output:
        subprocess.call(['genemark_gtf2gff3', gene_out], stdout = output, stderr = FNULL)
    GeneMarkTemp = os.path.join(args.out, 'genemark.temp.gff')
    with open(GeneMarkTemp, 'w') as output:
        subprocess.call(['perl', Converter, GeneMarkGFF3], stdout = output, stderr = FNULL)
    GeneMark = os.path.join(args.out, 'genemark.evm.gff3')
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
    #input params
    species = '--species=' + aug_species
    genome = '--genome=' + MaskGenome
    aug_out = os.path.join(args.out, 'augustus.gff3')
    #check for training data
    if lib.CheckAugustusSpecies(aug_species):
        lib.log.error("%s as already been trained, using existing parameters" % (aug_species))
        lib.log.info("Running Augustus gene prediction")
        aug_out = os.path.join(args.out, 'augustus.gff3')
        species = '--species=' + aug_species
        if lib.CheckAugustusSpecies(aug_species):
            lib.log.error("%s as already been trained, using existing parameters" % (aug_species))
            lib.log.info("Running Augustus gene prediction")
            if not os.path.isfile(aug_out):
                with open(aug_out, 'w') as output:
                    subprocess.call(['augustus', species, '--gff3=on', MaskGenome], stdout = output, stderr = FNULL)
            Augustus = os.path.join(args.out, 'augustus.evm.gff3')
            with open(Augustus, 'w') as output:
                subprocess.call(['perl', Converter, aug_out], stdout = output, stderr = FNULL)
            with open(aug_out, 'w') as output:
                subprocess.call(['augustus', species, '--gff3=on', genome], stdout = output, stderr = FNULL)
    else:
        lib.log.info("Training augustus using PASA data, this may take awhile")
        training = '--trainingset=' + args.pasa_gff
        aug_log = os.path.join(args.out, 'augustus.log')
        with open(aug_log, 'w') as logfile:
            subprocess.call(['autoAug.pl', '--noutr', '--singleCPU', species, genome, training], stdout=logfile, stderr=logfile)
    #get routine for formatting augustus results
    #
    #    
    
if not GeneMark:
    #now run GeneMark-ES, first check for gmhmm mod file, use if available otherwise run ES
    if not args.genemark_mod:
        GeneMarkGFF3 = os.path.join(args.out, 'genemark.gff')
        if not os.path.isfile(GeneMarkGFF3):
            lib.RunGeneMarkES(MaskGenome, args.cpus, args.out, GeneMarkGFF3)
        GeneMarkTemp = os.path.join(args.out, 'genemark.temp.gff')
        with open(GeneMarkTemp, 'w') as output:
            subprocess.call(['perl', Converter, GeneMarkGFF3], stdout = output, stderr = FNULL)
        GeneMark = os.path.join(args.out, 'genemark.evm.gff3')
        with open(GeneMark, 'w') as output:
            with open(GeneMarkTemp, 'rU') as input:
                lines = input.read().replace("Augustus","GeneMark")
                output.write(lines)
    else:   #have training parameters file, so just run genemark with
        GeneMarkGFF3 = os.path.join(args.out, 'genemark.gff')
        if not os.path.isfile(GeneMarkGFF3):
            lib.RunGeneMark(MaskGenome, args.genemark_mod, args.cpus, args.out, GeneMarkGFF3)
        GeneMarkTemp = os.path.join(args.out, 'genemark.temp.gff')
        with open(GeneMarkTemp, 'w') as output:
            subprocess.call(['perl', Converter, GeneMarkGFF3], stdout = output, stderr = FNULL)
        GeneMark = os.path.join(args.out, 'genemark.evm.gff3')
        with open(GeneMark, 'w') as output:
            with open(GeneMarkTemp, 'rU') as input:
                lines = input.read().replace("Augustus","GeneMark")
                output.write(lines)

if not Augustus: 
    if not args.augustus_species:
        aug_species = args.species.replace(' ', '_').lower()
    else:
        aug_species = args.augustus_species
    aug_out = os.path.join(args.out, 'augustus.gff3')
    species = '--species=' + aug_species
    if lib.CheckAugustusSpecies(aug_species):
        lib.log.error("%s as already been trained, using existing parameters" % (aug_species))
        lib.log.info("Running Augustus gene prediction")
        if not os.path.isfile(aug_out):
            with open(aug_out, 'w') as output:
                subprocess.call(['augustus', species, '--gff3=on', MaskGenome], stdout = output, stderr = FNULL)
        Augustus = os.path.join(args.out, 'augustus.evm.gff3')
        with open(Augustus, 'w') as output:
            subprocess.call(['perl', Converter, aug_out], stdout = output, stderr = FNULL)
 
    else: #run fCEGMA then train Augustus
        #run GAG to get protein sequences
        lib.log.info("Prepping data using GAG")
        subprocess.call(['gag.py', '-f', MaskGenome, '-g', GeneMarkGFF3, '--fix_start_stop', '-o', 'genemark_gag'], stdout = FNULL, stderr = FNULL)
        os.rename(os.path.join('genemark_gag', 'genome.proteins.fasta'), os.path.join(args.out, 'genemark.proteins.fasta'))
        #filter using fCEGMA models
        lib.log.info("Now filtering best fCEGMA models for training Augustus")
        fCEGMA_in = os.path.join(args.out, 'genemark.proteins.fasta')
        fCEGMA_out = os.path.join(args.out, 'fCEGMA_hits.txt')
        lib.fCEGMA(fCEGMA_in, args.cpus, 1e-100, args.out, GeneMarkGFF3, fCEGMA_out)

        #now run Augustus training based on training set
        fCEGMA_gff = os.path.join(args.out, 'training.gff3')
        fCEGMA_gff = os.path.abspath(fCEGMA_gff)
        #check number of models
        total = lib.countGFFgenes(fCEGMA_gff)
        if total < 100:
            lib.log.error("Number of training models %i is too low, need at least 100" % total)
            os._exit(1)
        lib.log.info("Now training Augustus with fCEGMA filtered dataset")
        if os.path.exists('autoAug'):
            shutil.rmtree('autoAug')
        species = '--species=' + aug_species
        genome = '--genome=' + MaskGenome
        training = '--trainingset=' + fCEGMA_gff
        aug_log = os.path.join(args.out, 'augustus.log')
        with open(aug_log, 'w') as logfile:
            subprocess.call(['autoAug.pl', '--noutr', '--singleCPU', species, genome, training], stdout=logfile, stderr=logfile)
        #get routine for formatting augustus results
        #
        #    

#just double-check that you've gotten here and both Augustus/GeneMark are finished
if not any([Augustus, GeneMark]):
    lib.log.error("Augustus or GeneMark prediction is missing, check log files for errors")
    os._exit(1)

#EVM related input tasks, find all predictions and concatenate together
if args.pasa_gff:
    pred_in = [Augustus, GeneMark, args.pasa_gff]
else:
    pred_in = [Augustus, GeneMark]
Predictions = os.path.join(args.out, 'predictions.gff3')
with open(Predictions, 'w') as output:
    for f in pred_in:
        with open(f) as input:
            output.write(input.read())

#set Weights file dependent on which data is present.
Weights = os.path.join(args.out, 'weights.evm.txt')
with open(Weights, 'w') as output:
    output.write("ABINITIO_PREDICTION\tAugustus\t1\n")
    output.write("ABINITIO_PREDICTION\tGeneMark\t1\n")
    if args.pasa_gff:
        output.write("OTHER_PREDICTION\ttransdecoder\t10\n")
    if exonerate_out:
        output.write("PROTEIN\tspliced_protein_alignments\t1\n")
    if Transcripts:
        output.write("TRANSCRIPT\tspliced_transcript_alignments\t1\n")

#total up Predictions
total = lib.countGFFgenes(Predictions)
lib.log.info('{0:,}'.format(total) + ' total gene models from all sources')
#setup EVM run
EVM_out = os.path.join(args.out, 'evm.round1.gff3')
EVM_script = os.path.join(script_path, 'funannotate-runEVM.py')
#get absolute paths for everything
Weights = os.path.abspath(Weights)
EVM_out = os.path.abspath(EVM_out)
Predictions = os.path.abspath(Predictions)
#parse entire EVM command to script
if Exonerate and Transcripts:
    evm_cmd = [sys.executable, EVM_script, str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--protein_alignments', Exonerate, '--transcript_alignments', Transcripts, '--weights', Weights, '--min_intron_length', str(args.min_intron_length), EVM_out]
elif not Exonerate and Transcripts:
    evm_cmd = [sys.executable, EVM_script, str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--transcript_alignments', Transcripts, '--weights', Weights, '--min_intron_length', str(args.min_intron_length), EVM_out]
elif not Transcripts and Exonerate:
    evm_cmd = [sys.executable, EVM_script, str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--protein_alignments', Exonerate, '--weights', Weights, '--min_intron_length', str(args.min_intron_length), EVM_out]
elif not any([Transcripts,Exonerate]):
    evm_cmd = [sys.executable, EVM_script, str(args.cpus), '--genome', MaskGenome, '--gene_predictions', Predictions, '--weights', Weights, '--min_intron_length', str(args.min_intron_length), EVM_out]
#run EVM
if not os.path.isfile(EVM_out):
    subprocess.call(evm_cmd)
total = lib.countGFFgenes(EVM_out)
lib.log.info('{0:,}'.format(total) + ' total gene models from EVM')
#run tRNAscan
lib.log.info("Predicting tRNAs")
tRNAscan = os.path.join(args.out, 'trnascan.gff3')
lib.runtRNAscan(MaskGenome, args.out, tRNAscan)
total = lib.countGFFgenes(tRNAscan)
lib.log.info('{0:,}'.format(total) + ' tRNAs found')
#combine tRNAscan with EVM gff
lib.log.info("Merging EVM output with tRNAscan output")
gffs = [tRNAscan, EVM_out]
GFF = os.path.join(args.out, 'evm.trnascan.gff')
with open(GFF, 'w') as output:
    for f in gffs:
        with open(f) as input:
            output.write(input.read())
#run GAG to get gff and proteins file for screening
lib.log.info("Reformatting GFF file using GAG")
subprocess.call(['gag.py', '-f', MaskGenome, '-g', GFF, '-o', 'gag1','--fix_start_stop', '-ril', str(args.max_intron_length)], stdout = FNULL, stderr = FNULL)
GAG_gff = os.path.join('gag1', 'genome.gff')
GAG_proteins = os.path.join('gag1', 'genome.proteins.fasta')
total = lib.countGFFgenes(GAG_gff)
lib.log.info('{0:,}'.format(total) + ' total gene models')
#filter bad models
lib.log.info("Filtering out bad gene models (internal stops, transposable elements, etc).")
CleanGFF = os.path.join(args.out, 'cleaned.gff3')
lib.RemoveBadModels(GAG_proteins, GAG_gff, args.min_protein_length, RepeatMasker, args.out, CleanGFF) 
total = lib.countGFFgenes(CleanGFF)
lib.log.info('{0:,}'.format(total) + ' gene models remaining')
#now we can rename gene models
lib.log.info("Re-naming gene models")
MAP = os.path.join(script_path, 'util', 'maker_map_ids.pl')
MAPGFF = os.path.join(script_path, 'util', 'map_gff_ids.pl')
mapping = os.path.join(args.out, 'mapping.ids')
with open(mapping, 'w') as output:
    subprocess.call(['perl', MAP, '--prefix', args.name, '--justify', '5', '--suffix', '-T', '--iterate', '1', CleanGFF], stdout = output, stderr = FNULL)
subprocess.call(['perl', MAPGFF, mapping, CleanGFF], stdout = FNULL, stderr = FNULL)   
#run GAG again with clean dataset, fix start/stops
subprocess.call(['gag.py', '-f', MaskGenome, '-g', CleanGFF, '-o', 'tbl2asn', '--fix_start_stop'], stdout = FNULL, stderr = FNULL)
#setup final output files
base = args.species.replace(' ', '_').lower()
final_fasta = base + '.scaffolds.fa'
final_gff = base + '.gff3'
final_gbk = base + '.gbk'
final_tbl = base + '.tbl'
final_proteins = base + '.proteins.fa'
final_smurf = base + '.smurf.txt'
#run tbl2asn in new directory directory
shutil.copyfile(os.path.join('tbl2asn', 'genome.fasta'), os.path.join('tbl2asn', 'genome.fsa'))
lib.log.info("Converting to Genbank format")
SBT = os.path.join(script_path, 'lib', 'test.sbt')
if args.isolate:
    ORGANISM = "[organism=" + args.species + "] " + "[isolate=" + args.isolate + "]"
else:
    ORGANISM = "[organism=" + args.species + "]"
subprocess.call(['tbl2asn', '-p', 'tbl2asn', '-t', SBT, '-M', 'n', '-Z', 'discrepency.report.txt', '-a', 'r10u', '-l', 'paired-ends', '-j', ORGANISM, '-V', 'b', '-c', 'fx'], stdout = FNULL, stderr = FNULL)
shutil.copyfile(os.path.join('tbl2asn', 'genome.fasta'), final_fasta)
shutil.copyfile(os.path.join('tbl2asn', 'genome.gff'), final_gff)
shutil.copyfile(os.path.join('tbl2asn', 'genome.gbf'), final_gbk)
shutil.copyfile(os.path.join('tbl2asn', 'genome.tbl'), final_tbl)
lib.log.info("Collecting final annotation files")

#Create AGP and contigs
lib.log.info("Creating AGP file and corresponding contigs file")
agp2fasta = os.path.join(script_path, 'util', 'fasta2agp.pl')
AGP = base + '.agp'
with open(AGP, 'w') as output:
    subprocess.call(['perl', agp2fasta, final_fasta], stdout = output, stderr = FNULL)
#run gb2smurf here so user can run secondary metabolite prediction for annotation
lib.log.info("Creating input files for SMURF server")
lib.gb2smurf(final_gbk, final_proteins, final_smurf)
lib.log.info("Funannotate predict is finished, final output files have %s base name in this directory" % (base))
lib.log.info("Note, you should pay attention to any tbl2asn errors now before running functional annotation, there is likely some manual editing required.")

#clean up intermediate folders
if os.path.isfile('genemark_gag'):
    shutil.rmtree('genemark_gag')
if os.path.isfile('genemark'):
    os.rename('genemark', os.path.join(args.out, 'genemark'))
if os.path.isfile('gag1'):
    shutil.rmtree('gag1')
if os.path.isfile('RepeatModeler'):
    os.rename('RepeatModeler', os.path.join(args.out, 'RepeatModeler'))
if os.path.isfile('RepeatMasker'):
    os.rename('RepeatMasker', os.path.join(args.out, 'RepatMasker'))
os._exit(1)
