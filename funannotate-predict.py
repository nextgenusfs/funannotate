#!/usr/bin/env python

import sys, os, subprocess,inspect, multiprocessing, shutil, argparse, time, fileinput
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
    description='''Script that does it all.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('-s','--species', required=True, help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('-n','--name', help='Shortname for genes, perhaps assigned by NCBI, eg. VC83')
parser.add_argument('--pipeline', choices=['rnaseq', 'fCEGMA', 'no_train', 'no_augustus_train', 'only_annotate_proteins'], help='Method to employ, BRAKER1, GeneMark/fCEGMA, no_training=')
parser.add_argument('--augustus_species', help='Specify species for Augustus')
parser.add_argument('--genemark_mod', help='Use pre-existing Genemark training file (e.g. gmhmm.mod)')
parser.add_argument('--protein_evidence', default='uniprot_sprot.fa', help='Specify protein evidence (multiple files can be separaed by a comma)')
parser.add_argument('--transcript_evidence', help='Specify transcript evidence')
parser.add_argument('--pasa_gff', help='Pre-computed PASA/TransDecoder high quality models')
parser.add_argument('--augustus_gff', help='Pre-computed Augustus gene models (GFF3)')
parser.add_argument('--genemark_gtf', help='Pre-computed GeneMark gene models (GTF)')
parser.add_argument('--exonerate_proteins', help='Pre-computed Exonerate protein alignments')
parser.add_argument('--repeatmasker', help='Pre-computed RepeatMasker out file')
parser.add_argument('--rna_bam', help='BAM (sorted) of RNAseq aligned to reference for BRAKER1')
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use')
args=parser.parse_args()

#create log file
log_name = args.out + '.funannotate-predict.log'
if os.path.isfile(log_name):
    os.remove(log_name)

#initialize script, log system info and cmd issue at runtime
lib.setupLogging(log_name)
FNULL = open(os.devnull, 'w')
cmd_args = " ".join(sys.argv)+'\n'
lib.log.debug(cmd_args)
print "-------------------------------------------------------"
lib.log.info("Operating system: %s, %i cores, %i GB RAM" % (sys.platform, multiprocessing.cpu_count(), lib.MemoryCheck()))

programs = ['hmmscan','blastp','blastn','gag.py','tbl2asn','runiprscan','gmes_petap.pl', 'BuildDatabase', 'RepeatModeler', 'RepeatMasker', 'genemark_gtf2gff3','autoAug.pl', 'maker', 'bedtools']
lib.CheckDependencies(programs)

#create temp folder
if not os.path.exists(args.out):
    os.makedirs(args.out)

#do some checks and balances
try:
    EVM = os.environ["EVM_HOME"]
except KeyError:
    lib.log.error("$EVM_HOME enironmental variable not found, either Evidence Modler is not installed or variable not in PATH")
    os._exit(1)

#alter the pipeline based on input args
if args.augustus_gff and args.genemark_gtf and args.pasa_gff and args.exonerate_proteins and args.repeatmasker: #all heavy lifting done, so run EVM and filter
    lib.log.info("Provided Augustus, GeneMark, PASA, Exonerate, and RepeatMasking. Running FunAnnotate accordingly...")
    Converter = os.path.join(EVM, 'EvmUtils', 'misc', 'augustus_GFF3_to_EVM_GFF3.pl')
    ProtConverter = os.path.join(EVM, 'EvmUtils', 'misc', 'exonerate_gff_to_alignment_gff3.pl')
    Validator = os.path.join(EVM, 'EvmUtils', 'gff3_gene_prediction_file_validator.pl')
    lib.log.info("Formatting input for Evidence Modeler")
    Augustus = os.path.join(args.out, 'augustus.evm.gff3')
    with open(Augustus, 'w') as output:
        subprocess.call(['perl', Converter, args.augustus_gff], stdout = output, stderr = FNULL)
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
    pred_in = [Augustus, GeneMark, args.pasa_gff]
    Predictions = os.path.join(args.out, 'predictions.gff3')
    with open(Predictions, 'w') as output:
        for f in pred_in:
            with open(f) as input:
                output.write(input.read())
    Exonerate = os.path.join(args.out, 'exonerate.evm.gff3')
    with open(Exonerate, 'w') as output:
        subprocess.call([ProtConverter, args.exonerate_proteins], stdout = output, stderr = FNULL)
    RepeatMasker = os.path.join(args.out, 'repeatmasked.gff3')
    with open(RepeatMasker, 'w') as output:
        subprocess.call(['rmOutToGFF3.pl', args.repeatmasker], stdout = output, stderr = FNULL)
    Weights = os.path.join(args.out, 'weights.evm.txt')
    with open(Weights, 'w') as output:
        output.write("ABINITIO_PREDICTION\tAugustus\t1\n")
        output.write("ABINITIO_PREDICTION\tGeneMark\t1\n")
        output.write("OTHER_PREDICTION\ttransdecoder\t10\n")
        output.write("PROTEIN\tspliced_protein_alignments\t1\n")

    #run EVM
    EVM_out = os.path.join(args.out, 'evm.round1.gff3')
    EVM_script = os.path.join(script_path, 'funannotate-runEVM.py')
    subprocess.call([sys.executable, EVM_script, args.input, Predictions, Exonerate, Weights, str(args.cpus), EVM_out])
    #run tRNAscan
    lib.log.info("Predicting tRNAs")
    tRNAscan = os.path.join(args.out, 'trnascan.gff3')
    lib.runtRNAscan(args.input, args.out, tRNAscan)
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
    subprocess.call(['gag.py', '-f', args.input, '-g', GFF, '-o', 'gag1'], stdout = FNULL, stderr = FNULL)
    GAG_gff = os.path.join('gag1', 'genome.gff')
    GAG_proteins = os.path.join('gag1', 'genome.proteins.fasta')
    #filter bad models
    CleanGFF = os.path.join(args.out, 'cleaned.gff3')
    lib.RemoveBadModels(GAG_proteins, GAG_gff, 50, RepeatMasker, args.out, CleanGFF)
    os._exit(1)



#run some prelim checks on the input
if args.pipeline == 'rnaseq':
    if not args.augustus_species:
        aug_species = args.species.replace(' ', '_').lower()
    else:
        aug_species = args.augustus_species
    if lib.CheckAugustusSpecies(aug_species):
        lib.log.error("%s as already been trained, choose a unique species name" % (aug_species))
        os._exit(1)

    #check input, will need some RNAseq BAM to run this pipeline
    if not args.rna_bam:
        lib.log.error("You specified RNAseq as a pipeline, but did not provide aligned RNAseq reads in BAM format")
        os._exit(1)
    '''
    #can't reformat headers if already aligned to the reference - hopefully this isn't a problem later
    #start by reformatting fasta headers
    sorted_input = os.path.join(args.out, 'genome.sorted.fa')
    lib.SortRenameHeaders(args.input, sorted_input)
    '''
    '''
    #not supporting mapping here, do that before running funannotate
    #now map RNAseq reads to reference using hisat2
    #build reference database
    subprocess.call(['hisat2-build', '-p', str(args.cpus), sorted_input, sorted_input], stdout = FNULL, stderr = FNULL)
    bam_out = os.path.join(args.out, 'rnaseq.sorted')
    hisat_log = os.path.join(args.out, 'hisat2.log')
    #organize the input reads
    if args.forward:
        FORWARD = '-1 ' + args.forward
    if args.reverse:
        REVERSE = '-2 ' + args.reverse
    if args.single:
        SINGLE = '-U ' + args.single
    with open(hisat_log, 'w') as logfile:
        if not args.single:
            subprocess.call(['hisat2', '-p', str(args.cpus), FORWARD, REVERSE, '|', 'samtools', 'view', '-@', str(args.cpus), '-bS', '-', '|', 'samtools', 'sort', '-@', str(args.cpus), '-', bam_out], stdout = FNULL, stderr = logfile)
        elif not args.forward:
            subprocess.call(['hisat2', '-p', str(args.cpus), SINGLE, '|', 'samtools', 'view', '-@', str(args.cpus), '-bS', '-', '|', 'samtools', 'sort', '-@', str(args.cpus), '-', bam_out], stdout = FNULL, stderr = logfile)
        else:
            subprocess.call(['hisat2', '-p', str(args.cpus), FORWARD, REVERSE, SINGLE, '|', 'samtools', 'view', '-@', str(args.cpus), '-bS', '-', '|', 'samtools', 'sort', '-@', str(args.cpus), '-', bam_out], stdout = FNULL, stderr = logfile)
        RNAseqBAM = bam_out + '.bam'
    '''

    #now soft-mask the genome, so do that with RepeatModeler/RepeatMasker
    masked_genome = aug_species + '.softmasked.fa'
    if not os.path.isfile(masked_genome):
        lib.RepeatModelMask(args.input, args.cpus, args.out, masked_genome)
    else:
        lib.log.info("Soft-masked genome found, skipping repeat masking")

    #now need to run BRAKER1
    lib.log.info("Now launching BRAKER to train GeneMark and Augustus")
    species = '--species=' + aug_species
    genome = '--genome=' + os.path.abspath(masked_genome)
    bam = '--bam=' + os.path.abspath(args.rna_bam)
    with open('braker.log') as logfile:
        subprocess.call(['braker.pl', '--fungus', '--cores', str(args.cpus), '-gff3', '--softmasking', '1', genome, species, bam], stdout = logfile, stderr = logfile)




elif args.pipeline == 'fCEGMA':
    if not args.augustus_species:
        aug_species = args.species.replace(' ', '_').lower()
    else:
        aug_species = args.augustus_species
    if lib.CheckAugustusSpecies(aug_species):
        lib.log.error("%s as already been trained, choose a unique species name" % (aug_species))
        os._exit(1)

    #start by reformatting fasta headers
    sorted_input = os.path.join(args.out, 'genome.sorted.fa')
    lib.SortRenameHeaders(args.input, sorted_input)

    #first task is to soft-mask the genome, so do that with RepeatModeler/RepeatMasker
    masked_genome = aug_species + '.softmasked.fa'
    masked_genome = os.path.abspath(masked_genome)
    if not os.path.isfile(masked_genome):
        lib.RepeatModelMask(sorted_input, args.cpus, args.out, masked_genome)
    else:
        lib.log.info("Soft-masked genome found, skipping repeat masking")

    #now run GeneMark-ES to get models
    genemark = args.out + '.genemark.gff3'
    genemark = os.path.abspath(genemark)
    lib.RunGeneMarkES(masked_genome, args.cpus, args.out, genemark)

    #run GAG to get protein sequences
    lib.log.info("Prepping data using GAG")
    subprocess.call(['gag.py', '-f', masked_genome, '-g', genemark, '--fix_start_stop', '-o', 'genemark_gag'], stdout = FNULL, stderr = FNULL)
    os.rename(os.path.join('genemark_gag', 'genome.proteins.fasta'), os.path.join(args.out, 'genemark.proteins.fasta'))

    #filter using fCEGMA models
    lib.log.info("Now filtering best fCEGMA models for training Augustus")
    fCEGMA_in = os.path.join(args.out, 'genemark.proteins.fasta')
    fCEGMA_out = os.path.join(args.out, 'fCEGMA_hits.txt')
    lib.fCEGMA(fCEGMA_in, args.cpus, 1e-100, args.out, genemark, fCEGMA_out)

    #now run Augustus training based on training set
    fCEGMA_gff = os.path.join(args.out, 'training.gff3')
    fCEGMA_gff = os.path.abspath(fCEGMA_gff)
    lib.log.info("Now training Augustus with fCEGMA filtered dataset")
    if os.path.exists('autoAug'):
        shutil.rmtree('autoAug')
    species = '--species=' + aug_species
    genome = '--genome=' + masked_genome
    training = '--trainingset=' + fCEGMA_gff
    aug_log = os.path.join(args.out, 'augustus.log')
    with open(aug_log, 'w') as logfile:
        subprocess.call(['autoAug.pl', '--noutr', '--singleCPU', species, genome, training], stdout=logfile, stderr=logfile)



elif args.pipeline == 'no_train':
    if not args.augustus_species:
        lib.log.error("You must specifiy a valid species training set for Augustus, --augustus_species botrytis_cinerea")
        os._exit(1)
    if not args.genemark_mod or not os.path.isfile(args.genemark_mod):
        lib.log.error("You must specifiy a valid GeneMark mod file, e.g. --genemark_mod gmhmm.mod")
        os._exit(1)
    if not lib.CheckAugustusSpecies(args.augustus_species):
        lib.log.error("%s not found in Augustus/config/species folder" % (args.augustus_species))
        os._exit(1)

    #start by reformatting fasta headers
    sorted_input = os.path.join(args.out, 'genome.sorted.fa')
    lib.SortRenameHeaders(args.input, sorted_input)

    #first task is to soft-mask the genome, so do that with RepeatModeler/RepeatMasker
    masked_genome = aug_species + '.softmasked.fa'
    if not os.path.isfile(masked_genome):
        lib.RepeatModelMask(args.input, args.cpus, args.out, masked_genome)
    else:
        lib.log.info("Soft-masked genome found, skipping repeat masking")


elif args.pipeline == 'no_augustus_train':
    if not args.augustus_species:
        lib.log.error("You must specifiy a valid species training set for Augustus, --augustus_species botrytis_cinerea")
        os._exit(1)
    if not lib.CheckAugustusSpecies(args.augustus_species):
        lib.log.error("%s not found in Augustus/config/species folder")
        os._exit(1)

    #start by reformatting fasta headers
    sorted_input = os.path.join(args.out, 'genome.sorted.fa')
    lib.SortRenameHeaders(args.input, sorted_input)

    #first task is to soft-mask the genome, so do that with RepeatModeler/RepeatMasker
    masked_genome = aug_species + '.softmasked.fa'
    if not os.path.isfile(masked_genome):
        lib.RepeatModelMask(args.input, args.cpus, args.out, masked_genome)
    else:
        lib.log.info("Soft-masked genome found, skipping repeat masking")

elif args.pipeline == 'only_annotate_proteins':
    print "Taday!"




