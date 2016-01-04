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
parser=argparse.ArgumentParser(prog='funannotate.py', usage="%(prog)s [options] -i genome.fasta",
    description='''Script that does it all.''',
    epilog="""Written by Jon Palmer (2015) nextgenusfs@gmail.com""",
    formatter_class = MyFormatter)
parser.add_argument('-i','--input', required=True, help='Genome in FASTA format')
parser.add_argument('-g','--gff', help='GFF3 annotation of genome (if available)')
parser.add_argument('-o','--out', required=True, help='Basename of output files')
parser.add_argument('-s','--species', required=True, help='Species name (e.g. "Aspergillus fumigatus") use quotes if there is a space')
parser.add_argument('-n','--name', help='Shortname for genes, perhaps assigned by NCBI, eg. VC83')
parser.add_argument('-e','--email', required=True, help='Email address for IPRSCAN server')
parser.add_argument('--pipeline', required=True, choices=['rnaseq', 'fCEGMA', 'no_train', 'no_augustus_train', 'only_annotate_proteins'], help='Method to employ, BRAKER1, GeneMark/fCEGMA, no_training=')
parser.add_argument('--augustus_species', help='Specify species for Augustus')
parser.add_argument('--genemark_mod', help='Use pre-existing Genemark training file (e.g. gmhmm.mod)')
parser.add_argument('--protein_evidence', default='uniprot.fa', help='Specify protein evidence for Maker')
parser.add_argument('--transcript_evidence', help='Specify transcript evidence for Maker')
parser.add_argument('--forward', help='RNAseq forward reads')
parser.add_argument('--reverse', help='RNAseq reverse reads')
parser.add_argument('--single', help='RNAseq single end reads')
parser.add_argument('--cpus', default=1, type=int, help='Number of CPUs to use')
parser.add_argument('--skip_annotation', action='store_true', help='Skip functional annotation')
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
lib.log.info("Operating system: %s, %i cores, %i GB RAM" % (sys.platform, multiprocessing.cpu_count(), lib.MemoryCheck()))

programs = ['hmmscan','blastp','blastn','gag.py','tbl2asn','runiprscan','gmes_petap.pl', 'BuildDatabase', 'RepeatModeler', 'RepeatMasker', 'genemark_gtf2gff3','autoAug.pl', 'maker']

#create temp folder
if not os.path.exists(args.out):
    os.makedirs(args.out)

#run some prelim checks on the input
if args.pipeline == 'rnaseq':
    if not args.augustus_species:
        aug_species = args.species.replace(' ', '_').lower()
    else:
        aug_species = args.augustus_species
    if lib.CheckAugustusSpecies(aug_species):
        lib.log.error("%s as already been trained, choose a unique species name" % (aug_species))
        os._exit(1)

    #check dependencies
    programs = ['hmmscan','blastp','blastn','gag.py','tbl2asn','runiprscan','gmes_petap.pl', 'BuildDatabase', 'RepeatModeler', 'RepeatMasker', 'genemark_gtf2gff3', 'augustus', 'hisat2', 'hisat2-build', 'samtools','braker1.pl', 'rmOutToGFF3.pl']
    lib.CheckDependencies(programs)

    #check input, will need some RNAseq reads to run this pipeline
    if not args.forward or not args.single:
        lib.log.error("You specified RNAseq as a pipeline, but did not provide RNAseq reads")
        os._exit(1)

    #start by reformatting fasta headers
    sorted_input = os.path.join(args.out, 'genome.sorted.fa')
    lib.SortRenameHeaders(args.input, sorted_input)

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

    #now soft-mask the genome, so do that with RepeatModeler/RepeatMasker
    masked_genome = aug_species + '.softmasked.fa'
    if not os.path.isfile(masked_genome):
        lib.RepeatModelMask(args.input, args.cpus, args.out, masked_genome)
    else:
        lib.log.info("Soft-masked genome found, skipping repeat masking")



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


#now we can go several directions 1) no RNAseq data -> then GeneMark-ES/fCEMGA route or 2) with RNAseq -> BRAKER1 route, or 3) use pre-existing Augustus training set and/or genemark mod file, or 4) just annotate proteins from genome?

if args.skip_annotation:
    lib.log.info("Skipping functional annotation, script now finished")
    os._exit(1)
'''
#now you are ready to initialize maker and run genome annotation


#run interpro scan, in background hopefully....
if not os.path.exists('iprscan'):
    os.makedirs('iprscan')
#keep track of number of times you launched RunIprScan
IPRcount = 0
lib.log.info("Starting RunIprScan and running in background")
p = subprocess.Popen(['runiprscan', '-i', args.input, '-m', args.email, '-o', 'iprscan'], stdout=subprocess.PIPE, stderr=subprocess.PIPE, preexec_fn=preexec_function)
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
'''



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
        Need to expand the archives with tar
        then need to concatentation the hmm profiles
            find fuNOG_hmm/ -name *.hmm -type f -maxdepth 1 -exec cat '{}' \; > fuNOG_4.5.hmm
        then hmmpress the resulting file
            hmmpress fuNOG_4.5.hmm
        then cleanup mess
            rm -R fuNOG_hmm/
            rm fuNOG.hmm.tar.gz
            rm fuNOG.annotations.tsv.gz
        #note EggNog DB is about  10 GB in size

PFAM:   Download
        ftp://ftp.ebi.ac.uk/pub/databases/Pfam//current_release/Pfam-A.hmm.gz
        Gunzip
        hmmpress Pfam-A.hmm
        rm Pfam-A.hmm.gz
        Download the mapping file

dbCAN:  Download: http://csbl.bmb.uga.edu/dbCAN/download/dbCAN-fam-HMMs.txt
                  http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt
        rename: mv dbCAN-fam-HMMs.txt dbCAN.hmm
                mv FamInfo.txt dbCAN.info.txt
        reformat names, gsed -i 's/.hmm$//g' dbCAN.hmm
        hmmpress dbCAN.hmm
        filtering: for fungi, use E-value < 1e-17 and coverage > 0.45

MEROPS: Downlowd the merops_scan.lib from MEROPS website - you need a log in
        need to reformat the headers to something useful, this will do it in bash
        sed 's/ - /#/g' merops_scan.lib | while read line; do set -- "$line"; IFS="#"; declare -a Array=($*); if [[ "${Array[0]}" == ">"* ]]; then echo ${Array[0]} ${Array[2]}; else echo $line; fi; done > merops_formatted.fa
        | sed 's/{.*unit: //g' | sed 's/}//g' | sed 's/-/:/g' > merops_formatted.fa
        then makeblast db:
        makeblastdb -in merops_formatted.fa -input_type fasta -dbtype prot -title MEROPS -parse_seqids -out MEROPS




SwissProt:
            Download the uniprot database - easier to parse the names
            ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/complete/uniprot_sprot.fasta.gz
            Gunzip
            makeblastdb -in uniprot_sprot.fasta -input_type fasta -dbtype prot -title uniprot -parse_seqids -out uniprot


'''
