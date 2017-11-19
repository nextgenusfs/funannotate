#!/usr/bin/env python

#Wrapper script for Funannotate package.
import sys, os, subprocess, inspect
from natsort import natsorted
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
import lib.library as lib

def flatten(l):
    flatList = []
    for elem in l:
        # if an element of a list is a list
        # iterate over this list and add elements to flatList 
        if type(elem) == list:
            for e in elem:
                flatList.append(e)
        else:
            flatList.append(elem)
    return flatList

def fmtcols(mylist, cols):
    justify = []
    for i in range(0,cols):
        length = max(map(lambda x: len(x), mylist[i::cols]))
        length += 2
        ljust = map(lambda x: x.ljust(length), mylist[i::cols])
        justify.append(ljust)
    justify = flatten(justify)
    num_lines = len(mylist) / cols
    lines = (' '.join(justify[i::num_lines]) 
             for i in range(0,num_lines))
    return "\n".join(lines)

version = '0.8.0'

default_help = """
Usage:       funannotate <command> <arguments>
version:     %s

Description: Funannotate is a genome prediction, annotation, and comparison pipeline.
    
Command:     clean          Find/remove small repetitive contigs
             sort           Sort by size and rename contig headers (recommended)
             species        list pre-trained Augustus species
             
             train          RNA-seq mediated training of Augustus/GeneMark
             predict        Run gene prediction pipeline
             update         RNA-seq/PASA mediated gene model refinement
             remote         Partial functional annotation using remote servers
             annotate       Assign functional annotation to gene predictions
             compare        Compare funannotated genomes
             
             fix            Remove adapter/primer contamination from NCBI error report
             check          Check Python module versions installed
             
             setup          Setup/Install databases and check dependencies
             database       Manage databases
             outgroups      Manage outgroups for funannotate compare
             
Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % version

if len(sys.argv) > 1:
    if sys.argv[1] == 'clean':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: The script sorts contigs by size, starting with shortest contigs it uses Mummer 
             to find contigs duplicated elsewhere, and then removes duplicated contigs.
    
Arguments:   -i, --input    Multi-fasta genome file (Required)
             -o, --out      Cleaned multi-fasta output file (Required)
             -p, --pident   Percent identity of overlap. Default = 95
             -c, --cov      Percent coverage of overlap. Default = 95
             -m, --minlen   Minimum length of contig to keep. Default = 500
            
Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-contig_cleaner.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
    elif sys.argv[1] == 'sort':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: This script sorts the input contigs by size (longest->shortest) and then relabels
             the contigs with a simple name (e.g. scaffold_1).  Augustus can have problems with
             some complicated contig names.
    
Arguments:   -i, --input    Multi-fasta genome file. (Required)
             -o, --output   Sorted by size and relabeled output file. (Required)
             -b, --base     Base name to relabel contigs. Default: scaffold
             --minlen       Shorter contigs are discarded. Default: 0
            
Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-sort_rename.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
            
    elif sys.argv[1] == 'train':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script is a wrapper for genome-guided Trinity followed by PASA. Dependencies are
             Hisat2, Trinity, Samtools, Fasta, GMAP, Blat, MySQL, PASA. RapMap is optional but will
             speed up analysis if data is PE stranded.
    
Required:  -i, --input              Genome multi-fasta file.
           -o, --out                Output folder name.
           -l, --left               Left/Forward FASTQ Illumina reads (R1)
           -r, --right              Right/Reverse FASTQ Illumina reads (R2)
           -s, --single             Single ended FASTQ reads

Optional:  --stranded               If RNA-seq library stranded. [RF,FR,F,R,no]
           --left_norm              Normalized left FASTQ reads (R1)
           --right_norm             Normalized right FASTQ reads (R2)
           --single_norm            Normalized single-ended FASTQ reads
           --trinity                Pre-computed Trinity transcripts (FASTA)
           --jaccard_clip           Turn on jaccard clip for dense genomes [Recommended for fungi]
           --no_normalize_reads     Skip read Normalization
           --no_trimmomatic         Skip Quality Trimming of reads
           --no_antisense_filter    Skip anti-sense filtering.
           --memory                 RAM to use for Jellyfish. Default: 50G
           -c, --coverage           Depth to normalize reads. Default: 50
           --pasa_alignment_overlap PASA --stringent_alignment_overlap. Default: 30.0
           --max_intronlen          Maximum intron length. Default: 3000
           --species                Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
           --strain                 Strain name
           --isolate                Isolate name
           --cpus                   Number of CPUs to use. Default: 2
             
ENV Vars:  If not passed, will try to load from your $PATH. 
           --PASAHOME
           --TRINITYHOME
            
Written by Jon Palmer (2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-train.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)            
            
    elif sys.argv[1] == 'predict':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script takes genome multi-fasta file and a variety of inputs to do a comprehensive whole
             genome gene model prediction.  Uses AUGUSTUS, GeneMark, BUSCO, BRAKER1, EVidence Modeler,
             GAG, tbl2asn, tRNAScan-SE, RepeatModeler, RepeatMasker, Exonerate, GMAP
    
Required:  -i, --input            Genome multi-fasta file.
           -o, --out              Output folder name.
           -s, --species          Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
       or
           --masked_genome        Soft-masked genome (repeats lowercase)
           --repeatmasker_gff3    RepeatMasker derived GFF3 file
           -s, --species          Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"

Optional:  --isolate              Isolate name, e.g. Af293
           --strain               Strain name.           
           --name                 Locus tag name (assigned by NCBI?). Default: FUN_
           --maker_gff            MAKER2 GFF file. Parse results directly to EVM.
           --pasa_gff             PASA generated gene models. filename:weight
           --other_gff            Annotation pass-through to EVM. filename:weight
           --rna_bam              RNA-seq mapped to genome to train Augustus/GeneMark-ET
           --repeatmasker_species Taxonomy to use for RepeatMasker, will skip RepeatModeler.      
           --augustus_species     Augustus species config. Default: uses species name
           --genemark_mod         GeneMark ini mod file.
           --protein_evidence     Proteins to map to genome (prot1.fa,prot2.fa,uniprot.fa). Default: uniprot.fa
           --transcript_evidence  mRNA/ESTs to align to genome (trans1.fa,ests.fa,trinity.fa). Default: none
           --busco_seed_species   Augustus pre-trained species to start BUSCO. Default: anidulans
           --optimize_augustus    Run 'optimze_augustus.pl' to refine training (long runtime)
           --busco_db             BUSCO models. Default: dikarya. `funannotate outgroups --show_buscos`
           --organism             Fungal-specific options. Default: fungus. [fungus,other]
           --ploidy               Ploidy of assembly. Default: 1
           -t, --tbl2asn          Assembly parameters for tbl2asn. Example: "-l paired-ends"
           -d, --database         Path to funannotate database. Default: $FUNANNOTATE_DB
           
           --augustus_gff         Pre-computed AUGUSTUS GFF3 results (must use --stopCodonExcludedFromCDS=False)
           --genemark_gtf         Pre-computed GeneMark GTF results
           --exonerate_proteins   Pre-computed exonerate protein alignments (see docs for format)
           --gmap_gff             Pre-computed transcript alignments (GFF3 gmap output)
           --repeatmodeler_lib    Pre-computed RepeatModeler library (multi-fasta)
           --use_diamond          Use Diamond instead of tblastn for mapping protein evidence
           
           --min_intronlen        Minimum intron length. Default: 10
           --max_intronlen        Maximum intron length. Default: 3000
           --min_protlen          Minimum protein length. Default: 50
           --keep_no_stops        Keep gene models without valid stops.
           --cpus                 Number of CPUs to use. Default: 2
             
ENV Vars:  If not specified at runtime, will be loaded from your $PATH 
           --FUNANNOTATE_DB
           --EVM_HOME
           --AUGUSTUS_CONFIG_PATH
           --GENEMARK_PATH
           --BAMTOOLS_PATH
            
Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-predict.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
                  
    elif sys.argv[1] == 'update':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script will run PASA mediated update of gene models. It can directly update
             the annotation from an NCBI downloaded GenBank file using RNA-seq data or can be
             used after funannotate predict to refine UTRs and gene model predictions. Kallisto
             is used to evidence filter most likely PASA gene models. Dependencies are
             Hisat2, Trinity, Samtools, Fasta, GMAP, Blat, MySQL, PASA, Kallisto, BedTools, 
             GAG. RapMap is optional but will speed up analysis if data is PE stranded.
    
Required:  -i, --input              Funannotate folder or Genome in GenBank format (.gbk,.gbff).

Optional:  -o, --out                Output folder name.
           -l, --left               Left/Forward FASTQ Illumina reads (R1)
           -r, --right              Right/Reverse FASTQ Illumina reads (R2)
           -s, --single             Single ended FASTQ reads
           --sbt                    NCBI Submission file.
           --stranded               If RNA-seq library stranded. [RF,FR,F,R,no]
           --left_norm              Normalized left FASTQ reads (R1)
           --right_norm             Normalized right FASTQ reads (R2)
           --single_norm            Normalized single-ended FASTQ reads
           --trinity                Pre-computed Trinity transcripts (FASTA)
           --jaccard_clip           Turn on jaccard clip for dense genomes [Recommended for fungi]
           --pasa_gff               PASA/Transdecoder GFF
           --pasa_config            PASA assembly config file
           --kallisto               Kallisto abundance tsv table
           --no_antisense_filter    Skip anti-sense filtering.
           --no_normalize_reads     Skip read Normalization
           --no_trimmomatic         Skip Quality Trimming of reads
           --memory                 RAM to use for Jellyfish. Default: 50G
           -c, --coverage           Depth to normalize reads. Default: 50
           --pasa_alignment_overlap PASA --stringent_alignment_overlap. Default: 30.0
           --max_intronlen          Maximum intron length. Default: 3000
           --min_protlen            Minimum protein length. Default: 50
           --p2g                    NCBI p2g file (if updating NCBI annotation)
           -t, --tbl2asn            Assembly parameters for tbl2asn. Example: "-l paired-ends"           
           --name                   Locus tag name (assigned by NCBI?). Default: use existing          
           --species                Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
           --strain                 Strain name
           --isolate                Isolate name
           --cpus                   Number of CPUs to use. Default: 2
             
ENV Vars:  If not passed, will try to load from your $PATH. 
           --PASAHOME
           --TRINITYHOME
            
Written by Jon Palmer (2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-update.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)                        

    elif sys.argv[1] == 'annotate':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script functionally annotates the results from funannotate predict.  It pulls
             annotation from PFAM, InterPro, EggNog, UniProtKB, MEROPS, CAZyme, and GO ontology.
    
Required:    -i, --input        Folder from funannotate predict
          or
             --genbank          Genome in GenBank format
             -o, --out          Output folder for results
          or   
             --gff              Genome GFF3 annotation file
             --fasta            Genome in multi-fasta format
             --proteins         Genome proteins in multi-fasta format
             -s, --species      Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
             -o, --out          Output folder for results

Optional:    --sbt              NCBI submission template file. (Recommended)
             --eggnog           Eggnog-mapper annotations file.
             --antismash        antiSMASH secondary metabolism results, GBK file.
             --iprscan          InterProScan XML file
             --phobius          Phobius pre-computed results.
             --isolate          Isolate name, e.g. Af293
             --strain           Strain name
             --busco_db         BUSCO models. Default: dikarya
             -t, --tbl2asn      Additional parameters for tbl2asn. Example: "-l paired-ends"
             -d, --database     Path to funannotate database. Default: $FUNANNOTATE_DB
             --force            Force over-write of output folder
             --cpus             Number of CPUs to use. Default: 2

ENV Vars:  If not specified at runtime, will be loaded from your $PATH  
             --FUNANNOTATE_DB
             --AUGUSTUS_CONFIG_PATH
            
Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-functional.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)

    elif sys.argv[1] == 'remote':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script runs remote server functional annotation for Phobius, InterProScan5, and
             antiSMASH (fungi).  These searches are slow, if you can setup these services locally
             it will be much faster to do that.  PLEASE do not abuse services!  
    
Required:    -i, --input         Funannotate input folder.
          or
             -g, --genbank       GenBank file (must be annotated).
             -o, --out           Output folder name.
          and   
             -m, --methods       Which services to run. Choices: phobius, antismash, interproscan, all
             -e, --email         Email address to identify yourself to services.
             
Optional:    --force             Force query even if antiSMASH server looks busy

Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
       
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-remote.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)

    elif sys.argv[1] == 'compare':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script does light-weight comparative genomics between funannotated genomes.  Output
             is graphs, phylogeny, CSV files, etc --> visualized in web-browser.  
    
Required:    -i, --input         List of funannotate genome folders or GBK files

Optional:    -o, --out           Output folder name. Default: funannotate_compare
             -d, --database      Path to funannotate database. Default: $FUNANNOTATE_DB
             --cpus              Number of CPUs to use. Default: 2
             --run_dnds          Calculate dN/dS ratio on all orthologs. [estimate,full]
             --go_fdr            P-value for FDR GO-enrichment. Default: 0.05
             --heatmap_stdev     Cut-off for heatmap. Default: 1.0
             --num_orthos        Number of Single-copy orthologs to use for RAxML. Default: 500
             --bootstrap         Number of boostrap replicates to run with RAxML. Default: 100
             --outgroup          Name of species to use for RAxML outgroup. Default: no outgroup
             --proteinortho      ProteinOrtho5 POFF results.

Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
       
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-compare.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
    elif sys.argv[1] == 'fix':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script parses an NCBI FCSreport.txt, which identifies regions of the assembly
             that contain adapter/primer contamination.  These regions are then removed using
             GAG.  
    
Required:    -i, --input         funannotate output folder
             -e, --error_report  NCBI FCSreport.txt
             --sbt               NCBI template submission file 

Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
       
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-fix.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
    elif sys.argv[1] == 'species':
        try:
            AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
        except KeyError:
            print("Error: Augustus is not properly configured. Please review installation instructions")
            sys.exit(1)
        #get the possible species from augustus
        augustus_list = []
        for i in os.listdir(os.path.join(AUGUSTUS, 'species')):
            if not i.startswith('.'):
                augustus_list.append(i)
        augustus_list = set(augustus_list)
        d = flatten(natsorted(augustus_list))
        print '--------------------------'
        print 'AUGUSTUS species options:'
        print '--------------------------'
        print lib.list_columns(d, cols=3)
        sys.exit(1)
        
    elif sys.argv[1] == 'check':
        cmd = os.path.join(script_path, 'util', 'check_modules.py')
        subprocess.call(cmd) 
        sys.exit(1)
    elif sys.argv[1] == 'setup':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script will download/format necessary databases for funannotate. 
    
Options:     -i, --install    Download format databases. Default: all
                              [merops,uniprot,dbCAN,pfam,repeats,go,
                               mibig,interpro,busco_outgroups,curated_names]
             -d, --database   Path to funannotate databse
             --force          Force overwriting database

Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)     
        arguments = sys.argv[2:]
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-setup.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)

    elif sys.argv[1] == 'database':
        #setup funannotate DB path
        try:
            FUNDB = os.environ["FUNANNOTATE_DB"]
        except KeyError:
            print('$FUNANNOTATE_DB not found, run funannotate setup and export ENV variable')
            sys.exit(1)
        dbfile = os.path.join(FUNDB, 'funannotate-db-info.txt')
        db_list = ['Database', 'Type', 'Version', 'Date', 'Num_Records']
        if not os.path.isfile(dbfile):
            print('Database is not properly configured, re-run funannotate setup')
            sys.exit(1)
        with open(dbfile, 'rU') as infile:
            for line in infile:
                line = line.rstrip()
                cols = line.split('\t')
                del cols[2]
                db_list.append(cols)

        d = flatten(db_list)
        db_print = fmtcols(d, 5)
        msg="""
--------------------------------------------------------------
Funannotate Databases currently installed:
--------------------------------------------------------------"""
        print msg    
        print db_print
        print('--------------------------------------------------------------')
        print('\nTo update a database type:\n\tfunannotate setup -i DBNAME -d {:} --force\n'.format(FUNDB))
        sys.exit(1)
    elif sys.argv[1] == 'outgroups':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Managing the outgroups folder for funannotate compare
    
Arguments:   -i, --input            Proteome multi-fasta file. Required. 
             --species              Species name for adding a species. Required.
             --busco_db             BUSCO db to use for --add. Default. dikarya
             --cpus                 Number of CPUs to use for BUSCO search.
             --show_buscos          List the busco_db options
             --show_outgroups       List the installed outgroup species.
               
Written by Jon Palmer (2016-2017) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
        
        arguments = sys.argv[2:]
        if '--show_outgroups' in arguments:
            files = [f for f in os.listdir(os.path.join(script_path, 'DB', 'outgroups'))]
            files = [ x.replace('_buscos.fa', '') for x in files ]
            files = [ x for x in files if not x.startswith('.') ]
            print "-----------------------------"
            print "BUSCO Outgroups:"
            print "-----------------------------"
            print lib.list_columns(files, cols=3)
            sys.exit(1)
        elif  '--show_buscos' in arguments:
            print "-----------------------------"
            print "BUSCO DB tree: (# of models)"
            print "-----------------------------"
            print lib.buscoTree
            sys.exit(1)
        elif len(arguments) > 1:
            cmd = os.path.join(script_path, 'util', 'add2outgroups.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            sys.exit(1)
    elif sys.argv[1] == 'version':
        print "funannotate v%s" % version
    else:
        print "%s option not recognized" % sys.argv[1]
        print default_help
        sys.exit(1)  
    
else:
    print default_help
        