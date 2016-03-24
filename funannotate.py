#!/usr/bin/env python

#Wrapper script for Funannotate package.

import sys, os, subprocess, inspect
from natsort import natsorted
script_path = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))

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

version = '0.1.2'

default_help = """
Usage:       funannotate <command> <arguments>
version:     %s

Description: Funannotate is a genome prediction, annotation, and comparison pipeline.
    
Command:     clean          Find/remove small repetitive contigs
             sort           Sort by size and rename contig headers (recommended)
             species        list pre-trained Augustus species
             
             predict        Run gene prediction pipeline
             annotate       Assign functional annotation to gene predictions
             compare        Compare funannotated genomes
             
             fix            Remove adapter/primer contamination from NCBI error report
             
Written by Jon Palmer (2016) nextgenusfs@gmail.com
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
            
Written by Jon Palmer (2016) nextgenusfs@gmail.com
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
            os._exit(1)
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
            
Written by Jon Palmer (2016) nextgenusfs@gmail.com
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
            os._exit(1)
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
             
Optional:  --isolate              Strain isolate, e.g. Af293            
           --name                 Locus tag name (assigned by NCBI?). Default: FUN_
           --pasa_gff             PASA generated gene models
           --rna_bam              RNA-seq mapped to genome to train Augustus/GeneMark-ET       
           --augustus_species     Augustus species config. Default: uses species name
           --genemark_mod         GeneMark ini mod file.
           --protein_evidence     Proteins to map to genome (prot1.fa,prot2.fa,uniprot.fa). Default: uniprot.fa
           --transcript_evidence  mRNA/ESTs to align to genome (trans1.fa,ests.fa,trinity.fa). Default: none
           --busco_seed_species   Augustus pre-trained species to start BUSCO. Default: generic  
           
           --augustus_gff         Pre-computed AUGUSTUS GFF3 results
           --genemark_gtf         Pre-computed GeneMark GTF results
           --exonerate_proteins   Pre-computed exonerate protein alignments (see docs for format)
           --gmap_gff             Pre-computed transcript alignments (GFF3 gmap output)
           --repeatmodeler_lib    Pre-computed RepeatModeler library (multi-fasta)

           --min_intronlen        Minimum intron length. Default: 10
           --max_intronlen        Maximum intron length. Default: 3000
           --min_protlen          Minimum protein length. Default: 50
           --cpus                 Number of CPUs to use. Default: 2
             
ENV Vars:  By default loaded from your $PATH, however you can specify at run-time if not in PATH  
           --EVM_HOME
           --AUGUSTUS_CONFIG_PATH
           --GENEMARK_PATH
           --BAMTOOLS_PATH
            
Written by Jon Palmer (2016) nextgenusfs@gmail.com
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
            os._exit(1)
    elif sys.argv[1] == 'annotate':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script functionally annotates the results from funannotate predict.  It pulls
             annotation from PFAM, InterPro, EggNog, UniProtKB, MEROPS, CAZyme, and GO ontology.
    
Required:    -i, --input        Folder from funannotate predict
             -e, --email        Valid email address for InterProScan server identification.
          or
             --genbank          Genome in GenBank format
             -e, --email        Valid email address for InterProScan server identification.
             -o, --out          Output folder for results
          or   
             --gff              Genome GFF3 annotation file
             --fasta            Genome in multi-fasta format
             --proteins         Proteins in multi-fasta format
             -e, --email        Valid email address for InterProScan server identification.
             -o, --out          Output folder for results

Optional:    --sbt              NCBI submission template file. (Recommended)
             --antismash        antiSMASH secondary metabolism results, GBK file.
             --iprscan          Folder of pre-computed InterProScan results (1 xml file per protein)
             -s, --species      Species name, use quotes for binomial, e.g. "Aspergillus fumigatus" 
             --isolate          Strain isolate, e.g. Af293  
             --skip_iprscan     Do not run InterProScan remote search.
             --force            Force over-write of output folder
             --cpus             Number of CPUs to use. Default: 2
            
Written by Jon Palmer (2016) nextgenusfs@gmail.com
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
            os._exit(1)
    elif sys.argv[1] == 'compare':
        help = """
Usage:       funannotate %s <arguments>
version:     %s

Description: Script does light-weight comparative genomics between funannotated genomes.  Output
             is graphs, CSV files, etc --> visualized in web-browser.  
    
Required:    -i, --input         List of funannotate genome folders

Optional:    -o, --out           Output folder name. Default: funannotate_compare
             --cpus              Number of CPUs to use. Default: 2
             --go_fdr            P-value for FDR GO-enrichment. Default: 0.05
             --heatmap_stdev     Cut-off for heatmap. Default: 1.0
             --num_orthos        Number of Single-copy orthologs to use for RAxML. Default: 500
             --bootstrap         Number of boostrap replicates to run with RAxML. Default: 100
             --outgroup          Name of species to use for RAxML outgroup. Default: no outgroup
             --show_outgroups    Show a list of pre-computed genomes to use as outgroups

Written by Jon Palmer (2016) nextgenusfs@gmail.com
        """ % (sys.argv[1], version)
       
        arguments = sys.argv[2:]
        if '--show_outgroups' in arguments:
            files = [f for f in os.listdir(os.path.join(script_path, 'DB', 'outgroups'))]
            files = [ x.replace('_buscos.fa', '') for x in files ]
            files = [ x for x in files if not x.startswith('.') ]
            print natsorted(files)
            os._exit(1)
        if len(arguments) > 1:
            cmd = os.path.join(script_path, 'bin', 'funannotate-compare.py')
            arguments.insert(0, cmd)
            exe = sys.executable
            arguments.insert(0, exe)
            subprocess.call(arguments)
        else:
            print help
            os._exit(1)
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

Written by Jon Palmer (2016) nextgenusfs@gmail.com
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
            os._exit(1)
    elif sys.argv[1] == 'species':
        try:
            AUGUSTUS = os.environ["AUGUSTUS_CONFIG_PATH"]
        except KeyError:
            print("Error: Augustus is not properly configured. Please review installation instructions")
            os._exit(1)
        #get the possible species from augustus
        augustus_list = []
        for i in os.listdir(os.path.join(AUGUSTUS, 'species')):
            if not i.startswith('.'):
                augustus_list.append(i)
        augustus_list = set(augustus_list)
        d = flatten(natsorted(augustus_list))
        print fmtcols(d, 3)
        os._exit(1)
 
    elif sys.argv[1] == 'version':
        print "funannotate v.%s" % version
    else:
        print "%s option not recognized" % sys.argv[1]
        print default_help
        os._exit(1)
    
    
else:
    print default_help
        