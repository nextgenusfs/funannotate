#!/usr/bin/env python
# -*- coding: utf-8 -*-

from __future__ import (absolute_import, division,
                        print_function, unicode_literals)
import sys
import os
import importlib
import subprocess
import funannotate
from funannotate.__version__ import __version__

global package_name
package_name = 'funannotate'

default_help = """
Usage:       {:} <command> <arguments>
version:     {:}

Description: Funannotate is a genome prediction, annotation, and comparison pipeline.

Commands:
  clean       Find/remove small repetitive contigs
  sort        Sort by size and rename contig headers
  mask        Repeatmask genome assembly

  train       RNA-seq mediated training of Augustus/GeneMark
  predict     Run gene prediction pipeline
  fix         Fix annotation errors (generate new GenBank file)
  update      RNA-seq/PASA mediated gene model refinement
  remote      Partial functional annotation using remote servers
  iprscan     InterProScan5 search (Docker or local)
  annotate    Assign functional annotation to gene predictions
  compare     Compare funannotated genomes

  util        Format conversion and misc utilities
  setup       Setup/Install databases
  test        Download/Run funannotate installation tests
  check       Check Python, Perl, and External dependencies [--show-versions]
  species     list pre-trained Augustus species
  database    Manage databases
  outgroups   Manage outgroups for funannotate compare

Written by Jon Palmer (2016-2022) nextgenusfs@gmail.com with contributions by Jason Stajich jasonstajich.phd@gmail.com
        """.format(package_name, __version__)

cleanHelp = """
Usage:       {:} clean <arguments>
version:     {:}

Description: The script sorts contigs by size, starting with shortest contigs it uses minimap2
             to find contigs duplicated elsewhere, and then removes duplicated contigs.

Arguments:
  -i, --input    Multi-fasta genome file (Required)
  -o, --out      Cleaned multi-fasta output file (Required)
  -p, --pident   Percent identity of overlap. Default = 95
  -c, --cov      Percent coverage of overlap. Default = 95
  -m, --minlen   Minimum length of contig to keep. Default = 500
  --exhaustive   Test every contig. Default is to stop at N50 value.
  --cpus         Number of CPUs to use. Default: 2
  --tmpdir       TMP directory to use
  --debug        Debug the output
        """.format(package_name, __version__)

sortHelp = """
Usage:       {:} sort <arguments>
version:     {:}

Description: This script sorts the input contigs by size (longest->shortest) and then relabels
             the contigs with a simple name (e.g. scaffold_1).  Augustus can have problems with
             some complicated contig names.  Alternatively pass -s,--simplify in order
             to split fasta headers at first space.

Arguments:
  -i, --input    Multi-fasta genome file. (Required)
  -o, --out      Sorted by size and relabeled output file. (Required)
  -s, --simplify Try to simplify the FASTA headers, split at first space.
  -b, --base     Base name to relabel contigs. Default: scaffold
  -m, --minlen   Shorter contigs are discarded. Default: 0
        """.format(package_name, __version__)

maskHelp = """
Usage:       {:} mask <arguments>
version:     {:}

Description: This script is a wrapper for repeat masking. Default is to run very simple
             repeat masking with tantan. The script can also run RepeatMasker and/or
             RepeatModeler. It will generate a softmasked genome. Tantan is probably not
             sufficient for soft-masking an assembly, but with RepBase no longer being
             available RepeatMasker/Modeler may not be functional for many users.

Arguments:
  -i, --input                  Multi-FASTA genome file. (Required)
  -o, --out                    Output softmasked FASTA file. (Required)

Optional:
  -m, --method                 Method to use. Default: tantan [repeatmasker, repeatmodeler]
  -s, --repeatmasker_species   Species to use for RepeatMasker
  -l, --repeatmodeler_lib      Custom repeat database (FASTA format)
  --cpus                       Number of cpus to use. Default: 2
  --debug                      Keep intermediate files
             """.format(package_name, __version__)

trainHelp = """
Usage:       {:} train <arguments>
version:     {:}

Description: Script is a wrapper for de novo genome-guided transcriptome assembly using
             Trinity followed by PASA. Illumina and Long-read (nanopore/pacbio) RNA-seq
             is also supported. Dependencies are hisat2, Trinity, samtools, fasta,
             minimap2, PASA.

Required:
  -i, --input              Genome multi-fasta file
  -o, --out                Output folder name

RNA-seq input: at least one source is required, multiple may be combined.
  -l, --left               Left/Forward FASTQ Illumina reads (R1). Space-separated list ok
  -r, --right              Right/Reverse FASTQ Illumina reads (R2). Space-separated list ok
  -s, --single             Single ended FASTQ reads. Space-separated list ok
  --left_norm              Pre-normalized left FASTQ reads (R1)
  --right_norm             Pre-normalized right FASTQ reads (R2)
  --single_norm            Pre-normalized single-ended FASTQ reads
  --pacbio_isoseq          PacBio Iso-seq long-reads (single file)
  --nanopore_cdna          Nanopore cDNA long-reads (single file)
  --nanopore_mrna          Nanopore direct mRNA long-reads (single file)
  --trinity                Pre-computed Trinity transcripts (FASTA)

Optional:
  --stranded               If RNA-seq library stranded. Default: no [RF,FR,F,R,no]
  --jaccard_clip           Turn on jaccard clip for dense genomes [Recommended for fungi]
  --no_normalize_reads     Skip read Normalization
  --trimmer                Read trimming software. Default: trimmomatic [trimmomatic,fastp,none]
  --no_trimmomatic         Skip Quality Trimming of reads (deprecated: use --trimmer none)
  --memory                 RAM to use for Jellyfish/Trinity. Default: 50G
  -c, --coverage           Depth to normalize reads. Default: 50
  -m, --min_coverage       Min depth for normalizing reads. Default: 5
  --pasa_db                Database to use. Default: sqlite [mysql,sqlite]
  --pasa_alignment_overlap PASA --stringent_alignment_overlap. Default: 30.0
  --aligners               Aligners to use with PASA: Default: minimap2 blat [gmap]
  --pasa_min_pct_aligned   PASA --MIN_PERCENT_ALIGNED. Default: 90
  --pasa_min_avg_per_id    PASA --MIN_AVG_PER_ID. Default: 95
  --pasa_num_bp_splice     PASA --NUM_BP_PERFECT_SPLICE_BOUNDARY. Default: 3
  --max_intronlen          Maximum intron length. Default: 3000
  --species                Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
  --strain                 Strain name, e.g. CEA10
  --isolate                Isolate name, e.g. Af293
  --header_length          Maximum length of FASTA headers. Default: 16
  --cpus                   Number of CPUs to use. Default: 2
  --no-progress            Do not print progress to stdout for long sub jobs
  --stop_after_trinity     Stop pipeline after Trinity genome-guided assembly, before PASA

ENV Vars:  If not passed, will try to load from your $PATH.
  --PASAHOME
  --TRINITYHOME
           """.format(package_name, __version__)

predictHelp = """
Usage:       {:} predict <arguments>
version:     {:}

Description: Script takes genome multi-fasta file and a variety of inputs to do a comprehensive whole
             genome gene prediction.  Uses AUGUSTUS, GeneMark, Snap, GlimmerHMM, BUSCO, EVidence Modeler,
             tbl2asn, tRNAScan-SE, Exonerate, minimap2.
Required:
  -i, --input              Genome multi-FASTA file (softmasked repeats)
  -o, --out                Output folder name
  -s, --species            Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"

Protein evidence: if none given, defaults to $FUNANNOTATE_DB/uniprot_sprot.fasta
  --protein_evidence       Proteins to map to genome. Space-separated list ok
  --protein_alignments     Pre-computed protein alignments in GFF3 format
  --p2g_prefilter          Pre-filter hits software selection. Default: diamond [tblastn]
  --p2g_pident             Exonerate percent identity. Default: 80
  --p2g_diamond_db         Premade diamond genome database for protein2genome mapping

Transcript evidence: optional, multiple sources may be combined.
  --transcript_evidence    mRNA/ESTs to align to genome. Space-separated list ok
  --transcript_alignments  Pre-computed transcript alignments in GFF3 format
  --aligners               Transcript aligners to use. Default: minimap2 [minimap2,gmap,blat]
  --rna_bam                RNA-seq mapped to genome to train Augustus/GeneMark-ET
  --stringtie              StringTie GTF result

Gene model input: pre-computed models passed through to EVM.
  --pasa_gff               PASA generated gene models. filename:weight
  --maker_gff              MAKER2 GFF file. Parse results directly to EVM
  --other_gff              Annotation pass-through to EVM. filename:weight. Space-separated list ok
  --augustus_gff           Pre-computed AUGUSTUS GFF3 (must use --stopCodonExcludedFromCDS=False)
  --genemark_gtf           Pre-computed GeneMark GTF results
  --trnascan               Pre-computed tRNAscan-SE results

Ab initio training:
  -p, --parameters         Ab intio parameters JSON file to use for gene predictors
  -w, --weights            Ab-initio predictor and EVM weight. Example: augustus:2 or pasa:10
  --augustus_species       Augustus species config. Default: uses species name
                           (aliases: --trained_species, --species_parameters)
  --min_training_models    Minimum number of models to train Augustus. Default: 200
  --busco_fallback         Fall back to BUSCO training if too few models. Default: enabled
  --no_busco_fallback      Exit instead of falling back to BUSCO training
  --busco_seed_species     Augustus pre-trained species to start BUSCO. Default: anidulans
  --busco_db               BUSCO models. Default: dikarya. `funannotate outgroups --show_buscos`
  --augustus_BUSCO_timeout Per-job Augustus timeout during BUSCO, seconds. Default: 300
  --optimize_augustus      Run 'optimze_augustus.pl' to refine training (long runtime)
  -gm, --genemark_mode     GeneMark mode. Default: ES [ES,ET]
  --genemark_mod           GeneMark ini mod file
  --auto-skip-genemark     Skip GeneMark if assembly is fragmented (longest scaffold <50 kb)
  --soft_mask              Softmasked length threshold for GeneMark. Default: 2000

Optional:
  --isolate                Isolate name, e.g. Af293
  --strain                 Strain name, e.g. CEA10
  --name                   Locus tag name (assigned by NCBI?). Default: FUN_
  --numbering              Specify where gene numbering starts. Default: 1
  --organism               Fungal-specific options. Default: fungus [fungus,other]
  --ploidy                 Ploidy of assembly. Default: 1
  --table                  NCBI genetic code for nuclear CDS. Default: 1
  --mtable                 NCBI genetic code for mitochondrial CDS. Default: same as --table
  --min_intronlen          Minimum intron length. Default: 10
  --max_intronlen          Maximum intron length. Default: 3000
  --min_protlen            Minimum protein length. Default: 50
  --keep_no_stops          Keep gene models without valid stops
  --repeats2evm            Use repeats in EVM consensus model building
  --repeat_filter          Repetitive gene model filtering. Default: overlap blast [overlap,blast,none]
  --keep_evm               Keep existing EVM results (for rerunning pipeline)
  --evm-partition-interval Min length between genes to make a partition. Default: 1500
  --no-evm-partitions      Do not split contigs into partitions
  -t, --tbl2asn            Assembly parameters for tbl2asn. Default: "-l paired-ends"
  --SeqCenter              Sequencing facilty for NCBI tbl file. Default: CFMR
  --SeqAccession           Sequence accession number for NCBI tbl file. Default: 12345
  --force                  Annotated unmasked genome
  --header_length          Maximum length of FASTA headers. Default: 16
  --cpus                   Number of CPUs to use. Default: 2
  --no-progress            Do not print progress to stdout for long sub jobs
  --tmpdir                 Volume/location to write temporary files. Default: /tmp
  -d, --database           Path to funannotate database. Default: $FUNANNOTATE_DB

ENV Vars:  If not specified at runtime, will be loaded from your $PATH
  --EVM_HOME               (Perl version; optional if Rust evidence_modeler is in PATH)
  --AUGUSTUS_CONFIG_PATH
  --GENEMARK_PATH
           """.format(package_name, __version__)

updateHelp = """
Usage:       {:} update <arguments>
version:     {:}

Description: Script will run PASA mediated update of gene models. It can directly update
             the annotation from an NCBI downloaded GenBank file using RNA-seq data or can be
             used after funannotate predict to refine UTRs and gene model predictions. Kallisto
             is used to evidence filter most likely PASA gene models. Dependencies are
             hisat2, Trinity, samtools, fasta, minimap2, PASA, kallisto, bedtools.

Required: existing annotation, pass one of the following two options.
  -i, --input              Funannotate folder or Genome in GenBank format (.gbk,.gbff)
    or
  -f, --fasta              Genome in FASTA format
  -g, --gff                Annotation in GFF3 format
  --species                Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
                           Required with -f/-g; with -i it is parsed from the GenBank file

  -o, --out                Output folder name. Not required if -i is a funannotate folder

RNA-seq input: at least one source is required, multiple may be combined.
  -l, --left               Left/Forward FASTQ Illumina reads (R1). Space-separated list ok
  -r, --right              Right/Reverse FASTQ Illumina reads (R2). Space-separated list ok
  -s, --single             Single ended FASTQ reads. Space-separated list ok
  --left_norm              Pre-normalized left FASTQ reads (R1)
  --right_norm             Pre-normalized right FASTQ reads (R2)
  --single_norm            Pre-normalized single-ended FASTQ reads
  --pacbio_isoseq          PacBio Iso-seq long-reads (single file)
  --nanopore_cdna          Nanopore cDNA long-reads (single file)
  --nanopore_mrna          Nanopore direct mRNA long-reads (single file)
  --trinity                Pre-computed Trinity transcripts (FASTA)

Optional:
  --stranded               If RNA-seq library stranded. Default: no [RF,FR,F,R,no]
  --jaccard_clip           Turn on jaccard clip for dense genomes [Recommended for fungi]
  --no_normalize_reads     Skip read Normalization
  --no_trimmomatic         Skip Quality Trimming of reads
  --memory                 RAM to use for Jellyfish/Trinity. Default: 50G
  -c, --coverage           Depth to normalize reads. Default: 50
  -m, --min_coverage       Min depth for normalizing reads. Default: 5
  --pasa_config            PASA assembly config file, i.e. from previous PASA run
  --pasa_db                Database to use. Default: sqlite [mysql,sqlite]
  --pasa_alignment_overlap PASA --stringent_alignment_overlap. Default: 30.0
  --aligners               Aligners to use with PASA: Default: minimap2 blat [gmap]
  --pasa_min_pct_aligned   PASA --MIN_PERCENT_ALIGNED. Default: 90
  --pasa_min_avg_per_id    PASA --MIN_AVG_PER_ID. Default: 95
  --pasa_num_bp_splice     PASA --NUM_BP_PERFECT_SPLICE_BOUNDARY. Default: 3
  --max_intronlen          Maximum intron length. Default: 3000
  --min_protlen            Minimum protein length. Default: 50
  --alt_transcripts        Expression threshold (percent) to keep alt transcripts. Default: 0.10 [0-1]
  --cpus                   Number of CPUs to use. Default: 2
  --no-progress            Do not print progress to stdout for long sub jobs

NCBI/output options:
  --name                   Locus tag name (assigned by NCBI?). Default: use existing
  --strain                 Strain name
  --isolate                Isolate name
  --sbt                    NCBI Submission file. Default: bundled config/test.sbt
  --p2g                    NCBI p2g file (if updating NCBI annotation)
  -t, --tbl2asn            Assembly parameters for tbl2asn. Default: "-l paired-ends"
  --SeqCenter              Sequencing facilty for NCBI tbl file. Default: CFMR
  --SeqAccession           Sequence accession number for NCBI tbl file. Default: 12345
  --table                  NCBI genetic code (transl_table). Default: inherit from
                           predict_results/*.parameters.json, otherwise 1
  --mtable                 NCBI genetic code for mitochondrial CDS

Advanced (pass-through, skips the corresponding pipeline step):
  --pasa_gff               Pre-computed PASA GFF3, used in place of running PASA
  --kallisto               Pre-computed Kallisto abundance table (TSV)

ENV Vars:  If not passed, will try to load from your $PATH.
  --PASAHOME
  --TRINITYHOME
        """.format(package_name, __version__)

testHelp = """
Usage:       {:} test <arguments>
version:     {:}

Description: This is a script that runs several unit tests.  It will download data and run
             several different tests to determine if installion is functioning properly. If
             you cannot download from the machine funannotate is installed at - then download
             the 7 tar.gz files from https://osf.io/bj7v4/files/ and run script from directory

Arguments:
  -t, --tests    Test sets to run. [all,clean,mask,predict,busco,rna-seq,annotate,compare]
  --cpus         Number of cpus to use. Default: 2
  --debug        Keep output files
        """.format(package_name, __version__)

fixHelp = """
Usage:       {:} fix <arguments>
version:     {:}

Description: Script takes a GenBank genome annotation file and an NCBI tbl file to
             generate updated annotation. Script is used to fix problematic gene models
             after running funannotate predict or funannotate update.

Required:
  -i, --input    Annotated genome in GenBank format.
  -t, --tbl      NCBI tbl annotation file.

Optional:
  -d, --drop     Gene models to remove/drop from annotation. File with locus_tag 1 per line.
  -o, --out      Output folder
  --tbl2asn      Parameters for tbl2asn. Default: "-l paired-ends"
  --table        NCBI genetic code (transl_table). Default: inherit from sibling parameters.json
  --mtable       NCBI genetic code for mitochondrial CDS
        """.format(package_name, __version__)

remoteHelp = """
Usage:       {:} remote <arguments>
version:     {:}

Description: Script runs remote server functional annotation for Phobius and
             antiSMASH (fungi).  These searches are slow, if you can setup these services
             locally it will be much faster to do that.  PLEASE do not abuse services!

Required:
  -m, --methods       Which services to run, space separated [phobius,antismash,all]
  -e, --email         Email address to identify yourself to services.

  -i, --input         Funannotate input folder.
    or
  -g, --genbank       GenBank file (must be annotated).
  -o, --out           Output folder name.

Optional:
  -a, --antismash     antiSMASH server to use. Default: fungi [fungi,plants]
  --force             Force over-write of output folder
            """.format(package_name, __version__)

setupHelp = """
Usage:       {:} setup <arguments>
version:     {:}

Description: Script will download/format necessary databases for funannotate.

Options:
  -i, --install    Download format databases. Default: all
                     [merops,uniprot,dbCAN,pfam,repeats,go,mibig,
                      interpro,busco_outgroups,gene2product,busco]
  -b, --busco_db   Busco Databases to install. Default: dikarya [all,fungi,aves,etc]
  -d, --database   Path to funannotate database
  -u, --update     Check remote md5 and update if newer version found
  -l, --local      Use local json links instead of fetching from GitHub
  -f, --force      Force overwriting database
  -w, --wget       Use wget to download instead of python requests
  -l, --local      Use local resource JSON file instead of current on github
        """.format(package_name, __version__)

iprscanHelp = """
Usage:       {:} iprscan <arguments>
version:     {:}

Description: This script is a wrapper for running InterProScan5 using Docker or from a
             local installation. The script splits proteins into smaller chunks and then
             launches several interproscan.sh "processes". It then combines the results.

Arguments:
  -i, --input        Funannotate folder or FASTA protein file. (Required)
  -m, --method       Search method to use: [local, docker] (Required)
  -n, --num          Number of fasta files per chunk. Default: 1000
  -o, --out          Output XML InterProScan5 file
  --debug            Keep intermediate files
  --no-progress      Do not print progress to stdout for long sub jobs

Docker arguments:
  -c, --cpus         Number of CPUs (total). Default: 12
  --cpus_per_chunk   Number of cpus per Docker instance. Default: 4

Local arguments:
  --iprscan_path     Path to interproscan.sh. Default: which(interproscan.sh)
  -c, --cpus         Number of InterProScan instances to run
                     (configure cpu/thread control in interproscan.properties file)
        """.format(package_name, __version__)

annotateHelp = """
Usage:       {:} annotate <arguments>
version:     {:}

Description: Script functionally annotates the results from funannotate predict.  It pulls
             annotation from PFAM, InterPro, EggNog, UniProtKB, MEROPS, CAZyme, and GO ontology.

Required:
  -i, --input                 Folder from funannotate predict
    or
  --genbank                   Annotated genome in GenBank format
  -o, --out                   Output folder for results
    or
  --gff                       Genome GFF3 annotation file
  --fasta                     Genome in multi-fasta format
  -o, --out                   Output folder for results

             -s,--species is required unless the organism name can be parsed from --genbank
  -s, --species               Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"

Optional:
  --sbt                       NCBI submission template file. (Recommended) Default: bundled test.sbt
  -a, --annotations           Custom annotations (3 column TSV file)
  --eggnog                    Eggnog-mapper annotations file (if NOT installed)
  --iprscan                   InterProScan5 XML results
  --antismash                 antiSMASH secondary metabolism results (GBK file)
  --renumber_antismash        Renumber antiSMASH clusters to be contiguous
  --phobius                   Phobius pre-computed results (if phobius NOT installed)
  --signalp                   SignalP pre-computed results (-org euk -format short)
  --p2g                       NCBI p2g file (protein_id -> locus_tag) from previous annotation
  --fix                       Gene/Product names to fix (TSV: GeneID\\tName\\tProduct)
  --remove                    Gene/Product names to remove (TSV: Gene\\tProduct)
  --isolate                   Isolate name
  --strain                    Strain name
  --rename                    Rename GFF gene models with locus_tag from NCBI
  -m, --mito-pass-thru        Mitochondrial genome/contigs, append with :mcode
  --table                     NCBI genetic code (transl_table). Default: from parameters.json, else 1
  --mtable                    NCBI genetic code for mitochondrial CDS. Default: unset
                              (per-contig codes from -m,--mito-pass-thru take precedence)
  --allow_ec_without_genename Keep EC_number on genes with no gene name (tbl2asn rejects these)
  --busco_db                  BUSCO models. Default: dikarya
  -t, --tbl2asn               Additional parameters for tbl2asn. Default: "-l paired-ends"
  -d, --database              Path to funannotate database. Default: $FUNANNOTATE_DB
  --force                     Force over-write of output folder
  --cpus                      Number of CPUs to use. Default: 2
  --tmpdir                    Volume/location to write temporary files. Default: /tmp
  --header_length             Maximum length of FASTA headers. Default: 16
  --no-progress               Do not print progress to stdout for long sub jobs
         """.format(package_name, __version__)

compareHelp = """
Usage:       {:} compare <arguments>
version:     {:}

Description: Script does light-weight comparative genomics between funannotated genomes.  Output
             is graphs, phylogeny, CSV files, etc --> visualized in web-browser.

Required:
  -i, --input         List of funannotate genome folders or GBK files

Optional:
  -o, --out           Output folder name. Default: funannotate_compare
  -d, --database      Path to funannotate database. Default: $FUNANNOTATE_DB
  --cpus              Number of CPUs to use. Default: 2
  --run_dnds          Calculate dN/dS ratio on all orthologs. [estimate,full]
  --go_fdr            P-value for FDR GO-enrichment. Default: 0.05
  --heatmap_stdev     Cut-off for heatmap. Default: 1.0
  --num_orthos        Number of Single-copy orthologs to use for ML. Default: 500
  --bootstrap         Number of boostrap replicates to run with RAxML. Default: 100
  --outgroup          Name of species to use for ML outgroup. Default: no outgroup
  --eggnog_db         EggNog database. Default: fuNOG
  --proteinortho      Proteinortho POFF results. in TSV format.
  --ml_method         Maxmimum Likelihood method: Default: iqtree [raxml,iqtree]
  --ml_model          Substitution model for IQtree. Default: modelfinder
  --no-progress       Do not print progress to stdout for long sub jobs
          """.format(package_name, __version__)

outgroupHelp = """
Usage:       {:} outgroups <arguments>
version:     {:}

Description: Managing the outgroups folder for funannotate compare

Arguments:
  -i, --input            Proteome multi-fasta file. Required.
  -s, --species          Species name for adding a species. Required.
  -b, --busco_db         BUSCO db to use. Default. dikarya
  -c, --cpus             Number of CPUs to use for BUSCO search.
  -d, --database         Path to funannotate database. Default: $FUNANNOTATE_DB
          """.format(package_name, __version__)

utilHelp = """
Usage:       {:} util <arguments>
version:     {:}

Commands:
  stats              Generate assembly and annotation stats
  contrast           Compare annotations to reference (GFF3 or GBK annotations)
  tbl2gbk            Convert TBL format to GenBank format
  gbk2parts          Convert GBK file to individual components
  gff2prot           Convert GFF3 + FASTA files to protein FASTA
  gff2tbl            Convert GFF3 format to NCBI annotation table (tbl)
  bam2gff3           Convert BAM coord-sorted transcript alignments to GFF3
  prot2genome        Map proteins to genome generating GFF3 protein alignments
  stringtie2gff3     Convert GTF (stringTIE) to GFF3 format
  quarry2gff3        Convert CodingQuarry output to proper GFF3 format
  gff-rename         Sort GFF3 file and rename gene models
  predict-tRNA       Run tRNAscan-SE on a genome, output non-overlapping tRNA GFF3
          """.format(package_name, __version__)


statsHelp = """
Usage:       {:} util stats <arguments>
version:     {:}

Description: Generate JSON file with genome assembly and annotation stats.

Arguments:
  -f, --fasta              Genome FASTA file (Required)
  -o, --out                Output file (JSON format)
  -g, --gff3               Genome Annotation (GFF3 format)
  -t, --tbl                Genome Annotation (NCBI TBL format)
  --transcript_alignments  Transcript alignments (GFF3 format)
  --protein_alignments     Protein alignments (GFF3 format)
          """.format(package_name, __version__)


gff2tblHelp = """
Usage:       {:} util gff2tbl <arguments>
version:     {:}

Description: Convert GFF3 file into NCBI tbl format. Tbl output to stdout.

Arguments:
  -g, --gff3           Reference Annotation. GFF3 format
  -f, --fasta          Genome FASTA file.
  --table              NCBI genetic code (transl_table). Default: 1
          """.format(package_name, __version__)

prot2genomeHelp = """
Usage:       {:} util prot2genome <arguments>
version:     {:}

Description: Map proteins to genome using exonerate. Output is EVM compatible GFF3 file.

Arguments:   -g, --genome       Genome FASTA format (Required)
             -p, --proteins     Proteins FASTA format (Required)
             -o, --out          GFF3 output file (Required)
             -f, --filter       Pre-filtering method. Default: diamond [diamond,tblastn]
             -d, --filter_db    Premade diamond genome database for pre-filtering
             -t, --tblastn_out  Output to save tblastn results. Default: off
             --tblastn          Use existing tblastn results
             --exonerate_pident Exonerate percent identity. Default: 80
             --ploidy           Ploidy of assembly. Default: 1
             --maxintron        Max intron length. Default: 3000
             --contig_expand    Basepairs to expand alignments before running exonerate. Default: 3000
             --cpus             Number of cpus to use. Default: 2
             --tmpdir           Volume/location to write temporary files. Default: /tmp
             --logfile          Logfile output file. Default: funannotate-p2g.log
             --debug            Keep intermediate files
             --no-progress      Do not print progress to stdout for long sub jobs
           """.format(package_name, __version__)

gff2protHelp = """
Usage:       {:} util gff2prot <arguments>
version:     {:}

Description: Convert GFF3 file and genome FASTA to protein sequences. FASTA output to stdout.

Arguments:   -g, --gff3           Reference Annotation. GFF3 format
             -f, --fasta          Genome FASTA file.
             --no_stop            Dont print stop codons
             --table              NCBI genetic code (transl_table). Default: 1
            """.format(package_name, __version__)

gbk2partsHelp = """
Usage:       {:} util gbk2parts <arguments>
version:     {:}

Description: Convert GenBank file to its individual components (parts) tbl, protein
             FASTA, transcript FASTA, and contig/scaffold FASTA.

Arguments:   -g, --gbk          Input Genome in GenBank format
             -o, --output       Output basename
             --table            NCBI genetic code (transl_table). Default: from GBK, else 1
            """.format(package_name, __version__)

contrastHelp = """
Usage:       {:} util contrast <arguments>
version:     {:}

Description: Compare/constrast annotations to reference. Annotations in either GBK or GFF3 format.

Arguments:   -r, --reference            Reference Annotation. GFF3 or GBK format
             -f, --fasta                Genome FASTA. Required if GFF3 used
             -q, --query                Annotation query. GFF3 or GBK format
             -o, --output               Output basename
             -c, --calculate_pident     Measure protein percent identity between query and reference
             """.format(package_name, __version__)

tbl2gbkHelp = """
Usage:       {:} util tbl2gbk <arguments>
version:     {:}

Description: Convert NCBI TBL annotations + Genome FASTA to GenBank format.

Required:    -i, --tbl          Annotation in NCBI tbl format
             -f, --fasta        Genome FASTA file.
             -s, --species      Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
Optional:
             --isolate          Isolate name
             --strain           Strain name
             --sbt              NCBI Submission Template file
             -t, --tbl2asn      Assembly parameters for tbl2asn. Example: "-l paired-ends"
             -o, --output       Output basename
            """.format(package_name, __version__)

bam2gff3Help = """
Usage:       {:} util bam2gff3 <arguments>
version:     {:}

Description: Convert BAM coordsorted transcript alignments to GFF3 format.

Arguments:   -i, --bam           BAM file (coord-sorted)
             -o, --output        GFF3 output file
            """.format(package_name, __version__)

stringtieHelp = """
Usage:       {:} util stringtie2gff3 <arguments>
version:     {:}

Description: Convert StringTIE GTF format to GFF3 funannotate compatible format. Output
             to stdout.

Arguments:   -i, --input        GTF file from stringTIE
            """.format(package_name, __version__)

quarryHelp = """
Usage:       {:} util quarry2gff3 <arguments>
version:     {:}

Description: Convert CodingQuarry output GFF to proper GFF3 format. Output to stdout.

Arguments:   -i, --input        CodingQuarry output GFF file. (PredictedPass.gff3)
             -n, --numbering    Gene numbering starts at. Default: 1
            """.format(package_name, __version__)

gffrenameHelp = """
Usage:       {:} util gff-rename <arguments>
version:     {:}

Description: Sort GFF3 file by contigs and rename gene models.

Arguments:   -g, --gff3           Reference Annotation. GFF3 format
             -f, --fasta          Genome FASTA file.
             -o, --out            Output GFF3 file
             -l, --locus_tag      Locus tag to use. Default: FUN
             -n, --numbering      Start number for genes. Default: 1
             --table              NCBI genetic code (transl_table). Default: 1
            """.format(package_name, __version__)

predictTRNAHelp = """
Usage:       {:} util predict-tRNA <arguments>
version:     {:}

Description: Run tRNAscan-SE on a genome FASTA file and produce a
             non-overlapping tRNA GFF3 file, the same trnascan.no-overlaps.gff3
             style result produced in predict_misc by funannotate predict.

Arguments:   -i, --input          Genome FASTA file, should be soft-masked (Required)
             -o, --out            Output non-overlapping tRNA GFF3 file (Required)
             --genes              Existing gene models GFF3 to drop overlapping tRNAs
             --gaps               Assembly gaps GFF3/BED to drop overlapping tRNAs
             --trnascan           Pre-computed tRNAscan-SE tabular output
             --cpus               Number of CPUs for tRNAscan-SE. Default: 1
             --tmpdir             Directory for intermediate files. Default: directory of --out
             --keep_tmp           Keep intermediate/temporary files
             --logfile            Logfile output file
            """.format(package_name, __version__)


#  Add subcmds into dictionary
info = {'clean': {'cmd': 'clean', 'append': None, 'help': cleanHelp, 'dir': '.'},
        'sort': {'cmd': 'sort', 'append': None, 'help': sortHelp, 'dir': '.'},
        'mask': {'cmd': 'mask', 'append': None, 'help': maskHelp, 'dir': '.'},
        'train': {'cmd': 'train', 'append': None, 'help': trainHelp, 'dir': '.'},
        'predict': {'cmd': 'predict', 'append': None, 'help': predictHelp, 'dir': '.'},
        'update': {'cmd': 'update', 'append': None, 'help': updateHelp, 'dir': '.'},
        'fix': {'cmd': 'fix', 'append': None, 'help': fixHelp, 'dir': '.'},
        'remote': {'cmd': 'remote', 'append': None, 'help': remoteHelp, 'dir': '.'},
        'check': {'cmd': 'check', 'append': None, 'help': '', 'dir': '.'},
        'database': {'cmd': 'database', 'append': None, 'help': '', 'dir': '.'},
        'setup': {'cmd': 'setupDB', 'append': None, 'help': setupHelp, 'dir': '.'},
        'annotate': {'cmd': 'annotate', 'append': None, 'help': annotateHelp, 'dir': '.'},
        'outgroups': {'cmd': 'outgroups', 'append': None, 'help': outgroupHelp, 'dir': '.'},
        'compare': {'cmd': 'compare', 'append': None, 'help': compareHelp, 'dir': '.'},
        'iprscan': {'cmd': 'iprscan', 'append': None, 'help': iprscanHelp, 'dir': '.'},
        'species': {'cmd': 'species', 'append': None, 'help': '', 'dir': '.'},
        'util': {'cmd': None, 'append': None, 'help': utilHelp, 'dir': '.'},
        'stats': {'cmd': 'stats', 'append': None, 'help': statsHelp, 'dir': 'utilities'},
        'gff2tbl': {'cmd': 'gff2tbl', 'append': None, 'help': gff2tblHelp, 'dir': 'utilities'},
        'gff2prot': {'cmd': 'gff2prot', 'append': None, 'help': gff2protHelp, 'dir': 'utilities'},
        'gbk2parts': {'cmd': 'gbk2parts', 'append': None, 'help': gbk2partsHelp, 'dir': 'utilities'},
        'contrast': {'cmd': 'contrast', 'append': None, 'help': contrastHelp, 'dir': 'utilities'},
        'tbl2gbk': {'cmd': 'tbl2gbk', 'append': None, 'help': tbl2gbkHelp, 'dir': 'utilities'},
        'bam2gff3': {'cmd': 'bam2gff3', 'append': None, 'help': bam2gff3Help, 'dir': 'utilities'},
        'stringtie2gff3': {'cmd': 'stringtie2gff3', 'append': None, 'help': stringtieHelp, 'dir': 'utilities'},
        'quarry2gff3': {'cmd': 'quarry2gff3', 'append': None, 'help': quarryHelp, 'dir': 'utilities'},
        'prot2genome': {'cmd': 'funannotate-p2g.py', 'append': None, 'help': prot2genomeHelp,
                        'dir': 'aux_scripts', 'subprocess': True},
        'test': {'cmd': 'test', 'append': None, 'help': testHelp, 'dir': '.'},
        'gff-rename': {'cmd': 'gff_reformat', 'append': None, 'help': gffrenameHelp, 'dir': 'utilities'},
        'predict-tRNA': {'cmd': 'predict_trna', 'append': None, 'help': predictTRNAHelp, 'dir': 'utilities'}
        }

# Note, the first dict record would correspond to: package_name/example1.py script to import
# the append key is to pass a command silently to the script
# the help key is to reference the above help menu strings
# dir key is to deal with any nested folder structure
# if subprocess present, then run subprocess module so multithreading is functional

# main function: will display help menu for each subcommand and import it and run main() for that script


def main():
    # start here
    cmdName = None
    if len(sys.argv) < 2:
        print(default_help)
        sys.exit(1)
    elif sys.argv[1] == 'version' or sys.argv[1] == '--version' or sys.argv[1] == '-version' or sys.argv[1] == '-v':
        print("{:} v{:}".format(package_name, __version__))
        sys.exit(0)
    elif sys.argv[1] == 'util':
        try:
            cmdName = sys.argv[2]
        except IndexError:
            print(info['util']['help'])
            sys.exit(1)
        try:
            arguments = sys.argv[3:]
        except IndexError:
            print(info[cmdName]['help'])
            sys.exit(1)
    elif sys.argv[1] in info:
        cmdName = sys.argv[1]
        arguments = sys.argv[2:]
    else:
        print(default_help)
        sys.exit(1)
    if cmdName and cmdName in info:
        if len(arguments) > 0 or cmdName == 'check' or cmdName == 'database' or cmdName == 'species':
            if '-h' in arguments or '--help' in arguments:
                print(info[cmdName]['help'])
                sys.exit(0)
            if info[cmdName]['append']:
                arguments.append(info[cmdName]['append'])
            if 'subprocess' in info[cmdName]:
                # run as subprocess to be able to use multiprocessing
                cmd = [sys.executable, os.path.join(os.path.dirname(
                    funannotate.__file__), info[cmdName]['dir'], info[cmdName]['cmd'])]
                cmd += arguments
                subprocess.call(cmd)
            else:
                if info[cmdName]['dir'] != '.':
                    mod = importlib.import_module('{:}.{:}.{:}'.format(
                        package_name, info[cmdName]['dir'], info[cmdName]['cmd']))
                else:
                    mod = importlib.import_module(
                        '{:}.{:}'.format(package_name, info[cmdName]['cmd']))
                mod.main(arguments)
        else:
            print(info[cmdName]['help'])
            sys.exit(1)
    else:
        print(default_help)
        sys.exit(1)


if __name__ == "__main__":
    main()
