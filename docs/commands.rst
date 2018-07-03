
.. _commands:

Funannotate Commands
================
A description for all funannotate commands.

Funannotate wrapper script
-------------------------------------
Funannotate is a series of Python scripts that are launched from a Python wrapper script.  Each command has a help menu which you can print to the terminal by issuing the command without any arguments, i.e. ``funannotate`` yields the following.

.. code-block:: none
    
    $ funannotate
    Usage:       funannotate <command> <arguments>
    version:     1.4.0

    Description: Funannotate is a genome prediction, annotation, and comparison pipeline.
    
    Command:     clean          Find/remove small repetitive contigs
                 sort           Sort by size and rename contig headers
                 mask           Repeatmask genome assembly
             
                 train          RNA-seq mediated training of Augustus/GeneMark
                 predict        Run gene prediction pipeline
                 fix            Fix annotation errors (generate new GenBank file)
                 update         RNA-seq/PASA mediated gene model refinement
                 remote         Partial functional annotation using remote servers
                 iprscan        InterProScan5 search (Docker or local)
                 annotate       Assign functional annotation to gene predictions
                 compare        Compare funannotated genomes
             
                 setup          Setup/Install databases
                 util           Format conversion and misc utilities           
                 check          Check Python, Perl, and External dependencies
                 species        list pre-trained Augustus species
                 database       Manage databases             
                 outgroups      Manage outgroups for funannotate compare


    Preparing Genome for annotation
-------------------------------------
funannotate clean
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Script "cleans" an assembly by looking for duplicated contigs. The script first sorts the
contigs by size, then starting with the shortest contig it runs a "leave one out" alignment
using Mummer to determine if contig is duplicated elsewhere. This script is meant to be run
with a haploid genome, it has not been tested as a method to haplodize a polyploid assembly.

.. code-block:: none

    Usage:       funannotate clean <arguments>
    version:     1.4.0

    Description: The script sorts contigs by size, starting with shortest contigs it uses minimap2 or
                 Mummer to find contigs duplicated elsewhere, and then removes duplicated contigs.
    
    Arguments:   -i, --input    Multi-fasta genome file (Required)
                 -o, --out      Cleaned multi-fasta output file (Required)
                 -p, --pident   Percent identity of overlap. Default = 95
                 -c, --cov      Percent coverage of overlap. Default = 95
                 -m, --minlen   Minimum length of contig to keep. Default = 500
                 --exhaustive   Test every contig. Default is to stop at N50 value.
                 --method       Program to use to calculate overlaps. Default: minimap2 [minimap2,mummer]

funannotate sort
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Simple script to sort and rename a genome assembly. Often assemblers output contig/scaffold
names that are incompatible with NCBI submission rules. Use this script to rename and/or drop
scaffolds that are shorter than a minimum length.

.. code-block:: none

    Usage:       funannotate clean <arguments>
    version:     1.4.0

    Description: This script sorts the input contigs by size (longest->shortest) and then relabels
                 the contigs with a simple name (e.g. scaffold_1).  Augustus can have problems with
                 some complicated contig names.
    
    Arguments:   -i, --input    Multi-fasta genome file. (Required)
                 -o, --output   Sorted by size and relabeled output file. (Required)
                 -b, --base     Base name to relabel contigs. Default: scaffold
                 --minlen       Shorter contigs are discarded. Default: 0


funannotate species
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This function will output the current trained species in Augustus.

.. code-block:: none

    $ funannotate species
    --------------------------
    AUGUSTUS species options:
    --------------------------
    Conidiobolus_coronatus                      culex                                       parasteatoda                                
    E_coli_K12                                  debaryomyces_hansenii                       pchrysosporium                              
    Xipophorus_maculatus                        elegans                                     pea_aphid                                   
    adorsata                                    elephant_shark                              pfalciparum                                 
    aedes                                       encephalitozoon_cuniculi_GB                 phanerochaete_chrysosporium                 
    amphimedon                                  eremothecium_gossypii                       pichia_stipitis                             
    ancylostoma_ceylanicum                      fly                                         pneumocystis                                
    anidulans                                   fusarium                                    pseudogymnaoscus_verrucosos_1779            
    arabidopsis                                 fusarium_graminearum                        pseudogymnoascus_destructans                
    aspergillus_fumigatus                       galdieria                                   pseudogymnoascus_destructans_20631-21       
    aspergillus_nidulans                        generic                                     rhizopus_oryzae                             
    aspergillus_oryzae                          heliconius_melpomene1                       rhodnius                                    
    aspergillus_terreus                         histoplasma                                 rice                                        
    b_pseudomallei                              histoplasma_capsulatum                      s_aureus                                    
    bombus_impatiens1                           honeybee1                                   s_pneumoniae                                
    bombus_terrestris2                          human                                       saccharomyces                               
    botrytis_cinerea                            kluyveromyces_lactis                        saccharomyces_cerevisiae_S288C              
    brugia                                      laccaria_bicolor                            saccharomyces_cerevisiae_rm11-1a_1          
    c_elegans_trsk                              lamprey                                     schistosoma                                 
    cacao                                       leishmania_tarentolae                       schistosoma2                                
    caenorhabditis                              lodderomyces_elongisporus                   schizosaccharomyces_pombe                   
    camponotus_floridanus                       magnaporthe_grisea                          seahare                                     
    candida_albicans                            maize                                       sulfolobus_solfataricus                     
    candida_guilliermondii                      maize5                                      template_prokaryotic                        
    candida_tropicalis                          malassezia_vespertilionis                   tetrahymena                                 
    chaetomium_globosum                         nasonia                                     thermoanaerobacter_tengcongensis            
    chicken                                     neurospora                                  tomato                                      
    chlamy2011                                  neurospora_crassa                           toxoplasma                                  
    chlamydomonas                               ophidiomyces_44736-75                       tribolium2012                               
    chlorella                                   ophidiomyces_ophidicola                     trichinella                                 
    coccidioides_immitis                        ophidiomyces_ophiodiicola                   ustilago                                    
    coprinus                                    ophiostoma_novo-ulmi                        ustilago_maydis                             
    coprinus_cinereus                           paranannizziopsis_24878-07                  verticillium_albo_atrum1                    
    coyote_tobacco                              paranannizziopsis_26609-01                  verticillium_longisporum1                   
    cryptococcus                                paranannizziopsis_27236-01                  wheat                                       
    cryptococcus_neoformans_gattii              paranannizziopsis_44736-75                  yarrowia_lipolytica                         


funannotate mask
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Repetitive elements should be soft-masked from a genome assembly to help direct the ab-initio gene
predictors. This can be accomplished with the often used RepeatModeler/RepeatMasker programs. 
A wrapper for RepeatModeler/RepeatMasker is the :code:`funannotate mask` script. Note you can
use any other software to soft-mask your genome prior to running the gene prediction script.

.. code-block:: none

    Usage:       funannotate mask <arguments>
    version:     1.4.0-9711747

    Description: This script is a wrapper for RepeatMasker and RepeatModeler. It will generate
                 a softmasked genome. Dependencies are RepeatMasker and RepeatModeler. 
    
    Arguments:   -i, --input                    Multi-FASTA genome file. (Required)
                 -o, --out                      Output softmasked FASTA file. (Required)
                 -s, --repeatmasker_species     Species to use for RepeatMasker
                 -l, --repeatmodeler_lib        Custom repeat database (FASTA format)
                 --cpus                         Number of cpus to use. Default: 2
                 --debug                        Keep intermediate files
                 

Training Ab-initio Gene Predictors
-------------------------------------
funannotate train
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In order to use this script you will need RNA-seq data from the genome you are annotating, if
you don't have RNA-seq data then `funannotate predict` will train Augustus during runtime. This script
is a wrapper for genome-guided Trinity RNA-seq assembly followed by PASA assembly.  These methods
will generate the input data to `funannotate predict`, i.e. coord-sorted BAM alignments, trinity
transcripts, and high quality PASA GFF3 annotation. This script unfortunately has lots of dependencies
that include Hisat2, Trinity, Samtools, Fasta, GMAP, Blat, MySQL, PASA, and RapMap. The $PASAHOME
and $TRINITYHOME environmental variables need to be set or passed at runtime.

.. code-block:: none

    Usage:       funannotate train <arguments>
    version:     1.4.0

    Description: Script is a wrapper for Illumina RNA-seq data through genome-guided Trinity 
                 followed by PASA. Long-read (nanopore/pacbio) RNA-seq is also supported.
                 Dependencies are hisat2, Trinity, samtools, fasta, minimap2, blat, PASA.
    
    Required:  -i, --input              Genome multi-fasta file.
               -o, --out                Output folder name.
               -l, --left               Left/Forward FASTQ Illumina reads (R1)
               -r, --right              Right/Reverse FASTQ Illumina reads (R2)
               -s, --single             Single ended FASTQ reads

    Optional:  --stranded               If RNA-seq library stranded. [RF,FR,F,R,no]
               --left_norm              Normalized left FASTQ reads (R1)
               --right_norm             Normalized right FASTQ reads (R2)
               --single_norm            Normalized single-ended FASTQ reads
               --pacbio_isoseq          PacBio long-reads
               --nanopore_cdna          Nanopore cDNA long-reads
               --nanopore_mrna          Nanopore mRNA direct long-reads
               --trinity                Pre-computed Trinity transcripts (FASTA)
               --jaccard_clip           Turn on jaccard clip for dense genomes [Recommended for fungi]
               --no_normalize_reads     Skip read Normalization
               --no_trimmomatic         Skip Quality Trimming of reads
               --no_antisense_filter    Skip anti-sense filtering.
               --memory                 RAM to use for Jellyfish. Default: 50G
               -c, --coverage           Depth to normalize reads. Default: 50
               --pasa_db                Database to use. Default: sqlite [mysql,sqlite]
               --pasa_alignment_overlap PASA --stringent_alignment_overlap. Default: 30.0
               --max_intronlen          Maximum intron length. Default: 3000
               --species                Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
               --strain                 Strain name
               --isolate                Isolate name
               --cpus                   Number of CPUs to use. Default: 2

    ENV Vars:  If not passed, will try to load from your $PATH. 
               --PASAHOME
               --TRINITYHOME


Gene Prediction
-------------------------------------
funannotate predict
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This script is the "meat and potatoes" of funannotate. It will parse the data you provide
and choose the best method to train the ab-initio gene predictors Augustus and GeneMark. After
the predictors are trained, it runs Evidence Modeler to generate consensus gene models from
all of the data present. Finally, the GFF3 file is converted to NCBI GenBank format.

.. code-block:: none

    Usage:       funannotate predict <arguments>
    version:     1.4.0

    Description: Script takes genome multi-fasta file and a variety of inputs to do a comprehensive whole
                 genome gene prediction.  Uses AUGUSTUS, GeneMark, BUSCO, BRAKER, EVidence Modeler,
                 tbl2asn, tRNAScan-SE, Exonerate, minimap2.
    
    Required:  -i, --input              Genome multi-FASTA file (softmasked repeats).
               -o, --out                Output folder name.
               -s, --species            Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"

    Optional:  --isolate                Isolate name, e.g. Af293
               --strain                 Strain name, e.g. FGSCA4           
               --name                   Locus tag name (assigned by NCBI?). Default: FUN_
               --numbering              Specify where gene numbering starts. Default: 1
               --maker_gff              MAKER2 GFF file. Parse results directly to EVM.
               --pasa_gff               PASA generated gene models. filename:weight
               --other_gff              Annotation pass-through to EVM. filename:weight
               --rna_bam                RNA-seq mapped to genome to train Augustus/GeneMark-ET 
               --augustus_species       Augustus species config. Default: uses species name
               --genemark_mode          GeneMark mode. Default: ES [ES,ET]
               --genemark_mod           GeneMark ini mod file.
               --busco_seed_species     Augustus pre-trained species to start BUSCO. Default: anidulans
               --optimize_augustus      Run 'optimze_augustus.pl' to refine training (long runtime)
               --busco_db               BUSCO models. Default: dikarya. `funannotate outgroups --show_buscos`
               --organism               Fungal-specific options. Default: fungus. [fungus,other]
               --ploidy                 Ploidy of assembly. Default: 1
               -t, --tbl2asn            Assembly parameters for tbl2asn. Example: "-l paired-ends"
               -d, --database           Path to funannotate database. Default: $FUNANNOTATE_DB
           
               --protein_evidence       Proteins to map to genome (prot1.fa prot2.fa uniprot.fa). Default: uniprot.fa
               --protein_alignments     Pre-computed exonerate protein alignments (see docs for format)
               --transcript_evidence    mRNA/ESTs to align to genome (trans1.fa ests.fa trinity.fa). Default: none
               --transcript_alignments  Pre-computed transcript alignments in GFF3 format
               --augustus_gff           Pre-computed AUGUSTUS GFF3 results (must use --stopCodonExcludedFromCDS=False)
               --genemark_gtf           Pre-computed GeneMark GTF results
           
               --min_intronlen          Minimum intron length. Default: 10
               --max_intronlen          Maximum intron length. Default: 3000
               --soft_mask              Softmasked length threshold for GeneMark. Default: 5000
               --min_protlen            Minimum protein length. Default: 50
               --repeat_filter          Repetitive gene model filtering. Default: overlap blast [overlap,blast,none]
               --keep_no_stops          Keep gene models without valid stops.
               --SeqCenter              Sequencing facilty for NCBI tbl file. Default: CFMR
               --SeqAccession           Sequence accession number for NCBI tbl file. Default: 12345
               --force                  Annotated unmasked genome.
               --cpus                   Number of CPUs to use. Default: 2
             
    ENV Vars:  If not specified at runtime, will be loaded from your $PATH 
               --EVM_HOME
               --AUGUSTUS_CONFIG_PATH
               --GENEMARK_PATH
               --BAMTOOLS_PATH

funannotate fix
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
While funannotate predict does its best to generate gene models that will pass NCBI annotation
specs, occasionally gene models fall through the cracks (i.e. they are errors that the author
has not seen yet).  Gene models that generate submission errors are automatically flagged 
by funannotate predict and alerted to the user. The user must manually fix the .tbl annotation
file to fix these models. This script is a wrapper for archiving the previous genbank annotations
and generating a new set with the supplied .tbl annotation file.

.. code-block:: none

    Usage:       funannotate fix <arguments>
    version:     1.4.0

    Description: Script takes a GenBank genome annotation file and an NCBI tbl file to
                 generate updated annotation. Script is used to fix problematic gene models
                 after running funannotate predict.
    
    Required:    -i, --input    Annotated genome in GenBank format.
                 -t, --tbl      NCBI tbl annotation file.

    Optional:    -o, --out      Output folder
                 --tbl2asn      Parameters for tbl2asn. Default: "-l paired-ends"


funannotate update
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This script updates gene models from `funannotate predict` using RNA-seq data. The method relies
on RNA-seq --> Trinity --> PASA --> Kallisto. Using this script you can also update an NCBI
GenBank genome using RNA-seq data, i.e. you can update gene models on a pre-existing 
submission and the script will maintain proper annotation naming/updating in accordance with 
NCBI rules.

.. code-block:: none

    Usage:       funannotate update <arguments>
    version:     1.4.0

    Description: Script will run PASA mediated update of gene models. It can directly update
                 the annotation from an NCBI downloaded GenBank file using RNA-seq data or can be
                 used after funannotate predict to refine UTRs and gene model predictions. Kallisto
                 is used to evidence filter most likely PASA gene models. Dependencies are
                 hisat2, Trinity, samtools, fasta, minimap2, blat, PASA, kallisto, bedtools.
    
    Required:  -i, --input              Funannotate folder or Genome in GenBank format (.gbk,.gbff).
        or
               -f, --fasta              Genome in FASTA format
               -g, --gff                Annotation in GFF3 format
               --species                Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
           
    Optional:  -o, --out                Output folder name
               -l, --left               Left/Forward FASTQ Illumina reads (R1)
               -r, --right              Right/Reverse FASTQ Illumina reads (R2)
               -s, --single             Single ended FASTQ reads
               --stranded               If RNA-seq library stranded. [RF,FR,F,R,no]
               --left_norm              Normalized left FASTQ reads (R1)
               --right_norm             Normalized right FASTQ reads (R2)
               --single_norm            Normalized single-ended FASTQ reads
               --pacbio_isoseq          PacBio long-reads
               --nanopore_cdna          Nanopore cDNA long-reads
               --nanopore_mrna          Nanopore mRNA direct long-reads
               --trinity                Pre-computed Trinity transcripts (FASTA)
               --jaccard_clip           Turn on jaccard clip for dense genomes [Recommended for fungi]
               --no_antisense_filter    Skip anti-sense filtering
               --no_normalize_reads     Skip read Normalization
               --no_trimmomatic         Skip Quality Trimming of reads
               --memory                 RAM to use for Jellyfish. Default: 50G
               -c, --coverage           Depth to normalize reads. Default: 50
               --pasa_config            PASA assembly config file, i.e. from previous PASA run
               --pasa_db                Database to use. Default: sqlite [mysql,sqlite]
               --pasa_alignment_overlap PASA --stringent_alignment_overlap. Default: 30.0
               --max_intronlen          Maximum intron length. Default: 3000
               --min_protlen            Minimum protein length. Default: 50
               --alt_transcripts        Expression threshold (percent) to keep alt transcripts. Default: 0.1 [0-1]
               --p2g                    NCBI p2g file (if updating NCBI annotation)
               -t, --tbl2asn            Assembly parameters for tbl2asn. Example: "-l paired-ends"           
               --name                   Locus tag name (assigned by NCBI?). Default: use existing  
               --sbt                    NCBI Submission file        
               --species                Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
               --strain                 Strain name
               --isolate                Isolate name
               --SeqCenter              Sequencing facilty for NCBI tbl file. Default: CFMR
               --SeqAccession           Sequence accession number for NCBI tbl file. Default: 12345
               --cpus                   Number of CPUs to use. Default: 2
             
    ENV Vars:  If not passed, will try to load from your $PATH. 
               --PASAHOME
               --TRINITYHOME


Adding Functional Annotation
-------------------------------------
funannotate remote
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Some programs are Linux-only and not compatible on Mac OSX, to accomodate all users there are
a series of remote based searches that can be done from the command line. anitSMASH secondary metabolite
gene cluster prediction, Phobius, and InterProScan5 can be done from this interface. Note that
if you can install these tools locally, those searches will likely be much faster and thus preferred.

.. code-block:: none

    Usage:       funannotate remote <arguments>
    version:     1.4.0

    Description: Script runs remote server functional annotation for Phobius, InterProScan5, and
                 antiSMASH (fungi).  These searches are slow, if you can setup these services locally
                 it will be much faster to do that.  PLEASE do not abuse services!  
    
    Required:    -i, --input         Funannotate input folder.
              or
                 -g, --genbank       GenBank file (must be annotated).
                 -o, --out           Output folder name.
              and   
                 -m, --methods       Which services to run, space separated [phobius,antismash,interproscan,all]
                 -e, --email         Email address to identify yourself to services.
             
    Optional:    --force             Force query even if antiSMASH server looks busy


funannotate iprscan
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This script is a wrapper for a local InterProScan5 run or a local Docker-based IPR run.  The Docker build uses the blaxterlab/interproscan image. 

.. code-block:: none

    Usage:       funannotate iprscan <arguments>
    version:     1.4.0

    Description: This script is a wrapper for running InterProScan5 using Docker or from a 
                 local installation. The script splits proteins into smaller chunks and then
                 launches several interproscan.sh "processes". It then combines the results.
                 Note if you are on a large cluster, you probably don't want to use this script
                 as likely the "cluster" mode of InterProScan5 will be faster.
    
    Arguments:   -i, --input        Funannotate folder or FASTA protein file. (Required)
                 -m, --method       Search method to use: [local, docker] (Required)
                 -n, --num          Number of fasta files per chunk. Default: 1000
                 -o, --out          Output XML InterProScan5 file
                    
        Docker arguments:
                 -c, --cpus         Number of CPUs (total). Default: 12     
                 --cpus_per_chunk   Number of cpus per Docker instance. Default: 4
             
        Local arguments:
                 --iprscan_path     Full path to interproscan.sh (Required)
                 -c, --cpus         Number of InterProScan instances to run
                                    (configure cpu/thread control in interproscan.properties file) 


funannotate annotate
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This script is run after `funannotate predict` or `funannotate update` and assigns functional
annotation to the protein coding gene models. The best functional annotation is done when
InterProScan 5 is run on your protein prior to running this script.

.. code-block:: none

    Usage:       funannotate annotate <arguments>
    version:     1.4.0

    Description: Script functionally annotates the results from funannotate predict.  It pulls
                 annotation from PFAM, InterPro, EggNog, UniProtKB, MEROPS, CAZyme, and GO ontology.
    
    Required:    -i, --input        Folder from funannotate predict
              or
                 --genbank          Genome in GenBank format
                 -o, --out          Output folder for results
              or   
                 --gff              Genome GFF3 annotation file
                 --fasta            Genome in multi-fasta format
                 -s, --species      Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
                 -o, --out          Output folder for results

    Optional:    --sbt              NCBI submission template file. (Recommended)
                 -a, --annotations  Custom annotations (3 column tsv file)
                 --eggnog           Eggnog-mapper annotations file (if NOT installed)
                 --antismash        antiSMASH secondary metabolism results (GBK file from output)
                 --iprscan          InterProScan5 XML file
                 --phobius          Phobius pre-computed results (if phobius NOT installed)
                 --isolate          Isolate name
                 --strain           Strain name
                 --rename           Rename GFF gene models with locus_tag from NCBI.
                 --fix              Gene/Product names fixed (TSV: GeneID   Name    Product)
                 --remove           Gene/Product names to remove (TSV: Gene Product)
                 --busco_db         BUSCO models. Default: dikarya
                 -t, --tbl2asn      Additional parameters for tbl2asn. Example: "-l paired-ends"
                 -d, --database     Path to funannotate database. Default: $FUNANNOTATE_DB
                 --force            Force over-write of output folder
                 --cpus             Number of CPUs to use. Default: 2

    ENV Vars:  If not specified at runtime, will be loaded from your $PATH  
                 --AUGUSTUS_CONFIG_PATH


Comparative Genomics
-------------------------------------
funannotate compare
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This script takes "funannotate" genomes (output from multiple `funannotate annotate`) and runs
some comparative genomic operations. The script compares the annotation and generates graphs,
CSV files, GO enrichment, dN/dS ratios, orthology, etc --> the output is visualized HTML format
in a web browser.

.. code-block:: none

    Usage:       funannotate compare <arguments>
    version:     1.4.0

    Description: Script does light-weight comparative genomics between funannotated genomes.  Output
                 is graphs, phylogeny, CSV files, etc --> visualized in web-browser.  
    
    Required:    -i, --input         List of funannotate genome folders or GBK files

    Optional:    -o, --out           Output folder name. Default: funannotate_compare
                 -d, --database      Path to funannotate database. Default: $FUNANNOTATE_DB
                 --cpus              Number of CPUs to use. Default: 2
                 --run_dnds          Calculate dN/dS ratio on all orthologs. [estimate,full]
                 --go_fdr            P-value for FDR GO-enrichment. Default: 0.05
                 --heatmap_stdev     Cut-off for heatmap. Default: 1.0
                 --num_orthos        Number of Single-copy orthologs to use for ML. Default: 500
                 --bootstrap         Number of boostrap replicates to run with RAxML. Default: 100
                 --outgroup          Name of species to use for ML outgroup. Default: no outgroup
                 --proteinortho      ProteinOrtho5 POFF results.
                 --ml_method         Maxmimum Liklihood method: Default: raxml [raxml,iqtree]
     

Installation and Database Management
-------------------------------------
funannotate setup
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This command needs to be run to download required databases. It requires the user to specify
a location to save the database files.  This location can then be added to the ~/.bash_profile
so funannotate knows where to locate the database files. 

.. code-block:: none

    Usage:       funannotate setup <arguments>
    version:     1.4.0

    Description: Script will download/format necessary databases for funannotate. 
    
    Options:     -i, --install    Download format databases. Default: all
                                  [merops,uniprot,dbCAN,pfam,repeats,go,
                                   mibig,interpro,busco_outgroups,gene2product]
                 -d, --database   Path to funannotate databse
                 -u, --update     Check remote md5 and update if newer version found
                 -f, --force      Force overwriting database

funannotate database
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Simple script displays the currently installed databases.

.. code-block:: none

    $ funannotate database

    Funannotate Databases currently installed:

      Database          Type        Version      Date         Num_Records   Md5checksum                     
      pfam              hmmer3      31.0         2017-02            16712   3e47eec766a99b599cb24f28c4d393f8
      gene2product      text        1.5          2018-03-03         24964   5822b876cd0c91e7eb37dff77da52e14
      interpro          xml         67.0         2018-03-01         33707   72ea4eaad2d7007a51e1946ad085ca1f
      dbCAN             hmmer3      6.0          2017-09-12           585   3cb06f6f93c72a56c9fa12a6294b41d5
      busco_outgroups   outgroups   1.0          2018-04-14             7   6795b1d4545850a4226829c7ae8ef058
      merops            diamond     12.0         2017-10-04          4968   d923f0177c6d27c3d2886c705347adc0
      mibig             diamond     1.3          2018-03-26         24085   84b3cd16e0b3b074e4b7ee18c6aa31fd
      uniprot           diamond     2018_06      2018-06-20        557713   99662b979378958e02f53b0153f42972
      go                text        2018-04-13   2018-04-13         47128   f04d79c24b4a42226638779b17e75c80
      repeats           diamond     1.0          2018-05-02         11950   4e8cafc3eea47ec7ba505bb1e3465d21

    To update a database type:
        funannotate setup -i DBNAME -d /usr/local/share/funannotate --force

funannotate outgroups
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
This script is a helper function to manage and update outgroups for `funannotate compare`. Outgroup
species can be specified in `funannotate compare` to use as a reference for BUSCO-mediated
maximum likelihood phylogeny. This script allows the user to add a genome to the available outgroups
folder by running BUSCO and formatting it appropriately. 

.. code-block:: none

    Usage:       funannotate outgroups <arguments>
    version:     1.4.0

    Description: Managing the outgroups folder for funannotate compare
    
    Arguments:   -i, --input            Proteome multi-fasta file. Required. 
                 --species              Species name for adding a species. Required.
                 --busco_db             BUSCO db to use for --add. Default. dikarya
                 --cpus                 Number of CPUs to use for BUSCO search.
                 --show_buscos          List the busco_db options
                 --show_outgroups       List the installed outgroup species.
                 -d, --database         Path to funannotate database. Default: $FUNANNOTATE_DB


