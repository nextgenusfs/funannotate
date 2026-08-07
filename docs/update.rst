
.. _update:

Adding UTRs and refining predictions
================================
If you have RNA-seq data and would like to use the PASA-mediated "annotation comparison" to add UTRs and refine gene model predictions, this can be accomplished using the :code:`funannotate update` command. This script can also be run as a stand-alone to re-align RNA-seq data and/or update an existing GenBank genome.

If you have run :code:`funannotate train` and then :code:`funannotate predict`, this script will re-use those data and you can simply pass :code:`funannotate update -i folder --cpus 12`.  This will add the gene predictions to the SQL database and then walk through each gene comparing to existing PASA alignments, PASA will make some adjustments to the gene models. As recommended by PASA developers, this is run twice in :code:`funannotate update`.


Why is :code:`funannotate update` so slow??
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The default SQL database for PASA is set to use SQLite -- this is for compatibility.  However, the limitation is that SQLite database in PASA is single threaded due to SQLite database lock issue. Thus even if you pass multiple cpus to the script, it will run all of the PASA steps single threaded, which can take a long time depending on PASA alignments and genome size. If you `setup PASA to use MySQL <https://github.com/PASApipeline/PASApipeline/wiki/setting-up-pasa-mysql>`_, then the scripts can run PASA multi-threaded and :code:`funannotate update` will run much faster.


.. code-block:: none

	Usage:       funannotate update <arguments>
	version:     1.8.16

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
