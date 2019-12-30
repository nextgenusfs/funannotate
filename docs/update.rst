
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
	version:     1.7.2

	Description: Script will run PASA mediated update of gene models. It can directly update
	the annotation from an NCBI downloaded GenBank file using RNA-seq data or can be
	used after funannotate predict to refine UTRs and gene model predictions. Kallisto
	is used to evidence filter most likely PASA gene models. Dependencies are
	hisat2, Trinity, samtools, fasta, minimap2, PASA, kallisto, bedtools.
	
	Required:  
	  -i, --input              Funannotate folder or Genome in GenBank format (.gbk,.gbff).
		or
	  -f, --fasta              Genome in FASTA format
	  -g, --gff                Annotation in GFF3 format
	  --species                Species name, use quotes for binomial, e.g. "Aspergillus fumigatus"
		   
	Optional:  
	  -o, --out                Output folder name
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
	  --no_normalize_reads     Skip read Normalization
	  --no_trimmomatic         Skip Quality Trimming of reads
	  --memory                 RAM to use for Jellyfish. Default: 50G
	  -c, --coverage           Depth to normalize reads. Default: 50
	  -m, --min_coverage       Min depth for normalizing reads. Default: 5
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


