
.. _utilities:

Utilities
================================
There are several scripts that maybe useful to users to convert between different formats, these scripts are housed in the :code:`funannotate util` submenu.


.. code-block:: none

    $ funannotate util

	Usage:       funannotate util <arguments>
	version:     1.8.16

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

Generate genome assembly stats
------------------------------
To generate genome assembly stats in a JSON file.

.. code-block:: none

    $ funannotate util stats

	Usage:       funannotate util stats <arguments>
	version:     1.8.16

	Description: Generate JSON file with genome assembly and annotation stats.

	Arguments:
          -f, --fasta              Genome FASTA file (Required)
          -o, --out                Output file (JSON format)
          -g, --gff3               Genome Annotation (GFF3 format)
          -t, --tbl                Genome Annotation (NCBI TBL format)
          --transcript_alignments  Transcript alignments (GFF3 format)
          --protein_alignments     Protein alignments (GFF3 format)

Comparing/contrast annotations to a reference
---------------------------------------
To compare/contrast genome annotations between different GFF3 or GBK files.

.. code-block:: none

    $ funannotate util contrast

	Usage:       funannotate util contrast <arguments>
	version:     1.8.16

	Description: Compare/constrast annotations to reference. Annotations in either GBK or GFF3 format.

	Arguments: -r, --reference            Reference Annotation. GFF3 or GBK format
                   -f, --fasta                Genome FASTA. Required if GFF3 used
                   -q, --query                Annotation query. GFF3 or GBK format
                   -o, --output               Output basename
                   -c, --calculate_pident     Measure protein percent identity between query and reference

Format Conversion
---------------------------------------

.. code-block:: none

    $ funannotate util tbl2gbk

	Usage:       funannotate util tbl2gbk <arguments>
	version:     1.8.16

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


.. code-block:: none

    $ funannotate util gbk2parts

	Usage:       funannotate util gbk2parts <arguments>
	version:     1.8.16

	Description: Convert GenBank file to its individual components (parts) tbl, protein
				 FASTA, transcript FASTA, and contig/scaffold FASTA.

	Arguments:   -g, --gbk          Input Genome in GenBank format
				       -o, --output       Output basename
				       --table            NCBI genetic code (transl_table). Default: from GBK, else 1


.. code-block:: none

    $ funannotate util gff2prot

	Usage:       funannotate util gff2prot <arguments>
	version:     1.8.16

	Description: Convert GFF3 file and genome FASTA to protein sequences. FASTA output to stdout.

	Arguments: -g, --gff3           Reference Annotation. GFF3 format
                   -f, --fasta          Genome FASTA file.
                   --no_stop            Dont print stop codons
                   --table              NCBI genetic code (transl_table). Default: 1

.. code-block:: none

    $ funannotate util gff2tbl

	Usage:       funannotate util gff2tbl <arguments>
	version:     1.8.16

	Description: Convert GFF3 file into NCBI tbl format. Tbl output to stdout.

	Arguments:
	  -g, --gff3           Reference Annotation. GFF3 format
	  -f, --fasta          Genome FASTA file.
	  --table              NCBI genetic code (transl_table). Default: 1


.. code-block:: none

    $ funannotate util bam2gff3

	Usage:       funannotate util bam2gff3 <arguments>
	version:     1.8.16

	Description: Convert BAM coordsorted transcript alignments to GFF3 format.

	Arguments: -i, --bam           BAM file (coord-sorted)
                   -o, --output        GFF3 output file


.. code-block:: none

    $ funannotate util protein2genome

	Usage:       funannotate util prot2genome <arguments>
	version:     1.8.16

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

.. code-block:: none

    $ funannotate util stringtie2gff3

	Usage:       funannotate util stringtie2gff3 <arguments>
	version:     1.8.16

	Description: Convert StringTIE GTF format to GFF3 funannotate compatible format. Output
				 to stdout.

	Arguments:   -i, --input        GTF file from stringTIE

.. code-block:: none

    $ funannotate util quarry2gff3

	Usage:       funannotate util quarry2gff3 <arguments>
	version:     1.8.16

	Description: Convert CodingQuarry output GFF to proper GFF3 format. Output to stdout.

	Arguments:   -i, --input        CodingQuarry output GFF file. (PredictedPass.gff3)
                     -n, --numbering    Gene numbering starts at. Default: 1

  .. code-block:: none

    $ funannotate util gff-rename

	Usage:       funannotate util gff-rename <arguments>
	version:     1.8.16

	Description: Sort GFF3 file by contigs and rename gene models.

	Arguments:   -g, --gff3           Reference Annotation. GFF3 format
                     -f, --fasta          Genome FASTA file.
                     -o, --out            Output GFF3 file
                     -l, --locus_tag      Locus tag to use. Default: FUN
                     -n, --numbering      Start number for genes. Default: 1
                     --table              NCBI genetic code (transl_table). Default: 1
