
.. _utilities:

Utilities
================================
There are several scripts that maybe useful to users to convert between different formats, these scripts are housed in the :code:`funannotate util` submenu.


.. code-block:: none

    $ funannotate util

    Usage:       funannotate util <arguments>
    version:     1.4.0
    
    Commands:    compare        Compare annotations to reference (GFF3 or GBK annotations)
                 tbl2gbk        Convert TBL format to GenBank format
                 gbk2parts      Convert GBK file to individual components
                 gff2proteins   Convert GFF3 + FASTA files to protein FASTA
                 gff2tbl        Convert GFF3 format to NCBI annotation table (tbl)
                 bam2gff3       Convert BAM coord-sorted transcript alignments to GFF3 
                 

Comparing annotations to a reference
---------------------------------------
To compare genome annotations between different GFF3 or GBK files.

.. code-block:: none

    $ funannotate util compare

    Usage:       funannotate util <arguments>
    version:     1.4.0

    Description: Compare annotations to reference. Annotations in either GBK or GFF3 format.
    
    Arguments:   -r, --reference            Reference Annotation. GFF3 or GBK format
                 -f, --fasta                Genome FASTA. Required if GFF3 used
                 -q, --query                Annotation query. GFF3 or GBK format
                 -o, --output               Output basename
                 -c, --calculate_pident     Measure protein percent identity between query and reference
                 
                 
Format Conversion
---------------------------------------

.. code-block:: none

    $ funannotate util tbl2gbk

    Usage:       funannotate util <arguments>
    version:     1.4.0

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

    Usage:       funannotate util <arguments>
    version:     1.4.0

    Description: Convert GenBank file to its individual components (parts) tbl, protein
                 FASTA, transcript FASTA, and contig/scaffold FASTA.
    
    Arguments:   -g, --gbk          Input Genome in GenBank format
                 -o, --output       Output basename

.. code-block:: none

    $ funannotate util gff2proteins

    Usage:       funannotate util <arguments>
    version:     1.4.0

    Description: Convert GFF3 file into NCBI tbl format. FASTA output to stdout.
    
    Arguments:   -g, --gff3           Reference Annotation. GFF3 format
                 -f, --fasta          Genome FASTA file.
                 --no_stop            Dont print stop codons

.. code-block:: none

    $ funannotate util gff2tbl

    Usage:       funannotate util <arguments>
    version:     1.4.0

    Description: Convert GFF3 file into NCBI tbl format. Tbl output to stdout.
    
    Arguments:   -g, --gff3           Reference Annotation. GFF3 format
                 -f, --fasta          Genome FASTA file.


.. code-block:: none

    $ funannotate util bam2gff3

    Usage:       funannotate util <arguments>
    version:     1.4.0

    Description: Convert BAM coordsorted transcript alignments to GFF3 format.
    
    Arguments:   -i, --bam           BAM file (coord-sorted)
                 -o, --output        GFF3 output file
