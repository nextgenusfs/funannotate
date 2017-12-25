
.. _dependencies:

Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Funannotate has a lot of dependencies.  However, it also comes with a few tools to help you get everything installed.  The first is that of :code:`funannotate check`.  You'll see in the output below that the :code:`fasta` tool is missing, which is Bill Pearsons :code:`fasta36` a dependency of the PASA pipeline.  Also the :code:`$PASAHOME`` and :code:`$TRINITYHOME`` variables are not set, that is because on this particular machine they are not installed, i.e. funannotate will alert you at runtime if it is missing a dependency.

.. code-block:: none
    
    $funannotate check --show-versions
    -------------------------------------------------------
    Checking dependencies for funannotate v1.0.0
    -------------------------------------------------------
    You are running Python v 2.7.14. Now checking python packages...
    biopython: 1.68
    goatools: 0.7.11
    matplotlib: 2.1.0
    natsort: 5.0.2
    numpy: 1.13.3
    pandas: 0.20.3
    psutil: 5.3.1
    scikit-learn: 0.19.0
    scipy: 0.19.1
    seaborn: 0.8.1
    All 10 python packages installed


    You are running Perl v 5.018002. Now checking perl modules...
    Bio::Perl: 1.006923
    Carp: 1.38
    Clone: 0.36
    DBI: 1.631
    Data::Dumper: 2.154
    File::Basename: 2.84
    Getopt::Long: 2.48
    Hash::Merge: 0.200
    Logger::Simple: 2.0
    POSIX: 1.32
    Parallel::ForkManager: 1.17
    Pod::Usage: 1.68
    Scalar::Util::Numeric: 0.40
    Storable: 2.41
    Text::Soundex: 3.04
    Thread::Queue: 3.07
    Tie::File: 0.99
    YAML: 1.15
    threads: 2.02
    threads::shared: 1.48
    All 20 Perl modules installed


    Checking external dependencies...
    RepeatMasker: RepeatMasker 4.0.7
    RepeatModeler: RepeatModeler 1.0.8
    Trinity: 2.4.0
    augustus: 3.2.1
    bamtools: bamtools 2.4.1
    bedtools: bedtools v2.26.0
    blat: BLAT v36
    braker.pl: braker.pl 
    diamond: diamond 0.9.13
    emapper.py: emapper-1.0.3
    ete3: 3.1.1
    exonerate: exonerate 2.2.0
    gmap: 2017-01-14
    gmes_petap.pl: 4.30
    hisat2: 2.0.5
    hmmscan: HMMER 3.1b2 (February 2015)
    hmmsearch: HMMER 3.1b2 (February 2015)
    kallisto: 0.43.1
    makeblastdb: makeblastdb 2.6.0+
    nucmer: 3.1
    pslCDnaFilter: no way to determine
    rmblastn: rmblastn 2.2.27+
    samtools: samtools 1.5
    tbl2asn: unknown, likely 25.3
    tblastn: tblastn 2.6.0+
        ERROR: fasta not installed
    Checking Environmental Variables...
    $FUNANNOTATE_DB=/usr/local/share/funannotate
    $EVM_HOME=/usr/local/opt/evidencemodeler
    $AUGUSTUS_CONFIG_PATH=/opt/augustus-3.2.1/config
    $GENEMARK_PATH=/usr/local/gmes_petap
    $BAMTOOLS_PATH=/usr/local/opt/bamtools/bin
        ERROR: PASAHOME not set. export PASAHOME=/path/to/dir
        ERROR: TRINITYHOME not set. export TRINITYHOME=/path/to/dir
    -------------------------------------------------------



