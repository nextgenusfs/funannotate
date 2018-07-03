
.. _dependencies:

Dependencies
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Funannotate has a lot of dependencies.  However, it also comes with a few tools to help you get everything installed.  The first is that of :code:`funannotate check`.  You'll see in the output below that the :code:`fasta` tool is missing, which is Bill Pearsons :code:`fasta36` a dependency of the PASA pipeline.  Also the :code:`$PASAHOME`` and :code:`$TRINITYHOME`` variables are not set, that is because on this particular machine they are not installed, i.e. funannotate will alert you at runtime if it is missing a dependency.

.. code-block:: none
    
	$ funannotate check --show-versions
	-------------------------------------------------------
	Checking dependencies for funannotate v1.4.0
	-------------------------------------------------------
	You are running Python v 2.7.11. Now checking python packages...
	biopython: 1.70
	goatools: 0.7.11
	matplotlib: 2.1.1
	natsort: 5.2.0
	numpy: 1.12.1
	pandas: 0.22.0
	psutil: 5.4.3
	requests: 2.18.4
	scikit-learn: 0.19.0
	scipy: 0.19.1
	seaborn: 0.8.1
	All 11 python packages installed


	You are running Perl v 5.026001. Now checking perl modules...
	Bio::Perl: 1.007002
	Carp: 1.42
	Clone: 0.39
	DBD::SQLite: 1.56
	DBD::mysql: 4.046
	DBI: 1.641
	DB_File: 1.84
	Data::Dumper: 2.167
	File::Basename: 2.85
	File::Which: 1.22
	Getopt::Long: 2.5
	Hash::Merge: 0.300
	JSON: 2.97001
	LWP::UserAgent: 6.33
	Logger::Simple: 2.0
	POSIX: 1.76
	Parallel::ForkManager: 1.19
	Pod::Usage: 1.69
	Scalar::Util::Numeric: 0.40
	Storable: 2.62
	Text::Soundex: 3.05
	Thread::Queue: 3.12
	Tie::File: 1.02
	URI::Escape: 3.31
	YAML: 1.24
	threads: 2.21
	threads::shared: 1.58
	All 27 Perl modules installed


	Checking external dependencies...
	RepeatMasker: RepeatMasker 4.0.7
	RepeatModeler: RepeatModeler 1.0.11
	Trinity: 2.5.1
	augustus: 3.2.1
	bamtools: bamtools 2.4.0
	bedtools: bedtools v2.27.1
	blat: BLAT v35
	diamond: diamond 0.9.19
	emapper.py: emapper-1.0.3
	ete3: 3.1.1
	exonerate: exonerate 2.4.0
	fasta: no way to determine
	gmap: 2017-06-20
	gmes_petap.pl: 4.30
	hisat2: 2.1.0
	hmmscan: HMMER 3.1b2 (February 2015)
	hmmsearch: HMMER 3.1b2 (February 2015)
	java: 1.8.0_92
	kallisto: 0.43.1
	mafft: v7.313 (2017/Nov/15)
	makeblastdb: makeblastdb 2.7.1+
	minimap2: 2.10-r761
	nucmer: 3.1
	pslCDnaFilter: no way to determine
	rmblastn: rmblastn 2.2.27+
	samtools: samtools 1.8
	tRNAscan-SE: 1.23 (April 2002)
	tbl2asn: unknown, likely 25.3
	tblastn: tblastn 2.7.1+
	trimal: trimAl v1.4.rev15 build[2013-12-17]
	All 30 external dependencies are installed

	Checking Environmental Variables...
	$FUNANNOTATE_DB=/usr/local/share/funannotate
	$PASAHOME=/Users/jon/software/PASApipeline
	$TRINITYHOME=/usr/local/opt/trinity
	$EVM_HOME=/Users/jon/software/evidencemodeler
	$AUGUSTUS_CONFIG_PATH=/Users/jon/software/augustus/config
	$GENEMARK_PATH=/Users/jon/software/gmes_petap
	$BAMTOOLS_PATH=/Users/jon/software/bamtools-2.4.0/bin
	All 7 environmental variables are set
	-------------------------------------------------------



