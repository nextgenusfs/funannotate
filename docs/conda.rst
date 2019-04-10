
.. _conda:

Conda mediated Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

I'd really like to build a bioconda installation package, but would need some help.  You can however install quite a few of the dependencies with conda.

.. code-block:: none
    
    #download miniconda2 or miniconda3, miniconda2 shown
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh
    
    #install miniconda
    /bin/bash ~/miniconda.sh -b -p /conda/installation/path
    
    #now update conda and install dependencies
    conda update -y conda
    conda create -y -n funannotate python=2.7 numpy pandas scipy matplotlib \
        seaborn natsort scikit-learn psutil biopython requests
        
    #activate your environment and keep building
    conda activate funannotate
    
    #install ete3 toolkit and dependencies
    conda install -y -c etetoolkit ete3 ete_toolchain
    
    #install external dependencies with bioconda
    conda install -y -c bioconda blast rmblast goatools fisher \
        bedtools blat hmmer exonerate diamond tbl2asn hisat2 ucsc-pslcdnafilter  \
        samtools raxml trimal mafft kallisto bowtie2 infernal mummer minimap2 \
        trinity evidencemodeler pasa codingquarry
    
    
This will automatically install most of the dependencies. 

    1.  Download/install GeneMark-ES/ET: (gmes_petap.pl must be in PATH)
        http://exon.gatech.edu/GeneMark/license_download.cgi
    
    2. Download/install Bamtools and Augustus
    
        If you are on linux, you can install the usual way:
        
        .. code-block:: none 
        
            conda install augustus
        
        If you are on Mac, install v3.2.1 from here: https://github.com/nextgenusfs/augustus

    3.  Install RepeatMasker/RepeatModeler  http://www.repeatmasker.org
    
     
    3b. Download Repbase RepeatMasker Libraries if you have not done so already.

    .. code-block:: none 
      
        wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz
        tar zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C /path/to/repeatmasker/location
        cd /path/to/repeatmasker/location
        ./configure

        #Soft-link a repeatmasker utility script into the PATH:
        ln -s /path/to/repeatmasker/location/repeatmasker/util/rmOutToGFF3.pl /usr/local/bin/rmOutToGFF3.pl
        
    4. Install Perl modules, i.e. with cpanminus -- or could use conda for many of these
    
     .. code-block:: none
     
         cpanm Getopt::Long Pod::Usage File::Basename threads threads::shared \
            Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager \
            DBI Text::Soundex Scalar::Util::Numeric Clone JSON LWP::UserAgent DBD::mysql URI::Escape
   
    5. Clone the funannotate repo and add to PATH
    
     .. code-block:: none
     
        git clone https://github.com/nextgenusfs/funannotate.git
        
        #add to PATH
        ln -s /path/to/funannotate/funannotate /path/to/conda/envs/funannotate/bin/funannotate
        
    6.  Setup funannotate databases, specify any location you have read/write access to to `-d`

    .. code-block:: none
        
        funannotate setup -d /path/to/DB

    7.  Export required ENV variables (your paths might differ slightly):
    
    .. code-block:: none

        export EVM_HOME=/path/to/evidencemodeler
        export AUGUSTUS_CONFIG_PATH=/path/to/augustus/config
        export BAMTOOLS_PATH=/path/to/bamtools/bin
        export GENEMARK_PATH=/path/to/gmes_petap_dir
        export FUNANNOTATE_DB=/path/to/DB
