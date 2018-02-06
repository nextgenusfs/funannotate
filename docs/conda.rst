
.. _conda:

Conda mediated Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

I'd really like to build a bioconda installation package, but would need some help.  You can however install quite a few of the dependencies with conda.

.. code-block:: none
    
    #download miniconda2
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O ~/miniconda.sh
    
    #install miniconda2
    /bin/bash ~/miniconda.sh -b -p /conda/installation/path
    
    #now update conda and install dependencies
    conda update -y conda
    conda install -y numpy pandas scipy matplotlib \
        seaborn natsort scikit-learn psutil biopython requests
    
    #install ete3 toolkit and dependencies
    conda install -y -c etetoolkit ete3 ete_toolchain
    
    #install external dependencies with bioconda
    conda install -y -c bioconda blast rmblast goatools fisher \
        bedtools blat hmmer exonerate diamond tbl2asn hisat2 ucsc-pslcdnafilter  \
        samtools raxml trimal mafft kallisto bowtie2 infernal mummer
    
    
This will automatically install most of the dependencies. 

    1.  Download/install GeneMark-ES/ET: (gmes_petap.pl must be in PATH)
        http://exon.gatech.edu/GeneMark/license_download.cgi
    
    2. Download/install Bamtools and Augustus
    
        If you are on linux, you can install the usual way:
        
        .. code-block:: none 
        
            #install bamtools
            wget https://github.com/pezmaster31/bamtools/archive/v2.5.0.tar.gz && \
            tar -zxvf v2.5.0.tar.gz && \
            rm v2.5.0.tar.gz && mv bamtools-2.5.0 bamtools && \
            cd bamtools && mkdir build && cd build && \
            cmake .. && make && sudo make install && \
            cd /usr/include && sudo ln -f -s ../local/include/bamtools/ && \
            cd /usr/lib/ &&  sudo ln -f -s /usr/local/lib/bamtools/libbamtools.* .
            
            #install augustus
            wget http://bioinf.uni-greifswald.de/augustus/binaries/old/augustus-3.2.3.tar.gz && \
            tar -zxvf augustus-3.2.3.tar.gz && rm augustus-3.2.3.tar.gz && \
            mv augustus-3.2.3 augustus && cd augustus && make clean && make
        
        If you are on Mac, install using this version: https://github.com/nextgenusfs/augustus
     
    3. Download install EVM/PASA/Trinity

    .. code-block:: none
    
        #EVidence modeler
        wget https://github.com/nextgenusfs/EVidenceModeler/archive/0.1.3.tar.gz && \
        tar -zxvf 0.1.3.tar.gz && rm 0.1.3.tar.gz && \
        mv EVidenceModeler-0.1.3 evidencemodeler
        
        #Trinity
        wget https://github.com/trinityrnaseq/trinityrnaseq/archive/Trinity-v2.5.1.tar.gz && \
        tar -zxvf Trinity-v2.5.1.tar.gz && rm Trinity-v2.5.1.tar.gz && \
        mv trinityrnaseq-Trinity-v2.5.1 Trinity && cd Trinity && make && make plugins
       
        #PASA
        wget https://github.com/PASApipeline/PASApipeline/archive/pasa-v2.2.0.tar.gz && \
        tar -zxvf pasa-v2.2.0.tar.gz && rm pasa-v2.2.0.tar.gz && \
        mv PASApipeline-pasa-v2.2.0 PASApipeline && cd PASApipeline && make clean && make
        

    4.  Install RepeatMasker/RepeatModeler  http://www.repeatmasker.org
    
     
    4b. Download Repbase RepeatMasker Libraries if you have not done so already.

    .. code-block:: none 
      
        wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz
        tar zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C /path/to/repeatmasker/location
        cd /path/to/repeatmasker/location
        ./configure

        #Soft-link a repeatmasker utility script into the PATH:
        ln -s /path/to/repeatmasker/location/repeatmasker/util/rmOutToGFF3.pl /usr/local/bin/rmOutToGFF3.pl
        
    5. Install Perl modules, i.e. with cpanminus
    
     .. code-block:: none
     
         cpanm Getopt::Long Pod::Usage File::Basename threads threads::shared \
            Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager \
            DBI Text::Soundex Scalar::Util::Numeric Clone JSON LWP::UserAgent DBD::mysql URI::Escape
   
    
    6.  Setup funannotate databases:

    .. code-block:: none
        
        funannotate setup -d /path/to/DB

    7.  Export required ENV variables (your paths might differ slightly):
    
    .. code-block:: none

        export EVM_HOME=/path/to/evidencemodeler
        export AUGUSTUS_CONFIG_PATH=/path/to/augustus/config
        export BAMTOOLS_PATH=/path/to/bamtools/bin
        export GENEMARK_PATH=/path/to/gmes_petap.pl
        export FUNANNOTATE_DB=/path/to/DB