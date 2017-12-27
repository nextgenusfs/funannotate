
.. _conda:

Conda mediated Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
I'd really like to build a bioconda installation package, but would need some help.  You can however install nearly all of the dependencies with conda.

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

    2.  Install python modules via PIP or conda:

    .. code-block:: none

        pip install -U biopython natsort psutil goatools fisher \
            numpy pandas matplotlib seaborn scikit-learn ete3
        
        conda install -c biocondabiopython natsort psutil goatools fisher \
            numpy pandas matplotlib seaborn scikit-learn
            
        conda install -c etetoolkit ete3 ete_toolchain


    3.  Install RepeatMasker/RepeatModler and corresponding Libraries if you have not done so already.

    .. code-block:: none 
      
        wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz
        tar zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C #{HOMEBREW_PREFIX}/opt/repeatmasker/libexec

        cd #{HOMEBREW_PREFIX}/opt/repeatmasker/libexec
        ./configure <config.txt

        #Soft-link a repeatmasker utility script into the PATH:
        ln -s #{HOMEBREW_PREFIX}/opt/repeatmasker/util/rmOutToGFF3.pl #{HOMEBREW_PREFIX}/bin/rmOutToGFF3.pl
        
    4.  Setup funannotate databases:

    .. code-block:: none
        
        funannotate setup -d /path/to/DB

    5.  Export required ENV variables (your paths might differ slightly):
    
    .. code-block:: none

        export EVM_HOME=#{HOMEBREW_PREFIX}/opt/evidencemodeler
        export AUGUSTUS_CONFIG_PATH=#{HOMEBREW_PREFIX}/opt/augustus/libexec/config
        export BAMTOOLS_PATH=#{HOMEBREW_PREFIX}/opt/bamtools/bin
        export GENEMARK_PATH=/path/to/gmes_petap.pl
        export FUNANNOTATE_DB=/path/to/DB