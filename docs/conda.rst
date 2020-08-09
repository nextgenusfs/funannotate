
.. _conda:

Conda mediated Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

I'd really like to build a bioconda installation package, but would need some help.  You can however install quite a few of the dependencies with conda.


**If you are on LINUX -- start here:**

.. code-block:: none
    
    #If you do not have conda, install: download miniconda2 or miniconda3, miniconda3 shown
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    /bin/bash ~/miniconda.sh -b -p /conda/installation/path
    
    #setup bioconda repository
    conda config --add channels defaults
    conda config --add channels etetoolkit
    conda config --add channels bioconda
    conda config --add channels conda-forge
    
    #now create a conda environment and install dependencies
    conda create -y -n funannotate python=2.7 numpy pandas scipy matplotlib seaborn \
        natsort scikit-learn psutil biopython requests blast rmblast goatools fisher  \
        bamtools augustus bedtools hmmer exonerate diamond>=0.9 tbl2asn ucsc-pslcdnafilter \
        samtools raxml trimal mafft>=7 iqtree kallisto>=0.46.0 bowtie2 infernal mummer minimap2 blat \
        trinity>=2.6.6 evidencemodeler pasa>=2.3 codingquarry stringtie gmap=2017.11.15 snap \
        ete3 salmon>=0.9 jellyfish>=2.2 htslib trnascan-se hisat2 glimmerhmm \
        trf perl-threaded perl-db-file perl-bioperl perl-dbd-mysql perl-dbd-sqlite \
        perl-text-soundex perl-scalar-util-numeric perl-data-dumper perl-dbi perl-clone \
        perl-json perl-logger-simple perl-hash-merge perl-yaml perl-pod-usage perl-getopt-long \
        perl-parallel-forkmanager perl-carp perl-soap-lite perl-class-inspector perl-app-cpanminus
    
    #if you are going to use remote search also need LWP module (not on conda)
    cpanm LWP
    
**If you are on MacOS X -- start here:**

.. code-block:: none
    
    #If you do not have conda, install: download miniconda2 or miniconda3, miniconda3 shown
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    /bin/bash ~/miniconda.sh -b -p /conda/installation/path
    
    #setup bioconda repository
    conda config --add channels defaults
    conda config --add channels etetoolkit
    conda config --add channels bioconda
    conda config --add channels conda-forge
    
    #now create a conda environment and install dependencies
    conda create -y -n funannotate python=2.7 numpy pandas scipy matplotlib seaborn \
        natsort scikit-learn psutil biopython requests blast rmblast goatools fisher \
        bedtools hmmer exonerate diamond>=0.9 tbl2asn ucsc-pslcdnafilter \
        samtools raxml trimal mafft>=7 iqtree kallisto>=0.46.0 bowtie2 infernal mummer \
        evidencemodeler  gmap=2017.11.15 hisat2 blat minimap2 snap glimmerhmm  \
        ete3 salmon>=0.9 jellyfish>=2.2 htslib trnascan-se codingquarry \
        trf perl-threaded perl-db-file perl-bioperl perl-dbd-mysql perl-dbd-sqlite \
        perl-text-soundex perl-scalar-util-numeric perl-data-dumper perl-dbi perl-clone \
        perl-json perl-logger-simple perl-hash-merge perl-yaml perl-pod-usage perl-getopt-long \
        perl-parallel-forkmanager perl-carp perl-soap-lite perl-class-inspector perl-app-cpanminus
    
    #if you are going to use remote search also need LWP module (not on conda)
    cpanm LWP

    
MacOSX: Need to install bamtools/augustus/trinity/pasa manually:

Install bamtools/Augustus from here: https://github.com/nextgenusfs/augustus

Trinity: https://github.com/trinityrnaseq/trinityrnaseq

PASA: https://github.com/PASApipeline/PASApipeline
    
    
**The above will automatically install most of the dependencies, below there are a few manual steps.**
        
    1.  Download/install GeneMark-ES/ET: (gmes_petap.pl must be in PATH)
        http://exon.gatech.edu/GeneMark/license_download.cgi
        
        * make sure to activate the license and move into proper location. you can test proper installation by running `gmes_petap.pl` in the terminal -- you should see help menu. Be careful of the shebang line, default is `/usr/bin/perl` which most likely is not what you want, more appropriate is `/usr/bin/env perl`
        
    2.  Install RepeatMasker/RepeatModeler  http://www.repeatmasker.org
    
     
    2b. Download Repbase RepeatMasker Libraries if you have not done so already.

    .. code-block:: none 
      
        wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz
        tar zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C /path/to/repeatmasker/location
        cd /path/to/repeatmasker/location
        ./configure

        #Soft-link a repeatmasker utility script into the PATH (may not need to do this depending on install)
        ln -s /path/to/repeatmasker/location/repeatmasker/util/rmOutToGFF3.pl /usr/local/bin/rmOutToGFF3.pl


    3. Setup Eggnog-mapper v4.5 or v5.0 [v5.0 is not being parsed properly yet in v1.5.3]
    
     .. code-block:: none
        
        #clone the eggnog mapper repo into a location you have read/write access
        git clone https://github.com/jhcepas/eggnog-mapper.git
        
        #move into folder and setup - this will put into eggnog-mapper/data location
        cd eggnog-mapper
        download_eggnog_data.py
        
        #finally add to your funannotate conda env so it is in path when env is activated
        ln -s /path/to/eggnog-mapper/emapper.py /path/to/conda/envs/funannotate/bin/emapper.py
        
	
	NOTE: MacOSX users -- the diamond version shipped with eggnog-mapper needs to be swapped 
	out as the binary provided is compiled on linux. Run a small test with emapper.py to check 
	functionality `emapper.py -m diamond -i test.fa -o test`
    
   
    4. Clone the funannotate repo and add to PATH
    
     .. code-block:: none
     
        git clone https://github.com/nextgenusfs/funannotate.git
        
        #add to PATH
        ln -s /path/to/funannotate/funannotate /path/to/conda/envs/funannotate/bin/funannotate
        
    5. Run funannotate check --show-versions, fix any issues. You will need to export some ENV variables.
    
    .. code-block:: none

        export EVM_HOME=/path/to/conda/envs/funannotate/opt/evidencemodeler-v1.1.1
        export TRINITYHOME=/path/to/conda/envs/funannotate/opt/trinity-2.6.6
        export PASAHOME=/path/to/conda/envs/funannotate/opt/pasa-2.3.3
        export AUGUSTUS_CONFIG_PATH=/path/to/augustus/config
        export GENEMARK_PATH=/path/to/gmes_petap_dir
        export FUNANNOTATE_DB=/path/to/funannotateDB
        
    6.  Setup funannotate databases, specify any location you have read/write access to to `-d` -- this is $FUNANNOTATE_DB

    .. code-block:: none
        
        funannotate setup -d /path/to/DB
        
    7.  If you want these ENV variables to be activated when you activate the conda environment, you can add them as a shell script to the the activate location of your environment, i.e. `/path/to/conda/envs/funannotate/etc/conda/activate.d/` and then you can put the corresponding `unset` commands in the deactivate directory, i.e. `/path/to/conda/envs/funannotate/etc/conda/deactivate.d/`
