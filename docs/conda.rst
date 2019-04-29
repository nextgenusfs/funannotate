
.. _conda:

Conda mediated Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

I'd really like to build a bioconda installation package, but would need some help.  You can however install quite a few of the dependencies with conda.

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
        bedtools blat hmmer exonerate diamond>=0.9 tbl2asn hisat2 ucsc-pslcdnafilter \
        samtools raxml trimal mafft>=7 iqtree kallisto bowtie2 infernal mummer minimap2 \
        trinity>=2.6.6 evidencemodeler pasa>=2.3 codingquarry stringtie gmap=2017.11.15 \
        ete3 salmon>=0.9 jellyfish>=2.2 htslib trnascan-se repeatmasker repeatmodeler \
        trf  perl-threaded perl-db-file perl-bioperl perl-dbd-mysql perl-dbd-sqlite \
        perl-text-soundex perl-scalar-util-numeric perl-data-dumper perl-dbi perl-clone \
        perl-json perl-logger-simple perl-hash-merge perl-yaml perl-pod-usage perl-getopt-long \
        perl-parallel-forkmanager perl-carp perl-app-cpanminus
    
The above will automatically install most of the dependencies, below there are a few manual steps.
    
    1. Download/install Augustus and bamtools
    
        If you are on linux, you can try to install the usual way with conda [can add to install above]:
        
        .. code-block:: none 
        
            conda install -n funannotate bamtools augustus
        
        However, this doesn't always work as there are several compilation issues with augustus depending on your flavor of linux. If the above results in segmentation faults, you'll have to install manually.  Check function outside of funannotate.
        
        If you are on Mac, install v3.2.1 from here: https://github.com/nextgenusfs/augustus
        
    2.  Download/install GeneMark-ES/ET: (gmes_petap.pl must be in PATH)
        http://exon.gatech.edu/GeneMark/license_download.cgi
        
        * make sure to activate the license and move into proper location. you can test proper installation by running `gmes_petap.pl` in the terminal -- you should see help menu
        
    3.  Install RepeatMasker/RepeatModeler  http://www.repeatmasker.org
    
     
    3b. Download Repbase RepeatMasker Libraries if you have not done so already.

    .. code-block:: none 
      
        wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz
        tar zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C /path/to/repeatmasker/location
        cd /path/to/repeatmasker/location
        ./configure

        #Soft-link a repeatmasker utility script into the PATH (may not need to do this depending on install)
        ln -s /path/to/repeatmasker/location/repeatmasker/util/rmOutToGFF3.pl /usr/local/bin/rmOutToGFF3.pl
        
    4. Install Perl modules that aren't avaialble on conda {need somebody to confirm these are still necessary}
    
     .. code-block:: none
        
        #make sure you are in funannotate environment
        conda activate funannotate
     
        #then using cpanminus to install these deps
        cpanm File::Basename Thread::Queue LWP::UserAgent
        


    5. Setup Eggnog-mapper [this is optional but recommended]
    
     .. code-block:: none
        
        #clone the eggnog mapper repo into a location you have read/write access
        git clone https://github.com/jhcepas/eggnog-mapper.git
        
        #move into folder and setup - this will put into eggnog-mapper/data location
        cd eggnog-mapper
        download_eggnog_data.py
        
        #finally add to your funannotate conda env so it is in path when env is activated
        ln -s /path/to/eggnog-mapper/emapper.py /path/to/conda/envs/funannotate/bin/emapper.py
        
	
	NOTE: MacOSX users -- the diamond version shipped with eggnog-mapper needs to be swapped out as the binary provided is compiled on linux. Run a small test with emapper.py to check functionality `emapper.py -m diamond -i test.fa -o test`
    
   
    6. Clone the funannotate repo and add to PATH
    
     .. code-block:: none
     
        git clone https://github.com/nextgenusfs/funannotate.git
        
        #add to PATH
        ln -s /path/to/funannotate/funannotate /path/to/conda/envs/funannotate/bin/funannotate
        
    7. Run funannotate check --show-versions, fix any issues. You will need to export some ENV variables.
    
    .. code-block:: none

        export EVM_HOME=/path/to/conda/envs/funannotate/opt/evidencemodeler-v1.1.1
        export TRINITYHOME=/path/to/conda/envs/funannotate/opt/trinity-2.6.6
        export PASAHOME=/path/to/conda/envs/funannotate/opt/pasa-2.3.3
        export AUGUSTUS_CONFIG_PATH=/path/to/augustus/config
        export GENEMARK_PATH=/path/to/gmes_petap_dir
        export FUNANNOTATE_DB=/path/to/funannotateDB
        
    8.  Setup funannotate databases, specify any location you have read/write access to to `-d` -- this is $FUNANNOTATE_DB

    .. code-block:: none
        
        funannotate setup -d /path/to/DB
        
    9.  If you want these ENV variables to be activated when you activate the conda environment, you can add them as a shell script to the the activate location of your environment, i.e. `/path/to/conda/envs/funannotate/etc/conda/activate.d/` and then you can put the corresponding `unset` commands in the deactivate directory, i.e. `/path/to/conda/envs/funannotate/etc/conda/deactivate.d/`
