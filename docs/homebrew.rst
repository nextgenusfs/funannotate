
.. _homebrew:

HomeBrew mediated Installation
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
While homebrew-science has been deprecated, the replacement is at brewsci, thus the homebrew (Mac) or LinuxBrew (Linux) should still be functional.  You can install most tools with Homebrew/Linuxbrew, however depending on your OS there are a few tools that need to be manually installed. Homebrew https://brew.sh and Linuxbrew http://linuxbrew.sh.

.. code-block:: none
    
    #setup homebrew taps
    brew tap brewsci/bio && brew tap brewsci/science && brew tap nextgenusfs/tap && brew update
    
    #install cpanminus (optional)
    brew install cpanminus
    
    #make sure Perl modules installed
    cpanm Getopt::Long Pod::Usage File::Basename threads threads::shared \
        Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager \
        DBI Text::Soundex Scalar::Util::Numeric Clone JSON LWP::UserAgent DBD::mysql URI::Escape
    
    #install funannotate   
    brew install funannotate
    
**NOTE:** if you are on Mac, Augustus in Homebrew does not compile correctly, you will need to install the version that is located https://github.com/nextgenusfs/augustus


This will automatically install most of the dependencies as well as the most current release of funannotate. Follow the instructions from homebrew, which are:

    1.  Download/install GeneMark-ES/ET: (gmes_petap.pl must be in PATH)
        http://exon.gatech.edu/GeneMark/license_download.cgi

    2.  Install python modules via PIP:

    .. code-block:: none

        pip install -U biopython natsort psutil goatools fisher \
            numpy pandas matplotlib seaborn scikit-learn ete3


    3.  Install RepeatMasker Libraries if you have not done so already.

    .. code-block:: none 
      
        wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz
        tar zxvf RepBaseRepeatMaskerEdition-20170127.tar.gz -C #{HOMEBREW_PREFIX}/opt/repeatmasker/libexec

        cd #{HOMEBREW_PREFIX}/opt/repeatmasker/libexec
        ./configure <config.txt

        #Soft-link a repeatmasker utility script into the PATH:
        ln -s #{HOMEBREW_PREFIX}/opt/repeatmasker/libexec/util/rmOutToGFF3.pl #{HOMEBREW_PREFIX}/bin/rmOutToGFF3.pl
        
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