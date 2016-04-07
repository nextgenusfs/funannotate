###Installation on Mac OSX

You will need working knowledge of terminal and command line utilities in order to install/run funannotate.  By far the most challenging aspect is installing all of the dependencies correctly.  On Mac OSX, to install these tools we will rely extensively on HomeBrew.  Note HomeBrew is not actually a dependency of funannotate, but for this installation guide it is required.  You will also need to have `sudo` privileges to get all of these tools installed.

####Step-by-step instructions:

1) Install HomeBrew
```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)”
```

2) Setup HomeBrew and Tap Repositories
```
#initial setup
brew doctor

#tap repositories
brew tap homebrew/science
brew tap nextgenusfs/tap
```

3) Install Python modules via pip
```
sudo pip install -U biopython natsort psutil goatools numpy pandas matplotlib seaborn scikit-learn
```

4) Install Perl modules via cpanm (brew install cpanm)
```
sudo cpanm BioPerl Getopt::Long Pod::Usage File::Basename threads threads::shared \
           Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager
```
5) Install funannotation via homebrew
```
#update brew to make sure everything is current
brew update

#now install funannotate and dependencies
brew install funannotate
```

6) Get RepBase data and reconfigure RepeatMasker/RepeatModeler. Register for [RepBase](http://www.girinst.org/repbase/)
```
#download RepeatMasker libraries and install
wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/repeatmaskerlibraries-20150807.tar.gz
tar zxvf repeatmaskerlibraries-20150807.tar.gz -C /usr/local/Cellar/repeatmasker/4.0.5/libexec

#now setup RepeatMasker
cd /usr/local/Cellar/repeatmasker/4.0.5/libexec
./configure <config.txt

#softlink GFF script to bin in path
ln /usr/local/Cellar/repeatmasker/4.0.5/libexec/util/rmOutToGFF3.pl /usr/local/bin
```

7) Download and install GeneMark-ES/ET [here](http://exon.gatech.edu/GeneMark/license_download.cgi)
```
#uncompress and then move gmes_petap subdirectory to /usr/local
tar -xvf $HOME/Downloads/gm_et_macosx.tar
mv gm_et_macosx/gmes_petap/ /usr/local

#download key, then move to home directory
cp $HOME/Downloads/gm_key ~/.gm_key
```

8) Add several components and ENV variables to `~/.bash_profile` which will get sourced by your BASHRC - you can use any text editor for this
```
#add folders to PATH
export PATH="/usr/local/gmes_petap:$PATH"

#add environmental variables
export AUGUSTUS_CONFIG_PATH=/usr/local/Cellar/opt/augustus/libexec/config
export EVM_HOME=/usr/local/Cellar/evidencemodeler/1.1.2
export GENEMARK_PATH=/usr/local/gmes_petap
export BAMTOOLS_PATH=/usr/local/Cellar/bamtools/2.4.0/bin
```

9) Re-launch a terminal window (or type `source ~/.bash_profile`). Finally run funannotate setup script to download databases and identify any problems.
```
#navigate into funannotate install directory
cd /usr/local/Cellar/funannotate/0.1.3/libexec

#run setup script, might need sudo here
./setup.sh
```
The script will download and format necessary databases and then check all of the dependencies of funannotate - any tool not properly installed will be flagged by the script.

####Python Dependencies:
* Python 2
* Biopython
* psutil
* natsort
* goatools
* numpy
* pandas
* matplotlib
* seaborn
* sklearn library

####Software Dependencies:
* HomeBrew
* Perl
* BioPerl
* Blast+
* Hmmer3
* EVidenceModeler
* RepeatModeler
* RepeatMasker
* GMAP
* Blat - if using PASA results to train Augustus
* pslCDnaFilter (kent tools) - if using PASA results to train Augustus
* BedTools
* Augustus
* GeneMark-ES/ET (gmes_petap.pl)
* BamTools
* Genome Annotation Generator (gag.py)
* tbl2asn
* BRAKER1 (optional if training Augustus with RNA-seq data BAM file)
* Mummer
* RAxML
* mafft
* trimal

####Environmental variables:
EVM_HOME, GENEMARK_PATH, BAMTOOLS_PATH, AUGUSTUS_CONFIG_PATH
