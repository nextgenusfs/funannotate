###Installation on Linux (Ubuntu)

You will need working knowledge of terminal and command line utilities in order to install/run funannotate.  By far the most challenging aspect is installing all of the dependencies correctly. You will also need to have `sudo` privileges to get all of these tools installed.  I'm going to use LinuxBrew to install a few of these packages, note there are other ways (notably sudo-apt), however Linux Brew provides up to date binaries for many science based programs and for a non-computer scientist (me) I find it easier to use.  Note, you may have most if not all of these dependencies installed on a server, thus you can skip ahead to install funannotate if you know what you are doing.  The setup script will help you determine what is not installed.

####Step-by-step instructions:


1) install dependencies (instructions for brand new Ubuntu box, i.e. Virutal box)
```
#now install libraries and necessary files
sudo apt-get update
sudo apt-get install -y build-essential
sudo apt-get upgrade -y
sudo apt-get dist-upgrade -y
sudo apt-get install -y git cmake
sudo apt-get install python-dev python-setuptools python-pip
sudo apt-get install libatlas-base-dev libfreetype6-dev libz-dev libboost-iostreams-dev libpng-dev
sudo apt-get install python-numpy python-scipy python-pandas python-matplotlib python-biopython python-psutil python-sklearn
sudo apt-get install bioperl cpanminus
```

2) Install python modules via PIP - some of the packages in apt-get are too old, so upgrade with `-U`
```
sudo pip install -U biopython matplotlib pandas numpy seaborn natsort goatools fisher scikit-learn
```

3) Install perl modules via cpanm or cpan or manually
```
sudo cpanm Getopt::Long Pod::Usage File::Basename threads threads::shared \
        Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager
```

4) Install LinuxBrew
```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Linuxbrew/linuxbrew/go/install)"

#add to ~/.bash_aliases
sudo gedit ~/.bash_aliases

#add the following and save, reload terminal
export PATH="$HOME/.linuxbrew/bin:$PATH"
export MANPATH="$HOME/.linuxbrew/share/man:$MANPATH"
export INFOPATH="$HOME/.linuxbrew/share/info:$INFOPATH"

#setup linuxbrew
brew doctor
brew tap homebrew/dupes
brew tap homebrew/science
brew tap nextgenusfs/tap
```

5) Install funannotate and dependencies using LinuxBrew
```
brew install funannotate
```

6) Download RepeatMasker libraries from [RepBase](http://www.girinst.org/repbase/) you will need to register
```
wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/repeatmaskerlibraries-20150807.tar.gz
tar zxvf repeatmaskerlibraries-20150807.tar.gz -C $HOME/.linuxbrew/Cellar/repeatmasker/4.0.5/libexec

#now setup RepeatMasker
cd $HOME/.linuxbrew/repeatmasker/4.0.5/libexec
./configure <config.txt

#softlink GFF script to bin
ln -s $HOME/.linuxbrew/Cellar/repeatmasker/4.0.5/libexec/util/rmOutToGFF3.pl $HOME/.linuxbrew/bin
```

7) Download and install GeneMark-ES/ET [here](http://exon.gatech.edu/GeneMark/license_download.cgi)
```
#uncompress and then move gmes_petap subdirectory to /usr/local
tar zxvf $HOME/Downloads/gm_et_linux_64.tar.gz
sudo mv gm_et_linux_64/gmes_petap/ /usr/local

#download key, then move to home directory
gunzip gm_key_64.gz
cp $HOME/Downloads/gm_key_64 ~/.gm_key
```

8) Add several components and ENV variables to `~/.bash_aliases` which will get sourced by bashrc
```
#example using gedit
sudo gedit ~/.bash_aliases

#add folder for GeneMark to PATH, HomeBrew will take care of other tools
export PATH="/usr/local/gmes_petap:$PATH"

#add environmental variables
export AUGUSTUS_CONFIG_PATH=$HOME/.linuxbrew/opt/augustus/libexec/config
export EVM_HOME=$HOME/.linuxbrew/Cellar/evidencemodeler/1.1.2
export GENEMARK_PATH=/usr/local/gmes_petap
export BAMTOOLS_PATH=$HOME/.linuxbrew/Cellar/bamtools/2.4.0/bin
```

9) Re-launch a terminal window (or type `source ~/.bash_aliases`). Finally run funannotate setup script to download databases and identify any problems.
```
#navigate into funannotate install directory
cd $HOME/.linuxbrew/Cellar/funannotate/version#/libexec

#run setup script, note you need sudo here to copy over the proper ProteinOrtho version
sudo ./setup.sh
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

