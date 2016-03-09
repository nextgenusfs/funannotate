###Installation on Linux (Ubuntu)

You will need working knowledge of terminal and command line utilities in order to install/run funannotate.  By far the most challenging aspect is installing all of the dependencies correctly. You will also need to have `sudo` privileges to get all of these tools installed.  I'm going to use LinuxBrew to install several of these packages, note there are other ways (notably sudo-apt), however Linux Brew provides up to date binaries for many science based programs and for a non-computer scientist (me) I find it easier to use.

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

####Step-by-step instructions:


1) install dependencies (instructions for brand new Ubuntu box, i.e. Virutal box)
```
#add repository to apt sources
sudo gedit /etc/apt/sources.list
#add following line, save
deb http://us.archive.ubuntu.com/ubuntu vivid main universe

#now install libraries and necessary files
sudo apt-get update
sudo apt-get install -y build-essential
sudo apt-get upgrade -y
sudo apt-get dist-upgrade -y
sudo apt-get install -y git cmake
sudo apt-get install python-dev python-setuptools python-pip
sudo apt-get install libatlas-base-dev libfreetype6-dev libz-dev
sudo apt-get install libboost-iostreams-dev bamtools libbamtools-dev
sudo apt-get install python-numpy python-scipy python-pandas python-matplotlib python-biopython python-psutil python-sklearn
sudo apt-get install bioperl cpanminus exonerate mummer bedtools ncbi-tools-bin
```

2) Install python modules via PIP
```
sudo pip install seaborn natsort goatools fisher
```

3) Install AUGUSTUS
```
wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus-3.2.1.tar.gz
sudo tar -xvfz augustus-3.2.1.tar.gz -C /usr/local
cd /usr/local/augustus-3.2.1
sudo make
sudo make install

#grant write access to config folder
sudo chown -R $(whoami) config/
```

4) Install Blat/pslCDnaFilter
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
sudo chmod +x blat
sudo mv blat /usr/local/bin/blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter
sudo chmod +x pslCDnaFilter
sudo mv pslCDnaFilter /usr/local/bin/pslCDnaFilter
```

5) Install LinuxBrew
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
```

6) Install dependencies using LinuxBrew, here can get most current version of BLAST+ and RAxML
```
brew install blast --without-check
brew install ncurses
brew install hmmer trimal mafft raxml repeatmasker
```

7) Download RepeatMasker libraries from [RepBase](http://www.girinst.org/repbase/) you will need to register
```
#move into repeatmasker folder
cd $HOME/.linuxbrew/repeatmasker/4.0.5/libexec
wget --user name --password pass http://www.girinst.org/server/RepBase/protected/repeatmaskerlibraries/repeatmaskerlibraries-20150807.tar.gz
tar -xzvf repeatmaskerlibraries-20150807.tar.gz

#now setup RepeatMasker, follow prompts
perl configure
```

8) Download install RepeatModeler


9) Install tRNAscan-SE
```
wget http://lowelab.ucsc.edu/software/tRNAscan-SE.tar.gz
sudo tar -xzvf tRNAscan-SE.tar.gz -C /usr/local
cd /usr/local/tRNAscan-SE-1.3.1
make
```

10) Install perl modules via cpanm
```
sudo cpanm Getopt::Long Pod::Usage File::Basename threads threads::shared \
        Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager
```

11) Download and install EVidence Modeler
```
sudo git clone https://github.com/EVidenceModeler/EVidenceModeler.git /usr/local/EVidenceModeler
```

12) Download and install Genome Annotation Generator
```
sudo git clone https://github.com/genomeannotation/GAG.git /usr/local/GAG
```

13) Download and install GeneMark-ES/ET [here](http://exon.gatech.edu/GeneMark/license_download.cgi)
```
#uncompress and then move gmes_petap subdirectory to /usr/local
tar -xzvf $HOME/Downloads/gm_et_linux_64.tar.gz
sudo mv gm_et_linux_64/gmes_petap/ /usr/local

#download key, then move to home directory
gunzip gm_key_64.gz
cp $HOME/Downloads/gm_key ~/.gm_key
```

14) Download and install BRAKER1
```
wget http://exon.gatech.edu/GeneMark/Braker/BRAKER1.tar.gz
sudo tar -xvzf BRAKER1.tar.gz -C /usr/local
```

15) Download and install funannotate
```
sudo git clone https://github.com/nextgenusfs/funannotate.git /usr/local/funannotate
```

16) Add several components and ENV variables to `~/.bash_aliases` which will get sourced by bashrc
```
#example using gedit
sudo gedit ~/.bash_aliases

#add folders to PATH
export PATH="/usr/local/funannotate:/usr/local/GAG:/usr/local/gmes_petap:/usr/local/BRAKER1:/usr/local/tRNAscan-SE-1.3.1:$PATH"

#add environmental variables
export AUGUSTUS_CONFIG_PATH=/usr/local/augustus-3.2.1/config
export EVM_HOME=/usr/local/EVidenceModeler
export GENEMARK_PATH=/usr/local/gmes_petap
export BAMTOOLS_PATH=/usr/local/Cellar/bamtools/2.4.0/bin
```

17) Re-launch a terminal window (or type `source ~/.bash_aliases`). Finally run funannotate setup script to download databases and identify any problems.
```
#navigate into funannotate install directory
cd /usr/local/funannotate

#run setup script, not you need sudo here to copy over the proper ProteinOrtho version
sudo ./setup.sh
```
The script will download and format necessary databases and then check all of the dependencies of funannotate - any tool not properly installed will be flagged by the script.




```

