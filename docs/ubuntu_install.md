###Installation on Linux (Ubuntu)

You will need working knowledge of terminal and command line utilities in order to install/run funannotate.  By far the most challenging aspect is installing all of the dependencies correctly. You will also need to have `sudo` privileges to get all of these tools installed.

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
sudo apt-get update
sudo apt-get install -y build-essential
sudo apt-get upgrade -y
sudo apt-get dist-upgrade -y
sudo apt-get install -y git cmake
sudo apt-get install python-dev python-setuptools python-pip
sudo apt-get install libatlas-base-dev libfreetype6-dev libz-dev
sudo apt-get install python-numpy python-scipy python-pandas python-matplotlib python-biopython python-psutil python-sklearn
sudo apt-get install bioperl ncbi-blast+ hmmer gmap bedtools exonerate mummer cpanminus trimal mafft raxml ncbi-tools bamtools
sudo pip install seaborn natsort goatools fisher
```

2) Install perl modules via cpanm
```
sudo cpanm Getopt::Long Pod::Usage File::Basename threads threads::shared \
        Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager
```

3) Install Blat and required kent-tools
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
sudo chmod +x blat
sudo mv blat /usr/local/bin/blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter
sudo chmod +x pslCDnaFilter
sudo mv pslCDnaFilter /usr/local/bin/pslCDnaFilter
```

4) Download and install EVidence Modeler
```
sudo git clone https://github.com/EVidenceModeler/EVidenceModeler.git /usr/local/EVidenceModeler
```

5) Download and install Genome Annotation Generator
```
sudo git clone https://github.com/genomeannotation/GAG.git /usr/local/GAG
```

6) Install Bamtools
```
git clone git://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir build
cd build
cmake ..
make
cd ..
sudo mv bin/bamtools-2.4.0 /usr/local/bin/bamtools
cd ..
rm -r bamtools
```

7) Download and install GeneMark-ES/ET [here](http://exon.gatech.edu/GeneMark/license_download.cgi)
```
#uncompress and then move gmes_petap subdirectory to /usr/local
tar -xzvf $HOME/Downloads/gm_et_linux_64.tar.gz
sudo mv gm_et_linux_64/gmes_petap/ /usr/local

#download key, then move to home directory
gunzip gm_key_64.gz
cp $HOME/Downloads/gm_key ~/.gm_key
```
8) Download and install AUGUSTUS
```
wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz
sudo tar -xvzf augustus.current.tar.gz -C /usr/local
cd /usr/local/augustus-3.2.1
sudo make
```

9) Download and install BRAKER1
```
wget http://exon.gatech.edu/GeneMark/Braker/BRAKER1.tar.gz
sudo tar -xvzf BRAKER1.tar.gz -C /usr/local
```

10) Download and install funannotate
```
sudo git clone https://github.com/nextgenusfs/funannotate.git /usr/local/funannotate
```

11) Add several components and ENV variables to `~/.bash_alias` which will get sourced by your BASHRC - you can use any text editor for this
```
#example using gedit
sudo gedit ~/.bash_alias

#add folders to PATH
export PATH="/usr/local/funannotate:/usr/local/GAG:/usr/local/gmes_petap:/usr/local/BRAKER1:$PATH"

#add environmental variables
export AUGUSTUS_CONFIG_PATH=/usr/local/opt/augustus/libexec/config
export EVM_HOME=/usr/local/EVidenceModeler
export GENEMARK_PATH=/usr/local/gmes_petap
export BAMTOOLS_PATH=/usr/local/Cellar/bamtools/2.4.0/bin
```

12) Re-launch a terminal window (or type `source ~/.bash_profile`). Finally run funannotate setup script to download databases and identify any problems.
```
#navigate into funannotate install directory
cd /usr/local/funannotate

#run setup script
./setup.sh
```
The script will download and format necessary databases and then check all of the dependencies of funannotate - any tool not properly installed will be flagged by the script.

