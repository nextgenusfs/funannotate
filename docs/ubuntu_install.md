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

####Environmental variables:
EVM_HOME, GENEMARK_PATH, BAMTOOLS_PATH, AUGUSTUS_CONFIG_PATH

1) installation instructions for funannotate on a brand new Ubuntu box, i.e. Virutal box
```
sudo apt-get update
sudo apt-get install -y build-essential
sudo apt-get upgrade -y
sudo apt-get dist-upgrade -y
sudo apt-get install -y git cmake
sudo apt-get install python-dev python-setuptools python-pip
sudo apt-get install libatlas-base-dev libfreetype6-dev libz-dev
sudo apt-get install python-numpy python-scipy python-pandas python-matplotlib python-biopython python-psutil python-sklearn
sudo apt-get install bioperl ncbi-blast+ hmmer gmap bedtools exonerate mummer cpanminus
sudo pip install seaborn natsort goatools fisher
```

2) Install perl modules via cpanm

#perl modules needed for utilities
sudo cpanm Getopt::Long Pod::Usage File::Basename threads threads::shared Thread::Queue Carp Data::Dumper YAML Hash::Merge Logger::Simple Parallel::ForkManager

3) Install Blat and required kent-tools
```
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/blat/blat
sudo chmod +x blat
mv blat /usr/local/bin/blat
wget http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/pslCDnaFilter
sudo chmod +x pslCDnaFilter
mv pslCDnaFilter /usr/local/bin/pslCDnaFilter
```

4) Install BamTools
```
git clone git://github.com/pezmaster31/bamtools.git
cd bamtools
mkdir build
cd build
cmake ..
make
```

#get funannotate




#bamtools


#GeneMark-ET
#Have to do this manually, accept license, etc
Download from: http://exon.gatech.edu/GeneMark/license_download.cgi
tar xzvf gm_et_linux_64.tar.gz
cp gm_key ~/.gm_key


#Augustus
wget http://bioinf.uni-greifswald.de/augustus/binaries/augustus.current.tar.gz




#do this with shell script???
#install EvidenceModeler
wget https://github.com/EVidenceModeler/EVidenceModeler/archive/v1.1.1.tar.gz

#install BRAKER1
wget http://exon.gatech.edu/genemark/Braker/BRAKER1.tar.gz

tbl2asn