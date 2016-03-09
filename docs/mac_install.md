###Installation on Mac OSX

You will need working knowledge of terminal and command line utilities in order to install/run funannotate.  By far the most challenging aspect is installing all of the dependencies correctly.  On Mac OSX, to install these tools we will rely extensively on HomeBrew.  Note HomeBrew is not actually a dependency of funannotate, but for this installation guide it is required.  You will also need to have `sudo` privileges to get all of these tools installed.

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

####Environmental variables:
EVM_HOME, GENEMARK_PATH, BAMTOOLS_PATH, AUGUSTUS_CONFIG_PATH

####Step-by-step instructions:

1) Install HomeBrew
```
ruby -e "$(curl -fsSL https://raw.githubusercontent.com/Homebrew/install/master/install)”
```
Then setup homebrew: type `brew doctor`, then type: `brew tap homebrew/science`

2) Install some tools via homebrew
```
#install python via homebrew
brew install blast hmmer bedtools bamtools blat gmap-gsnap repeatmodeler repeatmasker cpanm exonerate kent-tools bamtools mummer tbl2asn trnascan raxml

#use pip to install/update python modules
sudo pip install -U biopython natsort psutil goatools numpy pandas matplotlib seaborn scikit-learn
```
