[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/funannotate.svg)](https://github.com/nextgenusfs/funannotate/releases/latest)
[![DOI](https://zenodo.org/badge/48254740.svg)](https://zenodo.org/badge/latestdoi/48254740)
![Conda](https://img.shields.io/conda/dn/bioconda/funannotate)

funannotate is a pipeline for genome annotation (built specifically for fungi, but will also work with higher eukaryotes). Installation, usage, and more information can be found at [http://funannotate.readthedocs.io](http://funannotate.readthedocs.io)

#### Quickstart:

The pipeline can be installed with conda (via [bioconda](https://bioconda.github.io/)):
```
#add appropriate channels
conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

#then create environment
conda create -n funannotate funannotate
```
If you want to use GeneMark-ES/ET you will need to install that manually following developers instructions:
http://topaz.gatech.edu/GeneMark/license_download.cgi

Note that you will need to change the shebang line for all perl scripts in GeneMark to use `/usr/bin/env perl`. 
You will then also need to add `gmes_petap.pl` to the $PATH or set the environmental variable $GENEMARK_PATH to the gmes_petap directory.

To install just the python funannotate package, you can do this with pip:
```
pip install funannotate
```

To install the most updated code in master you can run:
```
python2 -m pip install git+https://github.com/nextgenusfs/funannotate.git
```
