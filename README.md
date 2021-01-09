[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/funannotate.svg)](https://github.com/nextgenusfs/funannotate/releases/latest)
[![DOI](https://zenodo.org/badge/48254740.svg)](https://zenodo.org/badge/latestdoi/48254740)
![Conda](https://img.shields.io/conda/dn/bioconda/funannotate)
![Docker Image Size (tag)](https://img.shields.io/docker/image-size/nextgenusfs/funannotate/latest)
![Docker Pulls](https://img.shields.io/docker/pulls/nextgenusfs/funannotate)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/5068)

![Alt text](funannotate-logo.png?raw=true "Funannotate")

funannotate is a pipeline for genome annotation (built specifically for fungi, but will also work with higher eukaryotes). Installation, usage, and more information can be found at [http://funannotate.readthedocs.io](http://funannotate.readthedocs.io)

#### Quickest start Docker:

You can use docker to run `funannotate`. Caveats are that GeneMark is not included in the docker image (see licensing below and you can complain to the developers for making it difficult to distribute/use). I've also written a bash script that can run the docker image and auto-detect/include the proper user/volume bindings.  This docker image is built off of the latest code in master, so it will be ahead of the tagged releases. The image includes the required databases as well, if you want just funannotate without the databases then that is located on docker hub as well `nextgenusfs/funannotate-slim`. So this route can be achieved with:

```
# download/pull the image from docker hub
$ docker pull nextgenusfs/funannotate

# download bash wrapper script (optional)
$ wget -O funannotate-docker https://raw.githubusercontent.com/nextgenusfs/funannotate/master/funannotate-docker

# might need to make this executable on your system
$ chmod +x /path/to/funannotate-docker

# assuming it is in your PATH, now you can run this script as if it were the funannotate executable script
$ funannotate-docker test -t predict --cpus 12
```

#### Quickstart Bioconda install:

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
python -m pip install funannotate
```

To install the most updated code in master you can run:
```
python -m pip install git+https://github.com/nextgenusfs/funannotate.git
```
