[![Latest Github release](https://img.shields.io/github/release/nextgenusfs/funannotate.svg)](https://github.com/nextgenusfs/funannotate/releases/latest)
[![DOI](https://zenodo.org/badge/48254740.svg)](https://zenodo.org/badge/latestdoi/48254740)
![Conda](https://img.shields.io/conda/dn/bioconda/funannotate)
![Docker Image Size (tag)](https://img.shields.io/docker/image-size/nextgenusfs/funannotate/latest)
![Docker Pulls](https://img.shields.io/docker/pulls/nextgenusfs/funannotate)
[![https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg](https://www.singularity-hub.org/static/img/hosted-singularity--hub-%23e32929.svg)](https://singularity-hub.org/collections/5068)
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat)](http://bioconda.github.io/recipes/funannotate/README.html)
[![European Galaxy server](https://img.shields.io/badge/usegalaxy-.eu-brightgreen?logo=data:image/png;base64,iVBORw0KGgoAAAANSUhEUgAAABgAAAASCAYAAABB7B6eAAAABGdBTUEAALGPC/xhBQAAACBjSFJNAAB6JgAAgIQAAPoAAACA6AAAdTAAAOpgAAA6mAAAF3CculE8AAAACXBIWXMAAAsTAAALEwEAmpwYAAACC2lUWHRYTUw6Y29tLmFkb2JlLnhtcAAAAAAAPHg6eG1wbWV0YSB4bWxuczp4PSJhZG9iZTpuczptZXRhLyIgeDp4bXB0az0iWE1QIENvcmUgNS40LjAiPgogICA8cmRmOlJERiB4bWxuczpyZGY9Imh0dHA6Ly93d3cudzMub3JnLzE5OTkvMDIvMjItcmRmLXN5bnRheC1ucyMiPgogICAgICA8cmRmOkRlc2NyaXB0aW9uIHJkZjphYm91dD0iIgogICAgICAgICAgICB4bWxuczp0aWZmPSJodHRwOi8vbnMuYWRvYmUuY29tL3RpZmYvMS4wLyI+CiAgICAgICAgIDx0aWZmOlJlc29sdXRpb25Vbml0PjI8L3RpZmY6UmVzb2x1dGlvblVuaXQ+CiAgICAgICAgIDx0aWZmOkNvbXByZXNzaW9uPjE8L3RpZmY6Q29tcHJlc3Npb24+CiAgICAgICAgIDx0aWZmOk9yaWVudGF0aW9uPjE8L3RpZmY6T3JpZW50YXRpb24+CiAgICAgICAgIDx0aWZmOlBob3RvbWV0cmljSW50ZXJwcmV0YXRpb24+MjwvdGlmZjpQaG90b21ldHJpY0ludGVycHJldGF0aW9uPgogICAgICA8L3JkZjpEZXNjcmlwdGlvbj4KICAgPC9yZGY6UkRGPgo8L3g6eG1wbWV0YT4KD0UqkwAAAn9JREFUOBGlVEuLE0EQruqZiftwDz4QYT1IYM8eFkHFw/4HYX+GB3/B4l/YP+CP8OBNTwpCwFMQXAQPKtnsg5nJZpKdni6/6kzHvAYDFtRUT71f3UwAEbkLch9ogQxcBwRKMfAnM1/CBwgrbxkgPAYqlBOy1jfovlaPsEiWPROZmqmZKKzOYCJb/AbdYLso9/9B6GppBRqCrjSYYaquZq20EUKAzVpjo1FzWRDVrNay6C/HDxT92wXrAVCH3ASqq5VqEtv1WZ13Mdwf8LFyyKECNbgHHAObWhScf4Wnj9CbQpPzWYU3UFoX3qkhlG8AY2BTQt5/EA7qaEPQsgGLWied0A8VKrHAsCC1eJ6EFoUd1v6GoPOaRAtDPViUr/wPzkIFV9AaAZGtYB568VyJfijV+ZBzlVZJ3W7XHB2RESGe4opXIGzRTdjcAupOK09RA6kzr1NTrTj7V1ugM4VgPGWEw+e39CxO6JUw5XhhKihmaDacU2GiR0Ohcc4cZ+Kq3AjlEnEeRSazLs6/9b/kh4eTC+hngE3QQD7Yyclxsrf3cpxsPXn+cFdenF9aqlBXMXaDiEyfyfawBz2RqC/O9WF1ysacOpytlUSoqNrtfbS642+4D4CS9V3xb4u8P/ACI4O810efRu6KsC0QnjHJGaq4IOGUjWTo/YDZDB3xSIxcGyNlWcTucb4T3in/3IaueNrZyX0lGOrWndstOr+w21UlVFokILjJLFhPukbVY8OmwNQ3nZgNJNmKDccusSb4UIe+gtkI+9/bSLJDjqn763f5CQ5TLApmICkqwR0QnUPKZFIUnoozWcQuRbC0Km02knj0tPYx63furGs3x/iPnz83zJDVNtdP3QAAAABJRU5ErkJggg==)](https://usegalaxy.eu/root?tool_id=funannotate_annotate)


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
conda create -n funannotate "python>=3.6,<3.9" funannotate
```
If `conda` is taking forever to solve the environment, I would recommend giving [mamba](https://github.com/mamba-org/mamba) a try:
```
#install mamba into base environment
conda install -n base mamba

#then use mamba as drop in replacmeent
mamba create -n funannotate funannotate
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
