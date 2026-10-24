
.. _install:

Installation
================================

.. toctree::
   :hidden:

   dependencies

Funannotate has a lot of dependencies and therefore installation is the most difficult part
of executing the pipeline. The  funannotate pipeline is written in python and can be installed
with pip, i.e. `pip install funannotate`.  You can see a list of :ref:`dependencies`,

Quickest start Docker:
----------------------

You can use docker to run `funannotate`. The image is published on the GitHub Container Registry and is rebuilt for each release: `ghcr.io/nextgenusfs/funannotate:v1.9.0` (the tags `1.9.0` and `latest` point to the same release). It contains funannotate, the Rust-optimized EVM, PASA and Trinity, and MariaDB for PASA. It does **not** contain the funannotate databases, the eggNOG database, GeneMark or SignalP (see licensing below). Download the databases once with `funannotate setup` and mount that folder into the container. I've also written a bash script that runs the image and detects the user and volume bindings; it uses the image above and mounts `$FUNANNOTATE_DB` when that variable is set. This route can be achieved with:

.. code-block:: none

    # download/pull the release image from the GitHub Container Registry
    $ docker pull ghcr.io/nextgenusfs/funannotate:v1.9.0

    # Apptainer/Singularity
    $ apptainer pull funannotate-1.9.0.sif docker://ghcr.io/nextgenusfs/funannotate:v1.9.0

    # download bash wrapper script (optional)
    $ wget -O funannotate-docker https://raw.githubusercontent.com/nextgenusfs/funannotate/master/funannotate-docker

    # might need to make this executable on your system
    $ chmod +x /path/to/funannotate-docker

    # assuming it is in your PATH, now you can run this script as if it were the funannotate executable script
    $ export FUNANNOTATE_DB=/path/to/funannotate_db   # made with `funannotate setup`
    $ funannotate-docker test -t predict --cpus 12


Quickstart Bioconda install:
----------------------------

.. note::

    **funannotate 1.9.0:** the bioconda package is built separately and may be an older release than this one.
    For 1.9.0 follow :ref:`conda` (a conda environment for the external tools, ``pip install funannotate``, and
    scripts that build the Rust-optimized EVM, PASA and Trinity), or use the container image (:ref:`containers`).
    The commands below install whatever version bioconda currently has.

The pipeline can be installed with conda (via [bioconda](https://bioconda.github.io/) ):

.. code-block:: none

    #add appropriate channels
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge

    #then create environment
    conda create -n funannotate "python>=3.10,<3.12" funannotate

If `conda` is taking forever to solve the environment, I would recommend giving [mamba](https://github.com/mamba-org/mamba) a try:

.. code-block:: none

    #install mamba into base environment
    conda install -n base mamba

    #then use mamba as drop in replacmeent
    mamba create -n funannotate funannotate


If you want to use GeneMark-ES/ET you will need to install that manually following developers instructions:
http://topaz.gatech.edu/GeneMark/license_download.cgi

Note that you will need to change the shebang line for all perl scripts in GeneMark to use `/usr/bin/env perl`.
You will then also need to add `gmes_petap.pl` to the $PATH or set the environmental variable $GENEMARK_PATH to the gmes_petap directory.

To install just the python funannotate package, you can do this with pip:

.. code-block:: none

    python -m pip install funannotate

To install the most updated code in master you can run:

.. code-block:: none

    python -m pip install git+https://github.com/nextgenusfs/funannotate.git



Please setup database and test your installation locally using the following:

.. code-block:: none

	#start up conda ENV
	conda activate funannotate

	#check that all modules are installed
	funannotate check --show-versions

	#download/setup databases to a writable/readable location
	funannotate setup -d $HOME/funannotate_db

	#set ENV variable for $FUNANNOTATE_DB
	echo "export FUNANNOTATE_DB=$HOME/funannotate_db" > /conda/installation/path/envs/funannotate/etc/conda/activate.d/funannotate.sh
	echo "unset FUNANNOTATE_DB" > /conda/installation/path/envs/funannotate/etc/conda/deactivate.d/funannotate.sh

	#run tests -- requires internet connection to download data
	funannotate test -t all --cpus X

