
.. _install:

Installation
================================

.. toctree::
   :hidden:
  
   dependencies

Funannotate has a lot of dependencies and therefore installation is the most difficult part
of executing the pipeline. The  funannotate pipeline is written in python and can be installed
with pip, i.e. `pip install funannotate`.  You can see a list of :ref:`dependencies`,

To provide an easier installation, funannotate can also be installed with conda.  The recommended
way of doing this is to create a conda environment, for example:


Download/install miniconda and configure the proper channels:

.. code-block:: none
    
    #If you do not have conda, install: miniconda3
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    /bin/bash ~/miniconda.sh -b -p /conda/installation/path
    
    #setup bioconda repository
    conda config --add channels defaults
    conda config --add channels bioconda
    conda config --add channels conda-forge
    
    #then create environment
    conda create -n funannotate python=2.7 funannotate
    
    
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

