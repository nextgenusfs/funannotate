
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

.. code-block:: none

	conda create -n funannotate funannotate
	

A more detailed installation of miniconda and the proper channels:

.. code-block:: none
    
    #If you do not have conda, install: miniconda3
    wget --quiet https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda.sh
    /bin/bash ~/miniconda.sh -b -p /conda/installation/path
    
    #setup bioconda repository
    conda config --add channels defaults
    conda config --add channels etetoolkit
    conda config --add channels bioconda
    conda config --add channels conda-forge
    
    #then create environment
    conda create -n funannotate python=2.7 funannotate

    


