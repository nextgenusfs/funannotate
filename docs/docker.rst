
.. _docker:

Docker Installation
================================
Docker is a solution where most of the dependencies are installed and you can start annotating
right away. Because some software and data require individual licensing, the core components
of funannotate are packaged into a docker container, but you must download a few things and
run a docker build locally to get a working container. Note that Eggnog-mapper is not installed
in Docker container as the databases were too large.

1) Download 64 bit Linux GeneMark-ET/ES Key (gm_key_64.gz) from http://exon.gatech.edu/Genemark/license_download.cgi


2) Download RepeatMasker libraries. Register for username at RepBase http://www.girinst.org/repbase/. You can then download the RepeatMasker Libraries most recent version, alternatively can download from command line like so:

.. code-block:: none

    wget --user name --password pass \
        https://www.girinst.org/server/archive/RepBase23.09/protected/repeatmaskerlibraries/RepBaseRepeatMaskerEdition-20170127.tar.gz
    
3) Get SignalP4.1 for linux 64 from CBS http://www.cbs.dtu.dk/cgi-bin/sw_request?signalp


4) Download Dockerfile:

.. code-block:: none

    wget https://raw.githubusercontent.com/nextgenusfs/funannotate/1.5.1/dockerbuild/Dockerfile

5) You should now have the following files in the same directory:

.. code-block:: none
    
    Dockerfile
    gm_key_64.gz
    RepBaseRepeatMaskerEdition-20170127.tar.gz
    signalp-4.1f.Linux.tar.gz

Now you can Build the docker container, which will setup the remaining tools and then download and format funannotate databases.:

.. code-block:: none

    docker build -t funannotate -f Dockerfile .
    

**Running the Docker container with your data:**

In order to run the docker container, you need to put all the files you will use for input to funannotate into the same folder, you can then launch the Docker container and mount your current folder with the following command:

.. code-block:: none

    #container is deleted after you exit
    docker run -it --rm -v $PWD:/home/linuxbrew/data funannotate
    
    #keep container, i.e. mysql databases generated, however will take a lot of HD space
    docker run -it -v $PWD:/home/linuxbrew/data funannotate

This will bring you to a bash prompt within the docker container where all dependencies are installed, so you can now issue the funannotate commands on your data. 

**Limitations with Docker:**

The funannotate docker image does not contain Eggnog-mapper because the databases sizes are too large (> 20 GB).  Eggnog-mapper is an important component of functional annotation, you can run this on the eggnog-mapper webserver and pass results to funannotate or perhaps set up an additional docker image running the eggnog-mapper software.

**Mac OSX users:**

The default storage-driver on docker for Mac is the overlay2 driver.  This driver seems to be incompatible with running/launching MySQL, thus if you are getting errors running funannotate you will need to change your storage-driver to "aufs".  This can be done in Docker preferences, Daemon tab, Advanced tab, and then change the storage-driver.  **Note this will delete all Docker images/containers on your virtual disk.**

.. code-block:: none

  {
  "storage-driver" : "aufs",
  "debug" : true,
  "experimental" : true
  }
