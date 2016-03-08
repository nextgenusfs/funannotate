# funannotate

###README and funannotate is still under construction....please stay tuned.

funannotate is a pipeline for genome annotation (built specifically for fungi, but could work with other eukaryotes).  Genome annotation is a complicated process that uses software from many sources, thus the hardest part about using funannotate will be getting all of the dependencies installed.  After that, funannotate requires only a few simple commands to go from genome assembly all the way to a functional annotated genome (InterPro, PFAM, MEROPS, CAZymes, GO ontology, etc) that is ready for submission to NCBI.  Moreover, funannotate incorporates a light-weight comparative genomics package that can get you started looking at differences between fungal genomes.

###Installation

funannotate will likely run on any POSIX system, although it has only been tested on Mac OSX and Ubuntu.

* [Mac OSX install instructions](docs/mac_install.md)
* [Ubuntu install instructions](docs/ubuntu_install.md)

###Setup

See installation instructions, but funannotate comes with a shell script to aid you in installation of the system, python, and perl dependencies.  The shell script will also download and format the databases that are required to run funannotate, the databases will occupy quite a bit of space, currently working (uncompressed) is ~ 24 GB.  Whenever possible, funannotate is configured to check external dependencies at runtime, however, `setup.sh` will help you during the initial setup.

To run the setup script, navigate to the funannotate home directory and type:
```
#to display help menu
$ ./setup.sh -h
    To download databases and check dependencies:   ./setup.sh
    To just download databases:  ./setup.sh db
    To just check dependencies:  ./setup.sh dep

#to check dependencies and download databases
$ ./setup.sh
```

###Funannotate help menu

To see the help menu, simply type `funannotate` in the terminal window.  Similarly, e.g `funannotate predict` without any arguments will give you the options available to pass to each script, this is consistent for all of the funanntoate commands.
```
$  funannotate

Usage:       funannotate <command> <arguments>
version:     0.1.0

Description: Funannotate is a genome prediction, annotation, and comparison pipeline.
    
Command:     clean          Find/remove small repetitive contigs
             sort           Sort by size and rename contig headers (recommended)
             species        list pre-trained Augustus species
             
             predict        Run gene prediction pipeline
             annotate       Assign functional annotation to gene predictions
             compare        Compare funannotated genomes
             
             fix            Remove adapter/primer contamination from NCBI error report
             
Written by Jon Palmer (2016) nextgenusfs@gmail.com
```

