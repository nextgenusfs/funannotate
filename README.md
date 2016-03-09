# funannotate

funannotate is a pipeline for genome annotation (built specifically for fungi, but could work with other eukaryotes).  Genome annotation is a complicated process that uses software from many sources, thus the hardest part about using funannotate will be getting all of the dependencies installed.  After that, funannotate requires only a few simple commands to go from genome assembly all the way to a functional annotated genome (InterPro, PFAM, MEROPS, CAZymes, GO ontology, etc) that is ready for submission to NCBI.  Moreover, funannotate incorporates a light-weight comparative genomics package that can get you started looking at differences between fungal genomes.

###Installation

funannotate will likely run on any POSIX system, although it has only been tested on Mac OSX and Ubuntu.

* [Mac OSX install instructions](docs/mac_install.md)
* [Ubuntu install instructions](docs/ubuntu_install.md)

###Setup

See installation instructions, but funannotate comes with a shell script to aid you in installation of the external, python, and perl dependencies.  The shell script will also download and format the databases that are required to run funannotate, the databases will occupy quite a bit of space, currently working (uncompressed) is ~ 24 GB.  Whenever possible, funannotate is configured to check external dependencies at runtime, however, `setup.sh` will help you during the initial setup.

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
             
```

###Using funannotate: a simple walkthrough

move into the `sample_data` directory of funannotate.

```
#for example, funannotate installed in /usr/local/funannotate
$ cd /usr/local/funannotate/sample_data

#run funannotate predict on genome 1
$ funannotate predict -i genome1.fasta -o genome1 -s "Genome one" \
    --isolate fun1 --name GN01_ --augustus_species botrytis_cinerea \
    --protein_evidence proteins.fa --transcript_evidence transcripts.fa --cpus 6

```
This command should complete in ~ 5 minutes, will produce an output folder named `genome1` which contains the results.  To save time, here we are using pre-trained botrytis_cinerea to run AUGUSTUS - normally funannotate will train AUGUSTUS for you depending on which input you give it.  

```
#generate functional annotation for genome 1
$ funannotate annotate -i genome1 -e youremail@domain.edu --cpus 6

```
The second command, will add functional annotation to your 66 protein models.  It should complete in ~ 15 minutes - 30 minutes (depending on how long remote query to InterProScan server takes.  The results are in the `annotate_results` folder and have all the necessary files for NCBI WGS submission (.tbl, .sqn, .contigs.fsa, .agp).  A GBK flatfile is also provided.


You can now run similar commands for genome2.fasta and genome3.fasta
```
#predict genome2
$ funannotate predict -i genome2.fasta -o genome2 -s "Genome two" \
    --isolate fun2 --name GN02_ --augustus_species botrytis_cinerea \
    --protein_evidence proteins.fa --transcript_evidence transcripts.fa --cpus 6

#predict genome3
$ funannotate predict -i genome3.fasta -o genome3 -s "Genome three" \
    --isolate fun3 --name GN03_ --augustus_species botrytis_cinerea \
    --protein_evidence proteins.fa --transcript_evidence transcripts.fa --cpus 6

#annotate genome2
$ funannotate annotate -i genome2 -e youremail@domain.edu --cpus 6

#annotate genome3
$ funannotate annotate -i genome3 -e youremail@domain.edu --cpus 6
```

You can now run some "lightweight" comparative genomics on these funannotated genomes:

```
#funannotate comparative genomics
$ funannotate compare

```