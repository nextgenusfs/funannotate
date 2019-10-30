
.. _prepare:

Preparing your Assembly
--------------------------------
There are a few things that you can do to your multi-FASTA assembly to get it "ready" to be annotated.  These steps include methods for removing small repetitive contigs from an assembly, sorting/renaming contig headers so they do not cause problems during prediction step, and repeatmasking your assembely (required).


Cleaning your Assembly
================================
When working with haploid assemblies, sometimes you want to remove some repetitive contigs that are contained in other scaffolds of the assembly. If the repeats are indeed unique, then we want to keep them in the assembly. Funannotate can help "clean" up repetitive contigs in your assembly.  This is done using a "leave one out" methodology using minimap2 or mummer (nucmer), where the the shortest contigs/scaffolds are aligned to the rest of the assembly to determine if it is repetitive. The script loops through the contigs starting with the shortest and workings its way to the N50 of the assembly, dropping contigs/scaffolds that are greater than the percent coverage of overlap (:code:`--cov`) and the percent identity of overlap (:code:`--pident`). 
 
.. code-block:: none

    $ funannotate clean

	Usage:       funannotate clean <arguments>
	version:     1.7.0

	Description: The script sorts contigs by size, starting with shortest contigs it uses minimap2
				 to find contigs duplicated elsewhere, and then removes duplicated contigs.
	
	Arguments:   
	  -i, --input    Multi-fasta genome file (Required)
	  -o, --out      Cleaned multi-fasta output file (Required)
	  -p, --pident   Percent identity of overlap. Default = 95
	  -c, --cov      Percent coverage of overlap. Default = 95
	  -m, --minlen   Minimum length of contig to keep. Default = 500
	  --exhaustive   Test every contig. Default is to stop at N50 value.


Sorting/Rename FASTA Headers    
================================
NCBI limits the number of characters in a FASTA header for submission to 16 characters and Augustus also has problems with longer contig/scaffold names. You can use this simple script to sort your assembly by length and then rename the FASTA headers.

.. code-block:: none

    $funannotate sort

	Usage:       funannotate sort <arguments>
	version:     1.7.0

	Description: This script sorts the input contigs by size (longest->shortest) and then relabels
				 the contigs with a simple name (e.g. scaffold_1).  Augustus can have problems with
				 some complicated contig names.
	
	Arguments:   
	  -i, --input    Multi-fasta genome file. (Required)
	  -o, --out      Sorted by size and relabeled output file. (Required)
	  -b, --base     Base name to relabel contigs. Default: scaffold
	  --minlen       Shorter contigs are discarded. Default: 0


.. _repeatmasking

RepeatMasking your Assembly
================================
This is an essential step in the annotation process. As of v1.4.0 repeatmasking has been decoupled from :code:`funannotate predict` in order to make it more flexible and accomodate those users that don't have access to the RepBase library (a requirement of RepeatMasker). The :code:`funannotate mask` command default is to run simple masking using tantan.  The script is a wrapper for RepeatModeler and RepeatMasker, however you can use any external program to softmask your assembly.  Softmasking is where repeats are represented by lowercase letters and all non-repetitive regions are uppercase letters. One alternative to RepeatMasker is RED (REpeat Detector) you can find a wrapper for this program `Redmask <https://github.com/nextgenusfs/redmask>`_.

.. code-block:: none
    
    $funannotate mask
    
	Usage:       funannotate mask <arguments>
	version:     1.7.0

	Description: This script is a wrapper for repeat masking. Default is to run very simple
				 repeat masking with tantan. The script can also run RepeatMasker and/or 
				 RepeatModeler. It will generate a softmasked genome. Tantan is probably not
				 sufficient for soft-masking an assembly, but with RepBase no longer being
				 available RepeatMasker/Modeler may not be functional for many users.
	
	Arguments:   
	  -i, --input                    Multi-FASTA genome file. (Required)
	  -o, --out                      Output softmasked FASTA file. (Required)

	Optional:
	  -m, --method                   Method to use. Default: tantan [repeatmasker, repeatmodeler]
	  -s, --repeatmasker_species     Species to use for RepeatMasker
	  -l, --repeatmodeler_lib        Custom repeat database (FASTA format)
	  --cpus                         Number of cpus to use. Default: 2
	  --debug                        Keep intermediate files
