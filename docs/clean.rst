
.. _clean:

Cleaning your Assembly
================================
 When working with haploid assemblies, sometimes you want to remove some repetitive contigs that are contained in other scaffolds of the assembly. If the repeats are indeed unique, then we want to keep them in the assembly. Funannotate can help "clean" up repetitive contigs in your assembly.  This is done using a "leave one out" methodology using Mummer (nucmer), where the the shortest contigs/scaffolds are aligned to the rest of the assembly to determine if it is repetitive. The script loops through the contigs starting with the shortest and workings its way to the N50 of the assembly, dropping contigs/scaffolds that are greater than the percent coverage of overlap (:code:`--cov`) and the percent identity of overlap (:code:`--pident`). 
 
.. code-block:: none

    $ funannotate clean
    Usage:       funannotate clean <arguments>
    version:     1.0.0

    Description: The script sorts contigs by size, starting with shortest contigs it uses Mummer 
                 to find contigs duplicated elsewhere, and then removes duplicated contigs.
    
    Arguments:   -i, --input    Multi-fasta genome file (Required)
                 -o, --out      Cleaned multi-fasta output file (Required)
                 -p, --pident   Percent identity of overlap. Default = 95
                 -c, --cov      Percent coverage of overlap. Default = 95
                 -m, --minlen   Minimum length of contig to keep. Default = 500
                 --exhaustive   Test every contig. Default is to stop at N50 value.
