
.. _tutorials:

Tutorials
================================
Funannotate can accommodate a variety of input data and depending on the data you have available you will use funannotate slightly differently, although the core modules are used in the following order:

clean --> sort --> train --> predict --> update --> annotate --> compare

The following sections will walk-through usage of funannotate for some common data types.

Genome and RNA sequencing data [link]
Genome and EST data [link]
Genome and proteomic data [link]
Genome only [link]


Genome assembly and RNA-seq 
-------------------------------------
This the "gold standard" in a sense and if the RNA-seq data is high quality should lead to a high quality annotation.  Paired-end stranded RNA-seq data will yield the best results, although unstranded and/or single ended reads can also be used.

For this walkthrough, lets assume you have stranded RNA-seq data from 3 different time points and you have a genome assembly built from Spades.  You have the following files in your folder:

.. code-block:: none

    Spades.genome.fa
    liquid_R1.fq.gz
    liquid_R2.fq.gz
    solid_R1.fq.gz
    solid_R2.fq.gz
    medium_R1.fq.gz
    medium_R2.fq.gz
    

1. Haploid fungal genome? This step is optional. Then run :code:`funannotate clean`. Will also run :code:`funannotate sort` to rename fasta headers.

.. code-block:: none

    funannotate clean -i Spades.genome.fa --minlen 1000 -o Spades.genome.cleaned.fa
    

2. Now sort your scaffolds by length and rename with a simple fasta header to avoid downstream problems.

.. code-block:: none

    funannotate sort -i Spades.genome.cleaned.fa -b scaffold -o MyAssembly.fa
    
    
3. Now you have a cleaned up/renamed assembly, run :code:`funannotate train` to align RNA-seq data, run Trinity, and then run PASA.

.. code-block:: none
    
    funannotate train -i MyAssembly.fa -o fun \
        --left liquid_R1.fq.gz solid_R1.fq.gz medium_R1.fq.gz \
        --right liquid_R2.fq.gz solid_R2.fq.gz medium_R2.fq.gz \
        --stranded RF --jaccard_clip --species "Pseudogenus specicus" \
        --strain JMP12345 --cpus 12

You'll notice that I flipped on the :code:`--jaccard_clip` option, since we have a fungal genome we are expected high gene density. This script will run and produce an output directory called :code:`fun` and sub-directory called :code:`training` where it will house the intermediate files. 

4. After training is run, the script will tell you what command to run next, in this case it is :code:`funannotate predict` with the following options:

.. code-block:: none   

    funannotate predict -i MyAssembly.fa -o fun \
        --species "Pseudogenus specicus" --strain JMP12345 \
        --transcript_evidence annotate/training/funannotate_train.trinity-GG.fasta \
        --pasa_gff annotate/training/funannotate_train.pasa.gff3 \
        --rna_bam annotate/training/funannotate_train.coordSorted.bam \
        --cpus 12

The script will run through the gene prediction pipeline. If some gene models are unable to be fixed automatically, it will warn you at the end of the script which gene models need to be manually fixed (there might be some errors in tbl2asn I've not seen yet or cannot be fixed without manual intervention).

5. Since we have RNA-seq data, we will use the :code:`funannotate update` command to add UTR data to the predictions and fix gene models that are in disagreement with the RNA-seq data. 

.. code-block:: none  

    funannotate update -i fun --cpus 12
    
Since we ran :code:`funannotate train` those data will be automatically parsed and used to update the UTR data using PASA comparison method. The script will then choose the best gene model at each locus using the RNA-seq data and pseudoalignment with Kallisto. The outputs from this script are located in the :code:`fun/update_results` folder. User will be alerted to any gene models that need to be fixed before moving onto functional annotation.

6. Now we have NCBI compatible gene models, we can now add functional annotation to the protein coding gene models. This is done with the :code:`funannotate annotate` command. But first we want to run InterProScan, Eggnog-mapper, and antiSMASH.

    1. Running InterProScan5.  You could install this locally and run with protein sequences. Otherwise I've built two other options, run from docker or run remotely using EBI servers.

    .. code-block:: none
    
        #run using docker
        funannotate iprscan -i fun -m docker --cpus 12
        
        #run locally (Linux only)
        funannotate iprscan -i fun -m local --iprscan_path /my/path/to/interproscan.sh
        
        #using remote search
        funannotate remote -i fun -m interproscan -e your-email@domain.edu

    2. Now we want to run Eggnog-mapper. You can run this on their webserver http://eggnogdb.embl.de/#/app/emapper or if you have it installed locally then :code:`funannotate annotate` will run it for you.
    
    3. If annotating a fungal genome and you are interested in secondary metabolism gene clusters you can run antiSMASH
    
    .. code-block:: none
    
        funannotate remote -i fun -m antismash -e your-email@domain.edu
    
    4. If you are on a Mac or you don't have phobius installed, you can also run this as a remote search

    .. code-block:: none
    
        funannotate remote -i fun -m phobius -e your-email@domain.edu
        
        #note you could run multiple searches at once
        funannotate remote -i fun -m phobius antismash -e your-email@domain.edu

Finally you can run the :code:`funannotate annotate` script incorporating the data you generated.

.. code-block:: none    

    funannotate annotate -i fun --iprscan Pseudogenus_specius.proteins.fa.xml \
        --cpus 12 --eggnog basename.emapper.annotations
    
Your results will be in the :code:`fun/annotate_results` folder.

    
Genome assembly only
------------------------------------- 
    
    
    