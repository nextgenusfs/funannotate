
.. _tutorials:

Tutorials
================================
Funannotate can accommodate a variety of input data and depending on the data you have available you will use funannotate slightly differently, although the core modules are used in the following order:

clean --> sort --> mask --> train --> predict --> update --> annotate --> compare

The following sections will walk-through usage of funannotate for some common data types.

`Genome and RNA sequencing data <http://funannotate.readthedocs.io/en/latest/tutorials.html#genome-assembly-and-rna-seq>`_.
`Genome only <http://funannotate.readthedocs.io/en/latest/tutorials.html#genome-assembly-only>`_.
`Options for non-fungal genomes <http://funannotate.readthedocs.io/en/latest/tutorials.html#non-fungal-genomes-higher-eukaryotes>`_.


Genome assembly and RNA-seq 
-------------------------------------
This the "gold standard" in a sense and if the RNA-seq data is high quality should lead to a high quality annotation.  Paired-end stranded RNA-seq data will yield the best results, although unstranded and/or single ended reads can also be used.  Funannotate can also handle long-read RNA-seq data from PacBio or nanopore sequencers.

For this walkthrough, lets assume you have stranded RNA-seq data from 3 different time points, you have Nanopore direct mRNA sequences, and you have a genome assembly built from Spades.  You have the following files in your folder:

.. code-block:: none

    Spades.genome.fa
    liquid_R1.fq.gz
    liquid_R2.fq.gz
    solid_R1.fq.gz
    solid_R2.fq.gz
    medium_R1.fq.gz
    medium_R2.fq.gz
    nanopore_mRNA.fq.gz
    

1. Haploid fungal genome? This step is optional. Then run :code:`funannotate clean`. Will also run :code:`funannotate sort` to rename fasta headers.

.. code-block:: none

    funannotate clean -i Spades.genome.fa --minlen 1000 -o Spades.genome.cleaned.fa
    

2. Now sort your scaffolds by length and rename with a simple fasta header to avoid downstream problems.

.. code-block:: none

    funannotate sort -i Spades.genome.cleaned.fa -b scaffold -o Spades.genome.cleaned.sorted.fa
    

3. Now we want to softmask the repetitive elements in the assembly.

.. code-block:: none

    funannotate mask -i Spades.genome.cleaned.fa --cpus 12 -o MyAssembly.fa
    
  
4. Now you have a cleaned up/renamed assembly where repeats have been softmasked, run :code:`funannotate train` to align RNA-seq data, run Trinity, and then run PASA.

.. code-block:: none
    
    funannotate train -i MyAssembly.fa -o fun \
        --left liquid_R1.fq.gz solid_R1.fq.gz medium_R1.fq.gz \
        --right liquid_R2.fq.gz solid_R2.fq.gz medium_R2.fq.gz \
        --nanopore_mrna nanopore_mRNA.fq.gz \
        --stranded RF --jaccard_clip --species "Pseudogenus specicus" \
        --strain JMP12345 --cpus 12

You'll notice that I flipped on the :code:`--jaccard_clip` option, since we have a fungal genome we are expected high gene density. This script will run and produce an output directory called :code:`fun` and sub-directory called :code:`training` where it will house the intermediate files. 

5. After training is run, the script will tell you what command to run next, in this case it is :code:`funannotate predict` with the following options:

.. code-block:: none   

    funannotate predict -i MyAssembly.fa -o fun \
        --species "Pseudogenus specicus" --strain JMP12345 \
        --cpus 12

The script will run through the gene prediction pipeline. Note that the scripts will automatically identify and reuse data from :code:`funannotate train`, including using the PASA gene models to train Augustus. If some gene models are unable to be fixed automatically, it will warn you at the end of the script which gene models need to be manually fixed (there might be some errors in tbl2asn I've not seen yet or cannot be fixed without manual intervention).

6. Since we have RNA-seq data, we will use the :code:`funannotate update` command to add UTR data to the predictions and fix gene models that are in disagreement with the RNA-seq data. 

.. code-block:: none  

    funannotate update -i fun --cpus 12
    
Since we ran :code:`funannotate train` those data will be automatically parsed and used to update the UTR data using PASA comparison method. The script will then choose the best gene model at each locus using the RNA-seq data and pseudoalignment with Kallisto. The outputs from this script are located in the :code:`fun/update_results` folder. User will be alerted to any gene models that need to be fixed before moving onto functional annotation.

7. Now we have NCBI compatible gene models, we can now add functional annotation to the protein coding gene models. This is done with the :code:`funannotate annotate` command. But first we want to run InterProScan, Eggnog-mapper, and antiSMASH.

    1. Running InterProScan5.  You could install this locally and run with protein sequences. Otherwise I've built two other options, run from docker or run remotely using EBI servers.

    .. code-block:: none
    
        #run using docker
        funannotate iprscan -i fun -m docker --cpus 12
        
        #run locally (Linux only)
        funannotate iprscan -i fun -m local --iprscan_path /my/path/to/interproscan.sh
        

    2. Now we want to run Eggnog-mapper. You can run this on their webserver http://eggnogdb.embl.de/#/app/emapper or if you have it installed locally then :code:`funannotate annotate` will run it for you.
    
    3. If annotating a fungal genome and you are interested in secondary metabolism gene clusters you can run antiSMASH
    
    .. code-block:: none
    
        funannotate remote -i fun -m antismash -e your-email@domain.edu
    
    4. If you are on a Mac or you don't have phobius installed, you can also run this as a remote search

    .. code-block:: none
    
        funannotate remote -i fun -m phobius -e your-email@domain.edu
        
        #note you could run multiple searches at once
        funannotate remote -i fun -m phobius antismash -e your-email@domain.edu

8. Finally you can run the :code:`funannotate annotate` script incorporating the data you generated.  Passing the funannotate folder will automatically incorporate the interproscan, antismash, phobius results. 

.. code-block:: none    

    funannotate annotate -i fun --cpus 12
    
Your results will be in the :code:`fun/annotate_results` folder.

    
Genome assembly only
------------------------------------- 
If you don't have any RNA-seq data that is okay as you can still generate a high quality annotation using funannotate.  If you are able to get some transcript evidence from closely related species this can also be helpful, if not, funannotate is flexible and can still generate annotation.


1. First we want to softmask the repetitive elements in the assembly.

.. code-block:: none

    funannotate mask -i Spades.assembly.fa --cpus 12 -o MyAssembly.fa
    
  
2. Now you have an assembly where repeats have been softmasked, run :code:`funannotate predict` to find genes.

.. code-block:: none
    
    funannotate predict -i MyAssembly.fa -o fun \
        --species "Pseudogenus specicus" --strain JMP12345 \
        --busco_seed_species botrytis_cinerea --cpus 12

The script will run through the gene prediction pipeline. It will use BUSCO2 to train Augustus and use self-training GeneMark-ES, note the :code:`--busco_seed_species` option which corresponds to a pre-trained parameters for Augustus (:code:`funannotate species` will display the local pre-trained options) - you want to pick a species that is close to the one you are annotating. If some gene models are unable to be fixed automatically, it will warn you at the end of the script which gene models need to be manually fixed (there might be some errors in tbl2asn I've not seen yet or cannot be fixed without manual intervention).

3. Now we have NCBI compatible gene models, we can now add functional annotation to the protein coding gene models. This is done with the :code:`funannotate annotate` command. But first we want to run InterProScan, Eggnog-mapper, and antiSMASH.

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

4. Finally you can run the :code:`funannotate annotate` script incorporating the data you generated.  Passing the funannotate folder will automatically incorporate the interproscan, antismash, phobius results. 

.. code-block:: none    

    funannotate annotate -i fun --cpus 12

Non-fungal genomes (higher eukaryotes)
------------------------------------- 
Since funannotate was originally written for fungal genomes, there are a few default values that you will want to pay attention to if you are not annotating a fungal genome.  

1. Maximum intron length, this parameter is set by default to 3000 bp throughout the scripts, to adjust you can use the :code:`--max_intronlen` flag. 

2. In the :code:`funannotate predict` menu there is a parameter for some fungal specific GeneMark options, these can be turned off by passing :code:`--organism other` at runtime. 

3. In larger genomes (i.e. > 100 MB?) you may get better results to pass the :code:`--repeats2evm` option to :code:`funannotate predict`, this will use the repeat GFF3 file in Evidence Modeler and will reduce the number of gene predictions.  Note you could run the pipeline once without this flag to see the results and then run it again adding the option to compare results. If you see a large discrepancy between GeneMark and Augustus predictions, this seems to be associated with repeat regions (where one of the ab initio predictors gets hung up on repeats), then adding the :code:`--repeats2evm` option will be beneficial.

4. Pay attention to the :code:`--busco_db` option in all scripts. The default is set for :code:`--busco_db dikarya` (default is specifically for dikaryotic fungi). Thus for other organisms :code:`--busco_db` needs to be properly set for each script where it is an option. You can see the available busco databases with the following command:

.. code-block:: none

	$ funannotate database --show_buscos
	-----------------------------
	BUSCO DB tree: (# of models)
	-----------------------------
	eukaryota (303)
		metazoa (978)
			nematoda (982)
			arthropoda (1066)
				insecta (1658)
				endopterygota (2442)
				hymenoptera (4415)
				diptera (2799)
			vertebrata (2586)
				actinopterygii (4584)
				tetrapoda (3950)
				aves (4915)
				mammalia (4104)
			euarchontoglires (6192)
				laurasiatheria (6253)
		fungi (290)
			dikarya (1312)
				ascomycota (1315)
					pezizomycotina (3156)
						eurotiomycetes (4046)
						sordariomycetes (3725)
						saccharomycetes (1759)
							saccharomycetales (1711)
				basidiomycota (1335)
			microsporidia (518)
		embryophyta (1440)
		protists (215)
			alveolata_stramenophiles (234)
