
.. _predict:

Gene Prediction
================================
 
Gene prediction in funannotate is dynamic in the sense that it will adjust based on the input parameters passed to the :code:`funannotate predict` script. At the core of the prediction algorithm is Evidence Modeler, which takes several different gene prediction inputs and outputs consensus gene models. The *ab initio* gene predictors are Augustus, snap, glimmerHMM, CodingQuarry and GeneMark-ES/ET (optional due to licensing). An important component of gene prediction in funannotate is providing "evidence" to the script, you can read more about :ref:`evidence`. To explain how :code:`funannotate predict` works, I will walk-through a few examples and describe step-by-step what is happening.

Note that as of funannotate v1.4.0, repeat masking is decoupled from :code:`funannotate predict`, thus predict is expecting that your genome input (:code:`-i`) is softmasked multi-FASTA file.  RepeatModeler/RepeatMasker mediated masking is now done with the :code:`funannotate mask` command. You can read more about `repeatmasking <prepare.html#repeatmasking-your-assembly>`__

Explanation of steps in examples:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**1. Genome fasta file, Trinity transcripts, RNAseq BAM file and PASA/transdecoder data.**

.. code-block:: none

    funannotate predict -i genome.fasta --species "Genome awesomenous" --isolate T12345 \
        --transcript_evidence trinity.fasta --rna_bam alignments.bam --pasa_gff pasa.gff3

- In this example, funannotate will run the following steps:
    1. Align Transcript Evidence to genome using minimap2
    2. Align Protein Evidence to genome using Diamond/Exonerate.
    3. Parse BAM alignments generating hints file
    4. Parse PASA gene models and use to train/run Augustus, snap, GlimmerHMM
    5. Extract high-quality Augustus predictions (HiQ)
    6. Run Stringtie on BAM alignments, use results to run CodingQuarry
    7. Pass all data to Evidence Modeler and run
    8. Filter gene models (length filtering, spanning gaps, and transposable elements)
    9. Predict tRNA genes using tRNAscan-SE
    10. Generate an NCBI annotation table (.tbl format)
    11. Convert to GenBank format using tbl2asn
    12. Parse NCBI error reports and alert user to invalid gene models


**2. Genome fasta file and EST transcripts.**

.. code-block:: none

    funannotate predict -i genome.fasta --species "Genome awesomenous" --isolate T12345 \
        --transcript_evidence my_ests.fa
        
- Funannotate will now run the following steps:
    1. Align Transcript Evidence (ESTs) to genome using minimap2
    2. Align Protein Evidence to genome using Diamond/Exonerate.
    3. Run GeneMark-ES (self-training) on genome fasta file
    4. Run a modified BUSCO2 script to identify conserved orthologs
    5. Combined GeneMark and BUSCO2 results, feed into Evidence Modeler
    6. Double-check EVM BUSCO2 consensus models are accurate --> use to train Augustus
    7. Run Augustus using training set derived from BUSCO2 orthologs
    8. Extract high-quality Augustus predictions (HiQ)
    9. Use BUSCO2 training set to train/run snap and GlimmerHMM
    10. Pass GeneMark, Augustus, HiQ, transcript align, protein align --> to EVM
    11. Filter gene models (length filtering, spanning gaps, and transposable elements)
    12. Predict tRNA genes using tRNAscan-SE
    13. Generate an NCBI annotation table (.tbl format)
    14. Convert to GenBank format using tbl2asn
    15. Parse NCBI error reports and alert user to invalid gene models
    

**3. Re-using a parameters JSON file containing training data for ab-initio prediction software. As of v1.7.0 `funannotate species` now saves training data for all of the ab-initio predictors, this can be re-used by using the parameters.json file.  Passing the `-p, --parameters` file will over-rule any existing training sets from `funannotate species`, ie:**

.. code-block:: none

    funannotate predict -i genome.fasta --species "Genome awesomenous" --isolate T12345 \
        -p genome_awesomenous.parameters.json
        
- Funannotate will now run the following steps:
    1. Align Protein Evidence to genome using Diamond/Exonerate.
    3. Run GeneMark (if found) using `mod` HMM file found in `--parameters`
    4. Run Augustus using training parameters found in `--parameters`
    5. Run snap using training parameters found in `--parameters`
    6. Run GlimmerHMM using training parameters found in `--parameters`
    7. Run CodingQuarry (if found) using training parameters found in `--parameters`
    8. Extract high-quality Augustus predictions (HiQ)
    9. Pass gene models and protein evidence --> to EVM
    10. Filter gene models (length filtering, spanning gaps, and transposable elements)
    11. Predict tRNA genes using tRNAscan-SE
    12. Generate an NCBI annotation table (.tbl format)
    13. Convert to GenBank format using tbl2asn
    14. Parse NCBI error reports and alert user to invalid gene models      
        
    
**4. Genome fasta file and Maker GFF output. Note this option is provided out of convenience for the user, however, it won't provide the best results.**

.. code-block:: none

    funannotate predict -i genome.fasta --species "Genome awesomenous" --isolate T12345 \
        --maker_gff my_previous_maker.gff


- Funannotate will now run the following steps:
    1. Parse --pasa_gff and/or --other_gff
    2. Extract gene models from Maker gff
    3. Pass Maker, pasa, other models --> to EVM
    4. Filter gene models (length filtering, spanning gaps, and transposable elements)
    5. Predict tRNA genes using tRNAscan-SE
    6. Generate an NCBI annotation table (.tbl format)
    7. Convert to GenBank format using tbl2asn
    8. Parse NCBI error reports and alert user to invalid gene models

How are repeats used/dealt with:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Repetitive regions are parsed from the softmasked genome fasta file -- these data are then turned into a BED file.  The softmasked genomes are then passed to the *ab initio* predictors Augustus and GeneMark which each have their internal ways of working with the data -- which according to the developers is preferential than hard masking the sequences. 

- :code:`--soft_mask` option controls how GeneMark deals with repetitive regions. By default this set to `2000` which means that GeneMark skips prediction on repeat regions shorter than 2 kb. 

- :code:`--repeats2evm` option passes the repeat GFF3 file to Evidence Modeler. This option is by default turned off this can too stringent for many fungal genomes that have high gene density. You might want to turn this option on for larger genomes or those that have a high repeat content.

- :code:`--repeat_filter` is an option that controls how funannotate filters out repetitive gene models. Default is to use both overlap and blast filtering -- overlap filtering uses the repeat BED file and drops gene models that are more than 90% contained within a repeat region while the blast filtering compares the amino acid sequences to a small database of known transposons.


Explanation of inputs and options:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**What are the inputs?**

The simplest way to run :code:`funannotate predict` is to provide a softmasked genome fasta file, an output folder, and a species name (binomial), i.e. this would look like:

.. code-block:: none

    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans"
           
**I already trained Augustus or training set is available.**

In this case you can use the pre-trained parameters directly which will save a lot of time. To use this option you can see which species are pre-trained on your system with the :code:`funannotate species` option.  Then you can specify which species parameters to use with the :code:`--augustus_species` option.

.. code-block:: none
    
    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans"
        --augustus_species anidulans
        
**I already have Augustus and/or GeneMark predictions.**

You can pass these predictions directly to funannotate using the :code:`--augustus_gff` and the :code:`--genemark_gtf` options. Note you need to run Augustus with the :code:`--stopCodonExcludedFromCDS=False` for it to be properly parsed by funannotate.

.. code-block:: none
    
    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans"
        --augustus_gff augustus.gff --genemark_gtf genemark.gtf

**How can I control the weights given to Evidence Modeler?**

Evidence Modeler builds consensus gene models and in addition to providing EVM with the predictions/evidence it also requires "weights" for each set of evidence. By default the inputs are set to 1 for *ab initio* predictions and transcript/protein alignments. If high quality gene models from PASA are passed :code:`--pasa_gff`, they default to a weight of 6. While if evidence from another GFF file is passed via :code:`--other_gff` those models are set to 1 by default.  You can control the weight of both the PASA evidence as well as the OTHER evidence by using a colon in the input. You now also control the weights for the ab-initio tools by utilizing the `-w, --weights` option i.e.

.. code-block:: none
    
    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans"
        --pasa_gff mypasamodels.gff3:8 --other_gff prediction.gff3:5
        
    #multiple GFF files can be passed to --other_gff
    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans"
        --pasa_gff mypasamodels.gff3:8 --other_gff prediction1.gff3:5 prediction2.gff3:1
        
    #controlling the weights directly
    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans"
    	--weights augustus:2 pasa:8 snap:1 
        
      
Submitting to NCBI, what should I know?
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Funannotate will produce NCBI/GeneBank-submission ready output, however, there are a few things you should do if planning on submitting to NCBI.

    1. Get a locus_tag number for your genome.
        You do this by starting a WGS genome submission and either specifying a locus tag or one will be assigned to you. The default in funannotate is to use "FUN". 
        
    2. Pre-submission inquiry of unannotated genome.
        If you are new to genome assembly/annotation submission, be aware that your assembly will have to undergo some quality checks before being accepted by NCBI. Sometimes this results in you have to update your assembly, i.e. remove contigs, split contigs where you have adapter contamination, etc. If you have already done your annotation and then have to make these changes it can be very difficult. Instead, you can start your WGS submission and request that the GenBank curators do a quality check on your assembly and fix any problems prior to generating annotation with funannotate. 
    
    3. Generated an SBT template file. https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
    
Explanation of the outputs:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
The output of :code:`funannotate predict` is written to the output/predict_results folder, which contains:

+---------------------------------+----------------------------------------------+
| **File Name**                   | **Description**                              |
+---------------------------------+----------------------------------------------+
| Basename.gbk                    | Annotated Genome in GenBank Flat File format |
+---------------------------------+----------------------------------------------+
| Basename.tbl                    | NCBI tbl annotation file                     |
+---------------------------------+----------------------------------------------+
| Basename.gff3                   | Genome annotation in GFF3 format             |
+---------------------------------+----------------------------------------------+
| Basename.scaffolds.fa           | Multi-fasta file of scaffolds                |
+---------------------------------+----------------------------------------------+
| Basename.proteins.fa            | Multi-fasta file of protein coding genes     |
+---------------------------------+----------------------------------------------+
| Basename.transcripts.fa         | Multi-fasta file of transcripts (mRNA)       |
+---------------------------------+----------------------------------------------+
| Basename.discrepency.report.txt | tbl2asn summary report of annotated genome   |
+---------------------------------+----------------------------------------------+
| Basename.error.summary.txt      | tbl2asn error summary report                 |
+---------------------------------+----------------------------------------------+
| Basename.validation.txt         | tbl2asn genome validation report             |
+---------------------------------+----------------------------------------------+
| Basename.parameters.json        | ab-initio training parameters                |
+---------------------------------+----------------------------------------------+


