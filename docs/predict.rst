
.. _predict:

Gene Prediction
================================
 
Gene prediction in funannotate is dynamic in the sense that it will adjust based on the input parameters passed to the :code:`funannotate predict` script. At the core of the prediction algorithm is Evidence Modeler, which takes several different gene prediction inputs and outputs consensus gene models. The two *ab initio* gene predictors are Augustus and GeneMark-ES/ET. An important component of gene prediction in funannotate is providing "evidence" to the script, you can read more about :ref:`evidence`. To explain how :code:`funannotate predict` works, I will walk-through a few examples and describe step-by-step what is happening.

**1. Genome fasta file, Trinity transcripts, RNAseq BAM file and PASA/transdecoder data.**

.. code-block:: none

    funannotate predict -i genome.fasta --species "Genome awesomenous" --isolate T12345 \
        --transcript_evidence trinity.fasta --rna_bam alignments.bam --pasa_gff pasa.gff3

- In this example, funannotate will run the following steps:
    1. Mask Repeats. Run RepeatModeler/RepeatMasker
    2. Align Transcript Evidence to genome using GMAP and Blat
    3. Align Protein Evidence to genome using Diamond/Exonerate.
    4. Launch BRAKER to train Augustus and GeneMark-ET
    5. Extract high-quality Augustus predictions (HiQ)
    6. Pass all data to Evidence Modeler and run
    7. Filter gene models (length filtering, spanning gaps, and transposable elements)
    8. Predict tRNA genes using tRNAscan-SE
    9. Generate an NCBI annotation table (.tbl format)
    10. Convert to GenBank format using tbl2asn
    11. Parse NCBI error reports and alert user to invalid gene models


**2. Genome fasta file and EST transcripts.**

.. code-block:: none

    funannotate predict -i genome.fasta --species "Genome awesomenous" --isolate T12345 \
        --transcript_evidence my_ests.fa
        
- Funannotate will now run the following steps:
    1. Mask Repeats. Run RepeatModeler/RepeatMasker
    2. Align Transcript Evidence (ESTs) to genome using GMAP and Blat
    3. Align Protein Evidence to genome using Diamond/Exonerate.
    4. Run GeneMark-ES (self-training) on genome fasta file
    5. Run a modified BUSCO2 script to identify conserved orthologs
    6. Combined GeneMark and BUSCO2 results, feed into Evidence Modeler
    7. Double-check EVM BUSCO2 consensus models are accurate --> use to train Augustus
    8. Run Augustus using training set derived from BUSCO2 orthologs
    9. Extract high-quality Augustus predictions (HiQ)
    10. Pass GeneMark, Augustus, HiQ, transcript align, protein align --> to EVM
    11. Filter gene models (length filtering, spanning gaps, and transposable elements)
    12. Predict tRNA genes using tRNAscan-SE
    13. Generate an NCBI annotation table (.tbl format)
    14. Convert to GenBank format using tbl2asn
    15. Parse NCBI error reports and alert user to invalid gene models
    
**3. Genome fasta file and Maker GFF output. Note this option is provided out of convenience for the user, however, it won't provide the best results.**

.. code-block:: none

    funannotate predict -i genome.fasta --species "Genome awesomenous" --isolate T12345 \
        --maker_gff my_previous_maker.gff


- Funannotate will now run the following steps:
    1. Mask Repeats. Run RepeatModeler/RepeatMasker
    2. Parse --pasa_gff and/or --other_gff
    3. Extract gene models from Maker gff
    4. Pass Maker, pasa, other models --> to EVM
    5. Filter gene models (length filtering, spanning gaps, and transposable elements)
    6. Predict tRNA genes using tRNAscan-SE
    7. Generate an NCBI annotation table (.tbl format)
    8. Convert to GenBank format using tbl2asn
    9. Parse NCBI error reports and alert user to invalid gene models


Explanation of inputs and options:
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
**What are the inputs?**

The simplest way to run :code:`funannotate predict` is to provide a genome fasta file, an output folder, and a species name (binomial), i.e. this would look like:

.. code-block:: none

    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans"
    

**Repeat Masking is too slow or I ran it already...**

You can also pass a :code:`--masked_genome` with a :code:`--repeatmasker_gff3` file to bypass repeat masking, i.e. you have run it already outside of funannotate.  Keep in mind that funannotate is expecting soft-masked genome, where repeats are represented as lowercase nucleotides.

The other way to bypass RepeatModeler generation of a repeat library is to specify a RepeatModeler library or a FASTA file of repeats generated elsewhere to :code:`--repeatmodeler_lib` option.  RepeatModeler can be bypassed by passing a RepeatMasker species option via the :code:`--repeatmasker_species`.

Putting these options to work:

.. code-block:: none

    funannotate predict --masked_genome mygenome.softmasked.fa --repeatmasker_gff3 repeats.gff3 \
        -o output_folder -s "Aspergillus nidulans"
        
    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans" \
        --repeatmodeler_lib myrepeats.lib
        
    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans" \
        --repeatmasker_species fungi
        
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

Evidence Modeler builds consensus gene models and in addition to providing EVM with the predictions/evidence it also requires "weights" for each set of evidence. By default the inputs are set to 1 for *ab initio* predictions and transcript/protein alignments. If high quality gene models from PASA are passed :code:`--pasa_gff`, they default to a weight of 10. While if evidence from another GFF file is passed via :code:`--other_gff` those models are set to 1 by default.  You can control the weight of both the PASA evidence as well as the OTHER evidence by using a semicolon in the input, i.e.

.. code-block:: none
    
    funannotate predict -i mygenome.fa -o output_folder -s "Aspergillus nidulans"
        --pasa_gff mypasamodels.gff3:8 --other_gff prediction.gff3:5
        
**Submitting to NCBI, what should I know?**

Funannotate will produce NCBI/GeneBank-submission ready output, however, there are a few things you should do if planning on submitting to NCBI.

    1. Get a locus_tag number for your genome.
        You do this by starting a WGS genome submission and either specifying a locus tag or one will be assigned to you. The default in funannotate is to use FUN_. 
        
    2. Pre-submission inquiry of unannotated genome.
        If you are new to genome assembly/annotation submission, be aware that your assembly will have to undergo some quality checks before being accepted by NCBI. Sometimes this results in you have to update your assembly, i.e. remove contigs, split contigs where you have adapter contamination, etc. If you have already done your annotation and then have to make these changes it can be very difficult. Instead, you can start your WGS submission and request that the GenBank curators do a quality check on your assembly and fix any problems prior to generating annotation with funannotate. 
    
    3. Generated an SBT template file. https://submit.ncbi.nlm.nih.gov/genbank/template/submission/
    
    

