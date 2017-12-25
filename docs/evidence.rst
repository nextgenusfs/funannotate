
.. _evidence:

Providing evidence to funannotate
==================================

Funannotate uses Evidence Modeler to combine *ab initio* gene model predictions with evidence (transcripts or proteins) aligned to the genome. Therefore, the evidence that you supply at runtime for :code:`--transcript_evidence` and :code:`--protein_evidence` are important. By default, funannotate will use the UniProtKb/SwissProt curated protein database for protein evidence.  However, you can specify other forms of protein evidence, perhaps from a well-annotated closely related species, using the :code:`--protein_evidence` option.  Multiple files can be passed to both :code:`--transcript_evidence` or :code:`--protein_evidence` by separating the files by spaces, for example:

.. code-block:: none

    funannotate predict -i genome.fa -s "Awesome species" --transcript_evidence trinity.fasta myESTs.fa \
        -o output --protein_evidence closely_related.fasta $FUNANNOTATE_DB/uniprot_sprot.fasta
        
You'll notice in this example, I also added the UniProt/SwissProt protein models located in the funannotate database. I should also note that adding protein evidence from ab initio predictors of closely related species should be avoided, this is because those models have not been validated.  What you are trying to do here is to provide the software with high-quality protein models so that information can be used to direct the *ab initio* gene prediction algorithms, so providing them with incorrect/truncated proteins isn't going to help your accuracy and in many cases it may hurt.  It is often okay to just stick with the default UniProtKb/SwissProt protein evidence.

**Sources of Evidence that work well:**

1. De-novo RNA-seq assemblies (i.e. output of Trinity)
2. ESTs (for fungal genomes ESTs from related species can be downloaded from JGI Mycocosm)
3. Curated Protein models from closely related species
