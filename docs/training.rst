.. _training:

Training data selection
================================

:code:`funannotate predict` trains the *ab initio* predictors Augustus and snap from a set of gene models. These models come from one of three sources: pre-trained parameters, PASA/TransDecoder models built from RNA-seq by :code:`funannotate train`, or BUSCO gene models. Wrong training models teach the predictors wrong gene structure, and the final gene set then loses genes. This page describes how funannotate chooses the training data, which thresholds it applies, and how to check what happened in a run.

The design and its validation against RefSeq annotations are described in `training_data_selection_methods.md <https://github.com/nextgenusfs/funannotate/blob/master/docs/training_data_selection_methods.md>`__ (in this docs folder).

Decision order
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

**funannotate train**

1. **RNA-seq concordance gate** (:code:`--min_rnaseq_map_rate`, default 10). The first :code:`--rnaseq_gate_reads` (default 200,000) reads are mapped to the genome with minimap2 (:code:`-x splice:sr`). If fewer than 10% map with MAPQ ≥ 1, the reads probably come from another organism (infected host tissue, a mislabeled species, an out-of-date file). train then stops before Trinity and PASA with **exit code 3**, so a workflow can predict without RNA-seq. Report: :code:`logfiles/train_rnaseq_gate.tsv`.
   The same reads also give the **median read identity** to the genome (1 − NM / aligned bases, over the mapped reads) and its 10th percentile. train only measures and records identity; it never stops on it. predict applies it (step 5). A copy of the report is written to :code:`training/funannotate_train.rnaseq_gate.tsv`, so predict finds it when a pipeline keeps only the :code:`training` folder.
2. **PASA options.** :code:`--pasa_unspliced_join_spliced` and :code:`--pasa_one_alignment_per_cdna` are passed to PASA only if the installed PASA supports them (PASApipeline ≥ v2.6.1-rc.2). Otherwise funannotate warns and continues without them.
3. **One PASA model per locus.** Same-strand models whose transcript spans overlap by at least :code:`--pasa_alignment_overlap` (30%) of the shorter model form one locus (transitive clusters). The kept model is ranked by: complete ORF, number of CDS exons, CDS length, transcript abundance (TPM), then gene ID.

**funannotate predict**

4. **Initial training source** for each predictor: pre-trained parameters (:code:`-p`, :code:`--augustus_species`), PASA models (:code:`--pasa_gff`), or BUSCO models.
5. **PASA training-set gate** (:code:`--min_pasa_complete_models`, default 500). predict counts complete-ORF models in the PASA GFF3: ATG start, one stop codon at the end, CDS length divisible by 3. Below the threshold, the predictors train from BUSCO instead, and the PASA models are still used as EVM evidence. The gate is skipped when no predictor trains from PASA, or when Augustus output already exists. Report: :code:`logfiles/predict_training_gate.tsv`.
   The **RNA-seq identity gate** (:code:`--min_rnaseq_identity`, default 0 = off) is checked at the same point. If the median read identity measured by train is below the threshold, the reads probably come from another strain or a related species. They map well, but PASA models built from them carry wrong gene structures. The predictors then train from BUSCO. The RNA-seq alignments, PASA models and transcripts are still used as hints and EVM evidence. The gate is skipped when train did not record identity (train not run in this output folder, an older train, or :code:`--min_rnaseq_map_rate 0`). It is off by default. In 40 RefSeq genomes, read identity did not predict whether PASA or BUSCO training was better (Spearman ρ = 0.18); the number of complete PASA models did (see :ref:`assessment_pasa2.6_fun1.9`). Set a threshold when you know the reads come from a related species.
6. **Training-set selection**, in this order:

   - partial ORFs are removed;
   - models whose introns match RNA-seq/protein hints ("keepers") are preferred when at least 200 exist;
   - multi-exon models are required when at least 200 exist;
   - single-exon models with protein support are admitted up to a cap (step 7);
   - redundant models are removed (DIAMOND, ≥ 80% identity and coverage);
   - one model is kept per overlapping cluster (either strand).

7. **Single-exon training genes** (on by default; :code:`--no_training_single_exon` turns this off). Complete single-exon models whose CDS is at least 80% covered by a same-strand protein2genome alignment are admitted. They are capped at share / (1 − share) × the number of multi-exon training models. The share is the fraction of single-CDS GeneMark-ES models (fallback: near-full-length protein alignments).
8. **Minimum training models** (:code:`--min_training_models`). If too few models remain, training falls back to BUSCO.

Reading the decision log
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Every decision above is written to :code:`logfiles/training_decisions.tsv`. train and predict both append to this file. Each decision is also written to the run log as a line starting with :code:`TRAINING-DECISION`, and predict logs a summary table (:code:`Training data decisions:`) before Augustus training.

Columns: :code:`command`, :code:`stage`, :code:`decision`, :code:`value`, :code:`threshold`, :code:`outcome`, :code:`reason`, :code:`timestamp`.

Example: a genome with too few complete PASA models trains from BUSCO:

.. code-block:: none

    stage                  decision                              value              threshold  outcome
    pasa_gate              complete-ORF models in the PASA GFF3  93                 >=500      BUSCO training
    training_mode_switch   predictors trained from PASA          switched to BUSCO             BUSCO training
    busco_training         BUSCO gene models validated           647                           BUSCO training models
    final_training_source  augustus trains from                  busco                         647 models in busco.final.gff3

Example: a genome that trains from PASA:

.. code-block:: none

    pasa_gate                  complete-ORF models in the PASA GFF3  2846  >=500           PASA training
    select_complete_orf        complete-ORF PASA models              2846  of 3,693 input  847 partial models removed
    select_keeper_filter       filterGeneMark keepers                1408  >=200           keeper filter ON
    select_single_exon_admit   single-exon models admitted           38    cap 289         38 of 38 supported admitted
    select_final               PASA training models                  1160                  written to final_training_models.gff3
    final_training_source      augustus trains from                  pasa                  1,160 models in final_training_models.gff3

Stage names: :code:`rnaseq_gate`, :code:`rnaseq_identity`, :code:`pasa_options` (train); :code:`training_mode_initial`, :code:`rnaseq_identity_gate`, :code:`pasa_gate`, :code:`pasa_gate_override`, :code:`training_mode_switch`, :code:`busco_training`, :code:`single_exon_share`, :code:`select_complete_orf`, :code:`select_keeper_filter`, :code:`select_multi_cds`, :code:`select_single_exon_support`, :code:`select_single_exon_admit`, :code:`select_redundancy`, :code:`select_overlap`, :code:`select_final`, :code:`min_training_models`, :code:`final_training_source` (predict).

Options
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. code-block:: none

    funannotate train
      --min_rnaseq_map_rate FLOAT         RNA-seq concordance gate, % mapped reads (default 10; 0 disables)
      --rnaseq_gate_reads INT             reads sampled for the gate (default 200000)
      --pasa_unspliced_join_spliced       PASA opt-in (needs PASApipeline >= v2.6.1-rc.2)
      --pasa_one_alignment_per_cdna       PASA opt-in (needs PASApipeline >= v2.6.1-rc.2)

    funannotate predict
      --min_pasa_complete_models INT      PASA training-set gate (default 500; 0 disables)
      --min_rnaseq_identity FLOAT         RNA-seq identity gate, median read identity % (default 0 = off)
      --no_training_single_exon           do not admit protein-supported single-exon training genes
      --min_training_models INT           minimum training models before BUSCO fallback

How the thresholds were chosen
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each choice was tested on five fungal genomes with RefSeq annotation (*Neurospora crassa* OR74A, *Aspergillus nidulans* FGSC A4, *Botrytis cinerea* B05.10, *Cryptococcus neoformans* H99, *Schizophyllum commune* H4-8). Augustus and snap were trained on half of the chromosomes, and the predictions were scored on the other half with gffcompare and BUSCO. Key results:

- On a genome with 93 complete PASA models, BUSCO training (the gate's choice) gave a holdout locus sensitivity/precision of 37.4/49.2, against 29.2/44.2 for the previous behavior.
- Single-exon training genes raised single-exon sensitivity by 5.2-7.6 points on four genomes. Multi-exon precision rose by 0.4-2.4 points, and multi-exon sensitivity changed by 0.0 to −0.9.
- The 10% RNA-seq threshold and the 500-model gate were first set from an 18-genome pilot.
- The 500-model gate was then tested on 40 RefSeq genomes (experiment B; see :ref:`assessment_pasa2.6_fun1.9`). In that test Augustus and snap were trained on half of the chromosomes, so the model counts were about half of the counts that a normal whole-genome run sees.

  - **How the gate works.** A genome with at least :code:`--min_pasa_complete_models` complete PASA models trains Augustus and snap from PASA. A genome with fewer trains them from BUSCO. A higher value therefore sends more genomes to BUSCO.
  - **With few models, PASA training often fails.** 4 of the 5 genomes with fewer than 500 complete models on half of the chromosomes lost 4.9-16.8 points of holdout locus F1 with PASA training. On the whole genome, these 4 genomes have 482-976 complete models, so a whole-genome gate at 500 would stop only 1 of them.
  - **From 500 to 1,500, the value makes almost no difference.** The genomes in this range are mixed: some do better with PASA and some with BUSCO. The mean gain over always training from PASA stays at about +1.1 points.
  - **Above 1,500, a higher value loses accuracy.** It sends genomes to BUSCO that do better with PASA. At 2,000, ten such genomes (9 of them better with PASA) move to BUSCO, and the mean gain falls to +0.9.
  - These points describe counts on half of the chromosomes. Whether 500 is also the right value for whole-genome counts is not yet measured (experiment C). The default stays at 500 until then. If a genome has fewer than about 1,000 complete PASA models, compare PASA and BUSCO training before you rely on either.
- Read identity of the RNA-seq did not predict whether PASA or BUSCO training was better in the same 40 genomes, so :code:`--min_rnaseq_identity` is off by default. RNA-seq evidence (hints and PASA models in EVM) raised holdout locus F1 in all 12 genomes tested, also for reads from another strain, so no gate removes it.

Full tables are in the methods document linked above.
