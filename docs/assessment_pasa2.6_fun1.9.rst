
.. _assessment_pasa2.6_fun1.9:

Training and evidence assessment (v1.9.0)
=========================================

This page records how well gene-prediction training and transcript evidence perform in
**funannotate v1.9.0-rc.3** with **PASApipeline v2.6.1-rc.2**, measured against RefSeq annotations in
September 2026. The data, figures (PNG and PDF), analysis scripts and a Markdown version of this page
are in ``docs/assessment_pasa2.6_fun1.9/``.

:Versions assessed: funannotate ``v1.9.0-rc.3`` (commit ``cd1b5ee``); PASApipeline ``v2.6.1-rc.2``
                    (hyphaltip/PASApipeline, commit ``1044957``)
:Baseline:          funannotate ``41a2fd7`` and the ``funannotate-1.9.0-rc.1`` image
:As of:             2026-09-26
:Related pages:     :ref:`training`, :ref:`evidence`; the selection methods are in
                    ``docs/training_data_selection_methods.md``

Summary
-------

1. Training Augustus and SNAP from RNA-seq-derived PASA models is rarely better than training from
   BUSCO genes by more than about 1 point. It is clearly worse when fewer than about 300 complete PASA
   models are available.
2. RNA-seq is valuable mainly as evidence for the final gene models, through PASA models in EVM
   (PASA weight 6).
3. funannotate therefore uses RNA-seq in both roles:

   - it trains from PASA only when enough complete models exist (``--min_pasa_complete_models``,
     default 500), and falls back to BUSCO otherwise;
   - it always passes the PASA models to EVM as evidence.

How it was measured
-------------------

- **Genomes.** Five fungal genomes with RefSeq annotation:

  - *Neurospora crassa* OR74A (divergent RNA-seq, 96.7% read identity);
  - *Aspergillus nidulans* FGSC A4 (same-strain RNA-seq);
  - *Botrytis cinerea* B05.10 (same-strain RNA-seq);
  - *Cryptococcus neoformans* H99 (same-strain RNA-seq);
  - *Schizophyllum commune* H4-8 (divergent RNA-seq).
- **Held-out chromosomes.** Chromosomes of at least 1 Mb were assigned alternately to training and
  holdout sets. Augustus and SNAP were trained on the training chromosomes only; only the holdout
  chromosomes were scored.
- **Scoring.** gffcompare (v0.12.10) at CDS level, plus proteome BUSCO (fungi_odb10).
- **Evidence settings.**

  - "Fixed" gives every arm the same EVM evidence, so only the training differs.
  - "Own" gives each arm its own PASA models as evidence.

PASA-trained versus BUSCO-trained predictors
--------------------------------------------

.. figure:: assessment_pasa2.6_fun1.9/figures/fig2_training_source_by_genome.png
   :width: 100%

   Holdout locus F1 with the same fixed evidence. Filled marker: PASA-trained; open marker:
   BUSCO-trained. PDF: ``assessment_pasa2.6_fun1.9/figures/fig2_training_source_by_genome.pdf``.

.. figure:: assessment_pasa2.6_fun1.9/figures/fig1_titration_pasa_vs_busco.png
   :width: 100%

   Experiment A (final). Complete PASA training models were subsampled within each genome and
   compared with BUSCO training on the same evidence. Whiskers: 95% bootstrap over subsample draws.

.. list-table:: Holdout locus F1, PASA-trained minus BUSCO-trained (points). PASA: mean over 3-5 subsets; BUSCO: mean of 3 repeat runs
   :header-rows: 1
   :widths: 24 19 19 19 19

   * - Complete training models
     - *A. nidulans*
     - *B. cinerea*
     - *C. neoformans* H99
     - *N. crassa*
   * - 50
     - −2.5
     - −3.0
     - −2.2
     - −4.6
   * - 100
     - −1.4
     - −1.5
     - −0.9
     - −2.6
   * - 200
     - −0.9
     - +0.1
     - +0.1
     - −1.3
   * - 300
     - −0.5
     - −0.4
     - +0.6
     - −1.5
   * - 500
     - −0.3
     - +0.7
     - +0.3
     - −0.5
   * - 750
     - −0.4
     - +0.2
     - +0.8
     - −1.0
   * - 1000
     - −0.2
     - +0.9
     - +1.1
     - 0.0
   * - 2000
     - +0.2
     - +1.7
     - +1.1
     - (pool 1,999)

- At 100 complete models and below, BUSCO training wins in all four genomes (by 0.9-4.6 points).
- At 500 models and above, the difference is between −1.0 and +1.7 points.
- Under the pre-agreed conservative rule (the lower 95% bound of the difference is at least 0 at that N
  and at every larger N tested), the locus-level crossover is 300 for *C. neoformans*, 1,000 for
  *B. cinerea* and 2,000 for *A. nidulans*. It is not reached for *N. crassa* up to 1,000 (the largest
  N possible). Exon and intron-chain levels give the same values, except *C. neoformans* exon (750).
- The most conservative single threshold over these four genomes is 2,000, where the gain over BUSCO
  is +0.2 to +1.7 points. The funannotate default gate is still 500. Whether to change it is a
  decision for after experiment B (about 40 RefSeq genomes).
- With only 93 complete models (*S. commune*), BUSCO training is ahead by +8.2 locus sensitivity and
  +5.0 precision.

RNA-seq as evidence
-------------------

.. figure:: assessment_pasa2.6_fun1.9/figures/fig3_evidence_vs_training_effect.png
   :width: 100%

   Change in holdout locus sensitivity and precision from the alignment fixes, split into the effect of
   the training set alone, the evidence alone, and both.

- On *N. crassa*, the fixes changed locus Sn / Pr as follows:

  - training set alone: +0.3 / +0.4;
  - evidence alone: +1.2 / +0.6;
  - both: +2.0 / +1.8.
- On *B. cinerea* the pattern is the same: training alone −0.2 / −0.3, evidence alone +0.5 / 0.0.
- EVM weights in these runs: PASA 6, Augustus HiQ 2, and 1 each for Augustus, GeneMark, SNAP, proteins
  and transcripts.
- Not measured: the separate contribution of RNA-seq-derived Augustus hints, and GeneMark with
  RNA-seq (these runs used GeneMark-ES).

Transcript alignment quality
----------------------------

.. figure:: assessment_pasa2.6_fun1.9/figures/fig4_minimap2_fix_validation.png
   :width: 80%

   Spliced minimap2 alignments that pass PASA validation, before and after the conversion fix in
   ``library.bam2gff3`` (coordinates from the CIGAR instead of the 2018 ``cs``-tag walk).

.. figure:: assessment_pasa2.6_fun1.9/figures/fig5_intron_accuracy_by_aligner.png
   :width: 100%

   Introns exactly matching a RefSeq intron, by aligner and alignment identity (*N. crassa*,
   divergent RNA-seq).

.. list-table:: Effect of the fixes on PASA training models (exact RefSeq CDS chains)
   :header-rows: 1
   :widths: 40 20 20 20

   * - Change
     - *N. crassa*
     - *B. cinerea*
     - *A. nidulans*
   * - minimap2 conversion fix (R1)
     - 1,481 → 2,096
     - 5,353 → 5,644
     - 3,171 → 3,220
   * - R1 + PASA ``--UNSPLICED_JOIN_SPLICED`` and ``--ONE_ALIGNMENT_PER_CDNA``
     - 2,096 → 2,213
     - 5,644 → 5,649
     - 3,220 → 3,221

- Fixed minimap2 placed 89.5% of introns exactly as RefSeq, against 61.0% for gmap and 60.8% for blat.
- Real minimap2 splice-placement errors were about 1.2% of introns, so no polishing step is needed.
- Relaxing PASA validation (identity 90%, no splice-boundary rule) added training models but did not
  improve predictions.

Single-exon training genes
--------------------------

.. figure:: assessment_pasa2.6_fun1.9/figures/fig8_single_exon_training_effect.png
   :width: 100%

   Change in exact CDS-chain accuracy when protein-supported single-exon genes are admitted to the
   training set (now the default; ``--no_training_single_exon`` turns it off).

Runtime
-------

- The minimap2 conversion fix increased ``funannotate train`` wall time by 12-20% (*N. crassa*
  1,266 → 1,520 s; *A. nidulans* 905 → 1,015 s), because PASA validates and assembles more alignments.
- gmap as PASA's aligner took 67% longer than blat, with no accuracy gain.

Production data (BFD, 8,007 PASA-trained genomes)
-------------------------------------------------

.. figure:: assessment_pasa2.6_fun1.9/figures/fig6_f1_production_timeline.png
   :width: 100%

   PASA training models with CDS length not divisible by 3, by file date. A PASA defect (duplicate
   GFF3 rows) caused frame errors in 5,599 genomes between 2026-07-06 and 2026-09-24. It is fixed in
   PASApipeline v2.6.1-rc.1 and later.

.. figure:: assessment_pasa2.6_fun1.9/figures/fig7_production_rnaseq_identity.png
   :width: 100%

   Median transcript-to-genome identity per genome against the read-identity categories (same strain
   ≥ 99%; divergent 90-99%; below 90% trains from BUSCO).

Open items
----------

- Experiment A: add *C. neoformans* H99, and repeat the BUSCO comparator on the same code snapshot.
- Experiment B: about 40 RefSeq genomes, to test whether the crossover holds across species.
- A hints-on versus hints-off arm; GeneMark-ET/EP; refitting EVM weights after the fixes.

Files
-----

All paths are relative to ``docs/assessment_pasa2.6_fun1.9/``.

- ``README.md``: this assessment in Markdown, with the full tables.
- ``evidence_and_alignment_methods.md``: methods and results draft for the evidence and alignment work.
- ``DECISIONS_snapshot_2026-09-26.md``: the decision log.
- ``data/``: the result tables every number comes from.
- ``figures/``: PNG and PDF figures, and ``make_figures.py``, which regenerates them from ``data/``.
- ``scripts/``: the analysis and scoring scripts.
