# Assessment: training from PASA vs BUSCO, and the value of RNA-seq evidence (funannotate v1.9.0-rc.3, PASA v2.6.1-rc.2)

This folder records how well funannotate's gene-prediction training and evidence perform after the 2026-09 review of PASApipeline and funannotate. It is written to be read without the original conversation.

- **Versions assessed:** funannotate **v1.9.0-rc.3** (tag `v1.9.0-rc.3`, commit `cd1b5ee`; image `funannotate-1.9.0-rc.3.sif`) and PASApipeline **v2.6.1-rc.2** (hyphaltip/PASApipeline tag `v2.6.1-rc.2`, commit `1044957`). Comparisons against earlier behavior use funannotate `41a2fd7` and the rc.1 image.
- **As of:** 2026-09-26.
- **Sphinx page:** `docs/assessment_pasa2.6_fun1.9.rst`.
- **Prepared by:** two Claude Opus 5.5 sessions working for J. Stajich. "REVIEW" did the PASA/evidence side; "SELECT" did the training-data selection side.
- **Measurement:** every number is from a result file in `data/` or from the BFD project directory named at the end.

## Answer in brief

1. **Training Augustus/SNAP from RNA-seq-derived PASA ORFs is rarely better than training from BUSCO genes by more than about 1 point.** It is clearly worse when few complete models exist.
2. **RNA-seq is valuable mainly as evidence for the final gene models,** through PASA models in EVM (PASA weight 6).
3. **funannotate now uses both roles:**
   - it trains from PASA only when enough complete models exist, and falls back to BUSCO otherwise (the "PASA gate", default 500 complete models);
   - it always uses the RNA-seq-derived PASA models as EVM evidence.

## 1. PASA-trained vs BUSCO-trained predictors

**Method.** Chromosomes of at least 1 Mb were split alternately into training and holdout sets. Augustus/SNAP were trained on the training chromosomes; the whole genome was predicted; only the holdout chromosomes were scored against RefSeq with gffcompare (CDS level).
- "fixed evidence": every arm gets the same EVM evidence, so only the training differs.
- "own evidence": each arm uses its own PASA models as evidence, which is closest to production.

**Holdout locus sensitivity / precision (%)** [`data/scorecard.tsv`; SELECT Table M2/M5/M6]:

| Genome (RNA-seq) | PASA-trained, fixed evidence | BUSCO-trained, fixed evidence | Note |
|---|---|---|---|
| N. crassa OR74A (divergent, 96.7% read identity) | 64.3 / 72.6 | 66.0 / 74.1 | BUSCO ahead |
| N. crassa, PASA with fixes R1+R2 + single-exon training | 66.4 / 74.9 | 66.0 / 74.1 | PASA ahead by <1 |
| A. nidulans FGSC A4 (same strain) | 54.6 / 55.8 | 54.8 / 56.8 | about equal |
| B. cinerea B05.10 (same strain) | 83.3 / 82.2 | 82.0 / 82.0 | PASA ahead |
| C. neoformans H99 (same strain) | 79.9 / 80.8 | 78.1 / 80.6 | PASA ahead |
| S. commune H4-8 (divergent, 93 complete PASA models) | 24.9 / 36.0 | 37.4 / 49.2 | BUSCO far ahead; the gate chooses BUSCO here |

**Dependence on the number of complete PASA training models (experiment A, interim).** This is holdout locus F1, PASA-trained minus BUSCO-trained, fixed evidence, mean over 3-5 random subsets [`data/titration_analysis_locus.tsv`]:

| Complete training models | A. nidulans | B. cinerea | N. crassa |
|---|---|---|---|
| 50 | −2.7 | −3.0 | −4.6 |
| 100 | −1.5 | −1.5 | −2.6 |
| 200 | −1.0 | +0.1 | −1.2 |
| 300 | −0.6 | −0.4 | −1.5 |
| 500 | −0.4 | +0.7 | −0.5 |
| 750 | −0.6 | +0.2 | −1.0 |
| 1000 | −0.3 | +0.9 | 0.0 |
| 2000 | +0.1 | +1.7 | pool is 1,999 |

![Figure 1](figures/fig1_titration_pasa_vs_busco.png)

**Reading:**
- Below about 300 complete models, BUSCO training wins by 1-5 points.
- Above 500 the two are within about ±1 point, except B. cinerea (+1.7 at 2,000).
- Under the pre-agreed conservative rule (the lower 95% bound of the difference must be at least 0), the crossover is 1,000 for B. cinerea and is not reached for A. nidulans or N. crassa.
- The gate threshold (500) is not too high. Its final value waits for the items in section 4.

## 2. RNA-seq as evidence

![Figure 3](figures/fig3_evidence_vs_training_effect.png)

**N. crassa, holdout locus Sn / Pr change from the alignment fixes R1 + R2** [`data/scorecard.tsv`]:

| Effect | Change |
|---|---|
| Training set only (fixed − fixed) | +0.3 / +0.4 |
| Evidence only, within R1 + R2 (own − fixed) | +1.2 / +0.6 |
| Both (own − own) | +2.0 / +1.8 |

- B. cinerea: training only −0.2 / −0.3; evidence only +0.5 / 0.0; both +0.7 / +0.5.
- EVM weights in these runs: PASA 6, Augustus HiQ 2, and 1 each for Augustus, GeneMark, SNAP, proteins and transcripts.
- **Not measured:** the contribution of RNA-seq-derived Augustus hints (no hints-off arm was run).
- **Not measured:** GeneMark with RNA-seq. These runs used GeneMark-ES (ab initio).

## 3. Performance of the fixes

**Accuracy of PASA evidence (training models against RefSeq; exact CDS-chain matches):**

| Change | N. crassa | B. cinerea | A. nidulans |
|---|---|---|---|
| R1: minimap2 → GFF3 conversion from the CIGAR (funannotate PR #1210) | 1,481 → 2,096 (+42%) | 5,353 → 5,644 (+5.4%) | 3,171 → 3,220 (+1.5%) |
| Spliced minimap2 alignments valid in PASA, before → after R1 | 0.8% → 64% | 12% → 88% | 75% → 95% |
| R1 + R2/F4 PASA flags, against R1 alone | 2,096 → 2,213 | 5,644 → 5,649 | 3,220 → 3,221 |

- **Splice-site accuracy on divergent reads** (introns exactly matching RefSeq): fixed minimap2 89.5%, gmap 61.0%, blat 60.8%.

![Figure 4](figures/fig4_minimap2_fix_validation.png)
![Figure 5](figures/fig5_intron_accuracy_by_aligner.png)
![Figure 8](figures/fig8_single_exon_training_effect.png)
- **Single-exon training genes** (R6 option b, now on by default): single-exon Sn +5.2 to +7.6 points; multi-exon Pr +0.4 to +2.4; multi-exon Sn 0 to −0.9 [`data/single_exon_scores.tsv`].
- **Training-set selection (R3/R5)** raised the share of training models exactly matching RefSeq [`data/benchmark.tsv`]:
  - N. crassa 38.9 → 78.7%;
  - A. nidulans 50.1 → 62.1%;
  - B. cinerea 75.5 → 89.4%.

**Runtime:**
- R1 increased funannotate train wall time by 12-20% (N. crassa 1,266 → 1,520 s; A. nidulans 905 → 1,015 s), because PASA has more valid alignments to process.
- gmap as PASA's aligner took 67% longer than blat, with no accuracy gain.

**Production impact of the PASA duplicate-row bug (F1)** [`data/production_f1_scan.tsv.gz`]:
- 5,599 of 8,007 PASA-trained BFD genomes (dated 2026-07-06 to 09-24) have frame-broken training models.
- Genomes dated earlier are clean.
- They are candidates for a rerun with v1.9.0-rc.3.

![Figure 6](figures/fig6_f1_production_timeline.png)
![Figure 7](figures/fig7_production_rnaseq_identity.png)

## 4. Open items

- Experiment A: add C. neoformans H99; repeat the BUSCO comparator 3 times on the same code snapshot (the interim comparator used an older snapshot).
- Experiment B: about 40 RefSeq BFD genomes, stratified by complete-model count and read identity, to test whether the crossover holds across species.
- A hints-on vs hints-off arm; GeneMark-ET/EP; refitting EVM weights after the fixes.

## Files in this folder

- `evidence_and_alignment_methods.md`: REVIEW's methods and results draft (evidence, alignment, F1, gate calibration).
- `../training_data_selection_methods.md`: SELECT's methods draft (gates, R3/R5, single-exon training, Tables M1-M8).
- `DECISIONS_snapshot_2026-09-26.md`: a snapshot of the shared decision log (D01-D97).
- `data/`:
  - `scorecard.tsv`: predict arms, holdout gffcompare scores.
  - `single_exon_scores.tsv`: single-exon vs multi-exon exact-match scores.
  - `benchmark.tsv`, `rank_benchmark.tsv`: training models and ranking variants against RefSeq.
  - `titration_scores.tsv`, `titration_analysis_{locus,exon,intron_chain}.tsv`: experiment A, interim.
  - `production_f1_scan.tsv.gz`, `production_identity.tsv.gz`: production scans (8,007 genomes).
- `figures/`: PNG (embedded above) and PDF versions of Figures 1-8. `figures/make_figures.py` regenerates all of them from `data/` (`/usr/bin/python3.12 docs/assessment_pasa2.6_fun1.9/figures/make_figures.py`; needs matplotlib).
- `scripts/`: the analysis scripts that produced the tables (`training_set_vs_refseq.py`, `titration_analysis.py`, `intron_discordance.py`, `production_f1_scan.py`, `production_identity.py`, `predict_scorer.py`, `diversity.py`, and the R13 scripts).
- Full code review of the PASA fork: hyphaltip/PASApipeline, `CODE_REVIEW_20260925.md` on branch `rust_optimize`.
- Working data (UCR HPCC): `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/`.
