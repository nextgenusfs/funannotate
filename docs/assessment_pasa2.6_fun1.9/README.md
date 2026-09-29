# Assessment: training from PASA vs BUSCO, and the value of RNA-seq evidence (funannotate v1.9.0-rc.3, PASA v2.6.1-rc.2)

This folder records how well funannotate's gene-prediction training and evidence perform after the 2026-09 review of PASApipeline and funannotate. It is written to be read without the original conversation.

- **Versions assessed:** funannotate **v1.9.0-rc.3** (tag `v1.9.0-rc.3`, commit `cd1b5ee`; image `funannotate-1.9.0-rc.3.sif`) and PASApipeline **v2.6.1-rc.2** (hyphaltip/PASApipeline tag `v2.6.1-rc.2`, commit `1044957`). Comparisons against earlier behavior use funannotate `41a2fd7` and the rc.1 image.
- **As of:** 2026-09-26; sections 4 and 5 added 2026-09-28.
- **Sphinx page:** `docs/assessment_pasa2.6_fun1.9.rst`.
- **Prepared by:** two Claude Opus 5.5 sessions working for J. Stajich. "REVIEW" did the PASA/evidence side; "SELECT" did the training-data selection side.
- **Measurement:** every number is from a result file in `data/` or from the BFD project directory named at the end.

## Answer in brief

1. **Training Augustus/SNAP from RNA-seq-derived PASA ORFs is rarely better than training from BUSCO genes by more than about 1 point.** It is clearly worse when few complete models exist.
2. **RNA-seq is valuable mainly as evidence for the final gene models,** through PASA models in EVM (PASA weight 6).
3. **funannotate now uses both roles:**
   - it trains from PASA only when enough complete models exist, and falls back to BUSCO otherwise (the "PASA gate", default 500 complete models);
   - it always uses the RNA-seq-derived PASA models as EVM evidence.
4. **Read identity of the RNA-seq does not predict the training outcome** (experiment B, 40 genomes; section 4). RNA-seq evidence helped in all 12 genomes tested, down to 94% identity. So a low-identity gate may change only the training source, never the evidence. The identity gate is off by default.
5. **The PASA gate of 500 complete models works in experiment B, but it is not yet calibrated for production** (section 5). Experiment B counted models on the training chromosomes (about half the genome). Production counts on the whole genome, which gives about 2 times more (1.5-15 times). On whole-genome counts, 500 would catch only 1 of the 5 genomes that lost badly with PASA training. Experiment C tests the whole-genome case.

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

**Dependence on the number of complete PASA training models (experiment A, final).** This is holdout locus F1, PASA-trained minus BUSCO-trained, fixed evidence. Each PASA value is the mean over 3-5 random subsets. The BUSCO comparator is the mean of 3 repeat runs, which differ by at most 0.1 point. Brackets: 95% bootstrap interval over subsets [`data/titration_analysis_locus.tsv`]:

| Complete training models | A. nidulans | B. cinerea | C. neoformans H99 | N. crassa |
|---|---|---|---|---|
| 50 | −2.5 [−2.8, −2.2] | −3.0 [−3.3, −2.6] | −2.2 [−3.0, −1.6] | −4.6 [−5.4, −3.8] |
| 100 | −1.4 [−1.6, −1.2] | −1.5 [−1.8, −1.4] | −0.9 [−1.4, −0.5] | −2.6 [−3.5, −1.9] |
| 200 | −0.9 [−1.0, −0.8] | +0.1 [0.0, +0.3] | +0.1 [−0.2, +0.5] | −1.3 [−1.5, −0.9] |
| 300 | −0.5 [−0.7, −0.3] | −0.4 [−1.0, +0.3] | +0.6 [+0.4, +0.8] | −1.5 [−2.1, −0.9] |
| 500 | −0.3 [−0.6, −0.1] | +0.7 [+0.2, +1.2] | +0.3 [+0.2, +0.4] | −0.5 [−1.0, 0.0] |
| 750 | −0.4 [−0.7, −0.2] | +0.2 [−0.8, +0.8] | +0.8 [+0.3, +1.1] | −1.0 [−1.5, −0.7] |
| 1000 | −0.2 [−0.3, 0.0] | +0.9 [+0.7, +1.2] | +1.1 [+1.0, +1.2] | 0.0 [−0.3, +0.4] |
| 2000 | +0.2 [+0.1, +0.4] | +1.7 [+1.7, +1.8] | +1.1 [+0.8, +1.4] | not tested (pool is 1,999) |

![Figure 1](figures/fig1_titration_pasa_vs_busco.png)

**Reading:**
- At 100 complete models and below, BUSCO training wins in all four genomes (by 0.9-4.6 points).
- At 500 models and above, the difference is between −1.0 and +1.7 points.
- Under the pre-agreed conservative rule (the lower 95% bound of the difference is at least 0 at that N and at every larger N tested), the crossover is:
  - C. neoformans H99: 300 (locus and intron chain), 750 (exon);
  - B. cinerea: 1,000 at all three levels;
  - A. nidulans: 2,000 at all three levels;
  - N. crassa: not reached up to 1,000 (the largest N possible).
- The most conservative single threshold over these four genomes is 2,000. At 2,000 the gain over BUSCO is small (+0.2 to +1.7 points). Between 500 and 2,000, the loss is at most 1.0 point (N. crassa, 750).
- Four genomes do not fix a threshold for all fungi. Experiment B (section 5) tested the crossover across 40 genomes. Its counts are on the training chromosomes, so its result does not transfer directly to the whole-genome counts that production uses (section 5.3).
- Exon and intron-chain results: `data/titration_analysis_exon.tsv`, `data/titration_analysis_intron_chain.tsv`.

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

## 4. RNA-seq read identity as a training gate (experiment B, 2026-09-28)

### 4.1 The question

- **External user report (not reproduced here).** A user mapped *P. brasiliensis* RNA-seq to three *Paracoccidioides* genomes. They used the method of the train gate: the first 200,000 R1 reads, `minimap2 -ax splice:sr --secondary=no`, primary alignments with MAPQ ≥ 1.

  | Genome | Reads mapped | Median read identity |
  |---|---|---|
  | Pb18 (same species) | 96.8% | 99.4% |
  | *P. lutzii* 7730 | 94.4% | 95.7% |
  | *P. lobogeorgii* | 79.8% | 93.6% |

  - The *P. lobogeorgii* reads passed the 10% map-rate gate.
  - Training from them left 7.0% of intact loci without a gene model. BUSCO training left 1.4%.
  - The user proposed an identity floor of about 97-98% in the same gate. The text says this value comes from "your titration".
- **Checks against the code and the data:**
  - Correct: the train gate measured only the map rate. minimap2 aligns divergent reads with high MAPQ, so the map rate cannot detect reads from a related species.
  - Not correct: experiment A titrated the number of complete training models, not read identity. No identity threshold had been tested before this section.
  - A design risk: the map-rate gate stops train (exit 3). If the identity floor also stopped train, a genome would lose all RNA-seq evidence (BAM hints, PASA models, transcripts), not only PASA training.

### 4.2 What was changed in funannotate (working tree, not committed as of 2026-09-28)

- **train** measures identity with the same 200,000 reads that the map-rate gate uses:
  - identity per read = 1 − NM / (aligned M/=/X bases + inserted bases), the D46 method; introns, deletions and clips are not in the denominator;
  - the median and the 10th percentile go to `logfiles/train_rnaseq_gate.tsv` (new columns `median_identity_pct`, `p10_identity_pct`, `identity_reads`) and to `training_decisions.tsv` (stage `rnaseq_identity`);
  - train never stops because of identity;
  - a copy of the report goes to `training/funannotate_train.rnaseq_gate.tsv`, because the BFD pipeline keeps only the `training/` folder.
- **predict** has a new option, `--min_rnaseq_identity` (0 = off; default changed from 95 to 0 on 2026-09-28 after this analysis):
  - It is checked at the same point as the PASA gate.
  - Below the threshold, predictors that would train from PASA train from BUSCO.
  - The RNA-seq BAM, PASA models and transcripts are still used as hints and EVM evidence.
  - It is skipped, with a recorded reason, when train did not record identity (older train folders, or `--min_rnaseq_map_rate 0`).
  - Stage `rnaseq_identity_gate`; the `training_mode_switch` record names the gate that failed.
- Tests: `tests/test_training_gates.py`, 49 tests (identity from CIGAR/NM, report parsing, the gate decision, minimap2 on synthetic reads with 5 substitutions per 150 bp → 96.67%). Full suite: 162 tests pass. No end-to-end train → predict run was done with this code.
- The first default, 95, was set before the measurement below. The user's words were "default threshold is probably 95% but we need to go through it empirically". After the measurement the user decided: "yes change default to 0".

### 4.3 Method

- **Genomes:** the 40 experiment B RefSeq genomes with PASA, BUSCO and holdout scores. C. siamense Cg363 and D. hansenii CBS767 have no reads (their RNA-seq failed the map-rate gate in the pilot).
- **Reads:** `training/normalize/left.norm.fq.gz` from the rc.3 pilot. This is the file that train's gate samples.
- **Measurement:** the funannotate functions above (`lib.sample_read_map_rate`, `lib.summarize_identity`), 200,000 reads, genome `experiment_B/<genome>/genome.fa`. One SLURM job (29190649) ran all 41 read sets in 3 minutes.
- **Check:** N. crassa gave 96.67% median identity and a 95.79% map rate. D46 gave 96.7% and 95.8% on raw reads.
- **Outcome:** holdout locus F1 (2·Sn·Pr / (Sn + Pr)), step B of experiment B:
  - PASA − BUSCO: `pasa.B` minus the mean of `busco_r1..r3.B`, same fixed evidence;
  - evidence effect: `busco` mean minus `norna.B`. Both are BUSCO-trained; `norna` has no RNA-seq evidence at any step (no BAM hints, no PASA models, no transcripts).
- **Complete models:** the `pasa_gate` record of `pasa.A` (complete-ORF PASA models on the train chromosomes).
- **Statistics:** Spearman rank correlation; genome-level bootstrap (2,000 resamples) for mean differences.
- Files: `data/expB_read_identity.tsv`, `data/expB_identity_vs_training.tsv`, `data/expB_identity_analysis.txt`; scripts `scripts/expB_run_identity.sh`, `scripts/expB_measure_identity.py`, `scripts/expB_analyze_identity.py`.

### 4.4 Results

**Read identity of the experiment B RNA-seq** [`data/expB_read_identity.tsv`]:
- 25 of 41 read sets have a median identity of 100%; 34 are at 98.5% or higher.
- Below 95%: P. antarcticum IBT 31339 (92.72%), S. commune H4-8 (94.04%), P. hubeiensis SY62 (94.04%).
- 95-97%: B. deweyae B1 (95.30%), N. crassa OR74A (96.67%), T. versicolor (97.00%).

**What predicts PASA − BUSCO holdout locus F1** [`data/expB_identity_analysis.txt`]:

| Predictor | Spearman ρ | Genomes |
|---|---|---|
| Median read identity | 0.18 | 40 |
| Median read identity, ≥ 500 complete models | 0.12 | 35 |
| Map rate, ≥ 500 complete models | 0.09 | 35 |
| log(complete PASA models) | 0.42 | 40 |

- Identity and log(complete models) are not correlated (ρ = 0.02). They are separate properties of a data set.

![Figure 9](figures/fig9_identity_vs_training_outcome.png)

**The genomes below 95% identity:**

| Genome | Identity | Map rate | Complete models | PASA − BUSCO | Gate that acts |
|---|---|---|---|---|---|
| P. antarcticum IBT 31339 | 92.72% | 63.4% | 1,730 | **+0.77** | identity gate only (would lose 0.77) |
| S. commune H4-8 | 94.04% | 68.3% | 354 | −4.90 | PASA gate already |
| P. hubeiensis SY62 | 94.04% | 61.2% | 141 | −11.30 | PASA gate already |

**Threshold sweep**, genomes with ≥ 500 complete models (n = 35), which the PASA gate does not already send to BUSCO:

| Threshold | Genomes below | BUSCO better, below | Mean PASA − BUSCO below [95% CI] | Mean above [95% CI] |
|---|---|---|---|---|
| 95.0% | 1 | 0 | +0.77 | +0.68 [−0.03, +1.33] |
| 97.0% | 3 | 1 | +0.05 [−0.66, +0.77] | +0.74 [−0.04, +1.39] |
| 98.0% | 4 | 2 | −0.05 [−0.51, +0.48] | +0.78 [+0.03, +1.48] |
| 99.0% | 11 | 3 | +0.32 [−1.27, +1.56] | +0.85 [+0.18, +1.56] |

- No threshold from 93% to 99.5% gives a net gain in mean locus F1 over the 35 genomes (range −0.18 to +0.01 points).

**Large losses from PASA training, and which gate catches them:**
- Fewer than 500 complete models (PASA gate catches them): A. niger CBS 101883 −16.79 (309 models), E. xenobiotica −14.20 (312), P. hubeiensis −11.30 (141), S. commune −4.90 (354). Exception: A. thermomutatus, 436 models, +2.85.
- Not caught by any gate:
  - M. bicuspidata: −6.56, with 98.67% identity, 1,456 complete models and a 33.9% map rate;
  - three yeasts with 100% identity: N. castellii −2.44, S. cerevisiae −2.04, K. capsulata −1.91.

**RNA-seq as evidence when training is from BUSCO** (`busco` − `norna`, holdout locus F1):

| Genome | Identity | BUSCO + RNA-seq evidence | BUSCO, no RNA-seq | Gain |
|---|---|---|---|---|
| S. commune H4-8 | 94.04% | 42.92 | 40.35 | +2.57 |
| N. crassa OR74A | 96.67% | 70.58 | 66.02 | +4.56 |
| T. versicolor FP-101664 | 97.00% | 46.94 | 42.01 | +4.93 |
| A. bisporus JB137-s8 | 98.55% | 45.00 | 31.47 | +13.53 |
| P. chlamydosporia 170 | 98.66% | 66.83 | 66.44 | +0.39 |
| C. militaris CM01 | 99.33% | 52.28 | 51.09 | +1.19 |
| A. nidulans FGSC A4 | 100% | 55.67 | 51.79 | +3.88 |
| B. cinerea B05.10 | 100% | 82.18 | 69.32 | +12.86 |
| C. neoformans H99 | 100% | 79.10 | 55.12 | +23.98 |
| E. mesophila CBS 40295 | 100% | 74.79 | 65.61 | +9.18 |
| K. lactis NRRL Y-1140 | 100% | 90.62 | 85.15 | +5.48 |
| R. microsporus ATCC 52813 | 100% | 49.10 | 41.51 | +7.59 |

- Mean gain +7.51 points [95% CI +4.48, +11.52]; positive in 12 of 12 genomes.
- Below 99% identity: mean +5.20, positive in 5 of 5, including S. commune at 94.0%.

### 4.5 Conclusions

1. **Read identity does not predict whether PASA or BUSCO training is better** in these 40 genomes. *Measured.* The complete-model count does, and the PASA gate already uses it.
2. **A 95% identity default is not supported by this data.** *Measured, small n.* It would change the training source for one genome (P. antarcticum), and PASA training was 0.77 points better there. The user-reported *P. lobogeorgii* case (93.6%, BUSCO better) points the other way. Two cases of similar identity with opposite outcomes do not support one threshold.
3. **RNA-seq evidence helps even when the reads are divergent.** *Measured, 12 genomes.* The gain is positive in every genome tested, down to 94.0% identity. So low identity must not remove RNA-seq evidence. The implemented gate therefore changes only the training source.
4. **Some losses are not explained by identity, map rate or complete-model count** (M. bicuspidata, three same-strain yeasts). *Measured; cause not investigated.*

### 4.6 Decision

- Keep the identity measurement and the `--min_rnaseq_identity` option, so users with a case like *P. lobogeorgii* can set a floor.
- **Decided (user, 2026-09-28): the default is 0 (report only).** train still measures and records identity for every run.
- To calibrate an identity gate, more genomes below 95% with ≥ 500 complete models are needed. The experiment B set has one.

### 4.7 Limitations

- Each genome has one PASA run. The BUSCO spread over 3 repeats is at most 0.31 points; the PASA run-to-run spread is not measured. Differences of about 1 point or less are within noise.
- Only 3 genomes are below 95% identity, and only 1 of them has ≥ 500 complete models.
- The reads are the normalized reads, which is what train's gate samples. The user report used raw R1 reads. For N. crassa the two gave the same median (96.67% here against 96.7% in D46).
- The user's outcome measure (intact loci without a model) is not the holdout locus F1 used here.
- The `norna` arm has no RNA-seq at any step, so the evidence gain includes both Augustus hints and EVM evidence. The two are not separated.
- The complete-model threshold across experiment B (D81) is in section 5.

## 5. The complete-model gate across 40 genomes (experiment B, D81, 2026-09-28)

### 5.1 Method

- **Genomes and arms:** the 40 experiment B genomes of section 4. PASA − BUSCO holdout F1 at locus, exon and intron-chain level, with BUSCO as the mean of 3 repeats.
- **Gate variable:** complete-ORF PASA models on the training chromosomes, counted with the PASA gate's own function (`lib.count_complete_orf_models`). In step A, the genome being trained is the training chromosomes, so this is the count the gate sees in experiment B. The counts match the `pasa_gate` records of `pasa.A` (for example N. crassa 2,380). **Correction (2026-09-28):** this is not the count production sees. Production trains on the whole genome and counts complete models there. The whole-genome count is 1.5-14.9 times the training-chromosome count (median 2.05; A. thermomutatus is 14.9 because its split fell back to 200 kb contigs) [`data/expB_whole_genome_gate_sweep.txt`]. An earlier version of this section said the two counts were the same kind. That was wrong.
- **Second variable:** final PASA training models after selection (`select_final` record of `pasa.A`).
- **Policy tested:** "train from PASA if complete ≥ T, else BUSCO". For each T: mean F1 gain over always-PASA and over always-BUSCO, with a genome-level bootstrap 95% CI (2,000 resamples).
- **Conservative rule (as in experiment A):** T* is the smallest T at which the lower 95% bound of the mean PASA − BUSCO among genomes with ≥ T models is at least 0, at T and at every larger T with at least 5 genomes.
- Files: `data/expB_complete_models.tsv`, `data/expB_complete_threshold_by_genome.tsv`, `data/expB_complete_threshold_analysis.txt`, `data/expB_final_models_sweep.txt`, `data/expB_combined_gate_sweep.txt`; scripts `scripts/expB_count_complete.py`, `scripts/expB_analyze_complete_threshold.py`, `scripts/expB_gate_sweeps.py`.

### 5.2 Results

**Mean holdout locus F1 over 40 genomes:** always-PASA 66.85, always-BUSCO 67.36.

**Policy gain by threshold, locus F1** [`data/expB_complete_threshold_analysis.txt`]:

| T (complete models) | Genomes ≥ T | PASA better, ≥ T | Mean PASA − BUSCO, ≥ T [95% CI] | Gain vs always-PASA [95% CI] | Gain vs always-BUSCO [95% CI] |
|---|---|---|---|---|---|
| 300 | 39 | 28 | −0.23 [−1.67, +0.93] | +0.28 [0.00, +0.85] | −0.23 [−1.64, +0.90] |
| 400 | 36 | 28 | +0.74 [+0.10, +1.35] | +1.18 [+0.12, +2.52] | +0.67 [+0.08, +1.21] |
| **500** | **35** | **27** | **+0.68 [+0.02, +1.31]** | **+1.11 [+0.05, +2.48]** | **+0.60 [+0.01, +1.14]** |
| 800 | 31 | 25 | +0.77 [+0.03, +1.43] | +1.11 [+0.05, +2.49] | +0.60 [+0.02, +1.14] |
| 1,250 | 29 | 24 | +0.85 [+0.06, +1.53] | +1.13 [+0.06, +2.51] | +0.61 [+0.04, +1.14] |
| 1,500 | 24 | 22 | +0.98 [+0.54, +1.39] | +1.10 [−0.08, +2.53] | +0.59 [+0.30, +0.90] |
| 2,000 | 14 | 13 | +1.07 [+0.42, +1.67] | +0.89 [−0.31, +2.36] | +0.38 [+0.12, +0.66] |
| 3,000 | 8 | 8 | +1.61 [+1.01, +2.20] | +0.83 [−0.39, +2.32] | +0.32 [+0.11, +0.58] |

- **The gain is flat from 350 to 1,500.** Against always-PASA it is +1.06 to +1.18 points; against always-BUSCO it is +0.54 to +0.67. Above 1,500 the gain falls, because genomes where PASA training is better are sent to BUSCO.
- **Below 500 complete models:** 5 genomes, mean −8.87. Four lose with PASA training: A. niger (309 models) −16.79, E. xenobiotica (312) −14.20, P. hubeiensis (141) −11.30, S. commune (354) −4.90. A. thermomutatus (436) gains +2.85.
- **No genome has 437-655 complete models.** Thresholds from 450 to 650 give the same result. This data cannot place the threshold within that range.
- **Exon and intron chain at T = 500:** gain vs always-PASA +1.04 [+0.06, +2.40] (exon) and +0.74 [−0.07, +1.84] (intron chain).
- **Conservative T*:** 800 (locus, intron chain), 1,250 (exon). The rule is sensitive to single genomes. At T = 700 the lower bound is −0.05 (locus), because W. jadinii (656 models, +2.19) leaves the ≥ T group. Experiment A, within 4 genomes, gave 300 to more than 1,000.

**Which count tracks the outcome** (Spearman ρ with PASA − BUSCO, locus / exon / intron chain):

| Variable | Locus | Exon | Intron chain |
|---|---|---|---|
| log complete models, training chromosomes (gate variable) | 0.42 | 0.43 | 0.33 |
| log final PASA training models, after selection | 0.53 | 0.45 | 0.38 |
| log complete models, whole genome | 0.56 | 0.54 | 0.46 |

- The whole-genome count is not what is trained in step A. Its higher ρ is *not explained*. One possibility (*inferred, not tested*) is that it reflects the overall depth and quality of the RNA-seq better than the count on half of the genome.

**The same policy with whole-genome counts as the gate variable** [`data/expB_whole_genome_gate_sweep.txt`]. The outcome is still the experiment B score, where training used the training chromosomes only.

| T | Gate variable: training chromosomes. Genomes below / gain vs always-PASA [95% CI] | Gate variable: whole genome. Genomes below / gain vs always-PASA [95% CI] |
|---|---|---|
| 500 | 5 / +1.11 [+0.05, +2.48] | 1 / +0.35 [0.00, +1.06] |
| 700 | 6 / +1.05 [0.00, +2.41] | 2 / +0.48 [0.00, +1.31] |
| 1,000 | 9 / +1.11 [+0.04, +2.49] | 4 / +1.18 [+0.12, +2.52] |
| 2,000 | 26 / +0.89 [−0.31, +2.36] | 6 / +1.19 [+0.12, +2.56] |

- The genomes that lost badly have these whole-genome counts: A. niger 976 (−16.79), P. hubeiensis 843 (−11.30), S. commune 696 (−4.90), E. xenobiotica 482 (−14.20). A. thermomutatus (+2.85) has 6,483. At 500, only E. xenobiotica is below the gate.
- These losses were measured with training on 141-354 models (the training chromosomes). In production these genomes would train on 482-976 models. **That case was not measured.**

**Final training models as an extra gate** [`data/expB_final_models_sweep.txt`, `data/expB_combined_gate_sweep.txt`]:
- Six genomes pass the complete ≥ 500 gate but keep fewer than 300 models after selection. All six are yeasts:

  | Genome | Complete | Final | PASA − BUSCO (locus) |
  |---|---|---|---|
  | M. bicuspidata | 1,456 | 253 | −6.56 |
  | N. castellii | 734 | 175 | −2.44 |
  | K. capsulata | 2,026 | 247 | −1.91 |
  | H. blattae | 1,506 | 169 | −0.49 |
  | W. jadinii | 656 | 281 | +2.19 |
  | M. sympodialis | 717 | 273 | +2.21 |

- Adding "final ≥ 300" to the complete ≥ 500 gate changes mean F1 by +0.18 [−0.17, +0.59] (locus), +0.19 [−0.10, +0.56] (exon) and −0.07 [−0.31, +0.14] (intron chain). All three intervals include 0.
- The existing `--min_training_models` check cannot be used for this. The same value applies after the BUSCO fallback, and predict exits when BUSCO gives fewer models. BUSCO training sets in experiment B hold 90-788 models (10th percentile 218, median 605), so a value of 300 would stop more than 10% of BUSCO-trained genomes. The BFD pipeline uses `--min_training_models 30`.

### 5.3 Conclusions

1. **In experiment B units, 500 works.** *Measured, 40 genomes.* With counts on the training chromosomes, the gate at 500 gains +1.11 locus F1 [+0.05, +2.48] over always training from PASA.
2. **In production units, 500 is not calibrated.** *Measured for the gate variable; the outcome at whole-genome training size is not measured.* Production counts on the whole genome, about 2 times more. With whole-genome counts, 500 catches 1 of the 5 losing genomes, and the gain falls to +0.35 [0.00, +1.06]. At 1,000 the gain is +1.18 [+0.12, +2.52]. But a genome that trains on 976 whole-genome models may do better than it did on 309. Experiment C measures this. Until then, the production default of 500 is uncalibrated, not shown to be wrong.
3. **How the threshold works (experiment B units).** A genome with at least T complete PASA models trains from PASA. A genome with fewer trains from BUSCO. A higher T therefore sends more genomes to BUSCO.
   - **Lower than 500 lets in genomes where PASA training fails.** At T = 300, the genomes with 309-354 models train from PASA and lose badly (A. niger −16.79, E. xenobiotica −14.20, S. commune −4.90).
   - **From 500 to 1,500 the mean barely changes.** Genomes in this range are mixed: some do better with PASA and some with BUSCO, and the differences roughly cancel. The gain stays at +1.10 to +1.13.
   - **Above 1,500 genomes that do better with PASA are sent to BUSCO.** Ten genomes have 1,500-1,999 complete models; 9 of them did better with PASA (mean +0.85). At T = 2,000 they train from BUSCO, and the mean gain falls from +1.10 to +0.89 (10 × 0.85 / 40 = 0.21 points).
   - So in experiment B units, 500 sits at the low end of a flat range. This reasoning applies to training-chromosome counts only (conclusion 2).
   - The conservative T* of 800-1,250 comes from single genomes near the threshold, not from a change in the mean.
4. **A second gate on final training models is not supported yet.** *Measured, 6 genomes.* All intervals include 0. The six genomes are all yeasts, so the effect may depend on selection in intron-poor genomes (*inferred, not tested*).

### 5.4 Limitations

- The thresholds were chosen and evaluated on the same 40 genomes, so the gains are optimistic.
- Only 5 genomes have fewer than 500 complete models, and none has 437-655.
- Each genome has one PASA run (section 4.7).
- The genome set is not a random sample of BFD genomes. It was stratified by the F1 defect and by identity (D111).
- The gate variable and the training set size differ between experiment B (training chromosomes) and production (whole genome). See conclusion 2.
- The correction in sections 5.1-5.3 came from an independent review by a second model (Fable 5.1) on 2026-09-28.

## 6. Open items

- PASA run-to-run noise: repeat the `pasa` arm (about 3 times) on 5-6 genomes near the threshold, including genomes with 437-655 complete models if any can be found.
- Why selection keeps so few training models in some yeasts, and whether a yeast-specific rule helps (section 5.2).
- Low map rate: the map-rate gate still stops train and removes all RNA-seq evidence. Test whether the reads that map help as evidence (C. siamense Cg363, D. hansenii CBS767 with `--min_rnaseq_map_rate 0`).
- Commit the identity-gate code and run one full train → predict test.
- More genomes below 95% read identity with ≥ 500 complete models, to calibrate an identity gate.
- A hints-on vs hints-off arm; GeneMark-ET/EP; refitting EVM weights after the fixes.

## Files in this folder

- `evidence_and_alignment_methods.md`: REVIEW's methods and results draft (evidence, alignment, F1, gate calibration).
- `../training_data_selection_methods.md`: SELECT's methods draft (gates, R3/R5, single-exon training, Tables M1-M8).
- `DECISIONS_snapshot_2026-09-26.md`: a snapshot of the shared decision log (D01-D97).
- `data/`:
  - `scorecard.tsv`: predict arms, holdout gffcompare scores.
  - `single_exon_scores.tsv`: single-exon vs multi-exon exact-match scores.
  - `benchmark.tsv`, `rank_benchmark.tsv`: training models and ranking variants against RefSeq.
  - `titration_scores.tsv`, `titration_analysis_{locus,exon,intron_chain}.tsv`: experiment A, final (149 runs; 4 genomes).
  - `production_f1_scan.tsv.gz`, `production_identity.tsv.gz`: production scans (8,007 genomes).
  - `expB_read_identity.tsv`, `expB_identity_vs_training.tsv`, `expB_identity_analysis.txt`: experiment B read identity and training outcome (section 4).
  - `expB_complete_models.tsv`, `expB_complete_threshold_by_genome.tsv`, `expB_complete_threshold_analysis.txt`, `expB_final_models_sweep.txt`, `expB_combined_gate_sweep.txt`: experiment B complete-model gate (section 5).
- `figures/`: PNG (embedded above) and PDF versions of Figures 1-9. `figures/make_figures.py` regenerates all of them from `data/` (`/usr/bin/python3.12 docs/assessment_pasa2.6_fun1.9/figures/make_figures.py`; needs matplotlib).
- `scripts/`: the analysis scripts that produced the tables (`training_set_vs_refseq.py`, `titration_analysis.py`, `intron_discordance.py`, `production_f1_scan.py`, `production_identity.py`, `predict_scorer.py`, `diversity.py`, the R13 scripts, and `expB_run_identity.sh`, `expB_measure_identity.py`, `expB_analyze_identity.py`, `expB_count_complete.py`, `expB_analyze_complete_threshold.py`, `expB_gate_sweeps.py`).
- Full code review of the PASA fork: hyphaltip/PASApipeline, `CODE_REVIEW_20260925.md` on branch `rust_optimize`.
- Working data (UCR HPCC): `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/` (renamed 2026-09-26 from `do_pasa_rust_vs_perl/`, which is now a symlink).
