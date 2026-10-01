# Assessment: training from PASA vs BUSCO, and the value of RNA-seq evidence (funannotate v1.9.0-rc.3, PASA v2.6.1-rc.2)

This folder records how well funannotate's gene-prediction training and evidence perform after the 2026-09 review of PASApipeline and funannotate. It is written to be read without the original conversation.

- **Versions assessed:** funannotate **v1.9.0-rc.3** (tag `v1.9.0-rc.3`, commit `cd1b5ee`; image `funannotate-1.9.0-rc.3.sif`) and PASApipeline **v2.6.1-rc.2** (hyphaltip/PASApipeline tag `v2.6.1-rc.2`, commit `1044957`). Comparisons against earlier behavior use funannotate `41a2fd7` and the rc.1 image.
- **As of:** 2026-09-26; sections 4 and 5 added 2026-09-28; sections 6 and 7 added 2026-09-29; sections 8-10 added 2026-09-30.
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
5. **The PASA gate of 500 complete models works in experiment B, but it is not yet calibrated for production** (section 5). Experiment B counted models on the training chromosomes (about half the genome). Production counts on the whole genome, which gives about 2 times more (1.5-15 times). On whole-genome counts, 500 would catch only 1 of the 5 genomes that lost badly with PASA training. Experiment C (section 6) tested the whole-genome case: the genomes with 482-976 whole-genome complete models still lost 6.6-18.2 points with PASA training, so 500 is too low in production and about 1,000 separates the 8 genomes tested.
6. **The new identity code passes a wiring test** (section 7): 33 PASS, 0 FAIL; no change to training decisions at default settings.
7. **A full train on the genome's own sequence rescues PASA training when its own Trinity assembly is long** (section 8). E. xenobiotica gains +9.4 and A. niger +6.6 locus F1 over the best production arm. It does not help S. commune or P. hubeiensis, whose own assemblies are also short.
8. **The EVM combination is the main accuracy limit found so far** (sections 9-10). GeneMark-ES alone beats the final models in 7 of 10 genomes, by 2.9-10.6 intron-chain F1 points. For 14-55% of wrong genes, an EVM input already had the exact RefSeq chain. Start/stop errors (4-14% of genes) are not visible in the gffcompare F1 used in sections 1-8.
9. **New EVM weights gain +3.8 to +4.9 F1 on 38 held-out genomes** (section 10). The weights were fitted on two curated references (B. cinerea, C. neoformans H99). The recommended set `augustus:1 hiq:3 genemark:2 pasa:4` is not worse in any of the 76 genome-arm runs. It can be set with `predict -w` today.

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
2. **In production units, 500 is not calibrated.** *Measured for the gate variable. Experiment C (section 6) then measured the outcome at whole-genome training size: 500 is too low.* Production counts on the whole genome, about 2 times more. With whole-genome counts, 500 catches 1 of the 5 losing genomes, and the gain falls to +0.35 [0.00, +1.06]. At 1,000 the gain is +1.18 [+0.12, +2.52]. But a genome that trains on 976 whole-genome models may do better than it did on 309. Experiment C measures this. Until then, the production default of 500 is uncalibrated, not shown to be wrong.
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

## 6. Experiment C: the PASA gate in production units (2026-09-29)

### 6.1 Question and method

- **Question:** section 5 showed that experiment B counted complete PASA models on the training chromosomes, while production counts on the whole genome. Do the genomes that lost with PASA training still lose when Augustus and SNAP train on all of their whole-genome models? If they do, the production gate of 500 is too low.
- **Genomes (8):** the 5 genomes with fewer than 500 complete models on the training chromosomes (A. niger CBS 101883, P. hubeiensis SY62, S. commune H4-8, E. xenobiotica CBS 118157, A. thermomutatus HMR AF 39) and 3 controls (A. nidulans FGSC A4, B. cinerea B05.10, N. crassa OR74A).
- **Arms:** PASA training (`--min_pasa_complete_models 0`) and BUSCO training (`1e9`), **3 repeats each**, on the whole genome. Inputs and flags are those of experiment B (rc.3 image, rc.3 pilot PASA models, BAM and transcripts, BFD predict flags).
- **Scoring:** training uses the whole genome, so no chromosome is held out. Instead, the gene spans of all training models of all 6 runs of a genome form a mask. RefSeq mRNAs and predictions that overlap the mask are removed, and the rest are scored with gffcompare as in experiment B (`scripts/expC_masked_score.py`). Every run of a genome is scored on the same genes, and no training gene is scored.
  - Check: without a mask, on the held-out chromosomes, the scorer reproduces the experiment B score of A. nidulans exactly (5,014 mRNAs, locus 55.8 / 56.7).
  - Bias: the mask removes genes with good evidence, so absolute F1 is lower than in experiment B. The comparison of interest is PASA against BUSCO on the same genes.
- **Jobs:** 16 jobs (3 repeats each; 2 on exfab, 14 on epyc because the exfab user limit is 32 CPUs), all exit 0; 48 runs scored. Files: `data/expC_scores.tsv`, `data/expC_summary.tsv`, `data/expC_analysis.txt`; scripts `scripts/expC_*`.

### 6.2 Results

**Holdout-free (masked) locus F1, mean of 3 repeats (SD), sorted by whole-genome complete models** [`data/expC_analysis.txt`]:

| Genome | Complete models, whole genome | Final PASA training models | RefSeq mRNAs scored / total | PASA | BUSCO | PASA − BUSCO [95% CI] | Experiment B PASA − BUSCO |
|---|---|---|---|---|---|---|---|
| E. xenobiotica | 482 | 343 | 11,360 / 13,187 | 54.78 (0.05) | 69.20 (0.10) | **−14.41** [−14.59, −14.23] | −14.20 |
| S. commune | 696 | 266 | 14,561 / 16,193 | 34.08 (0.00) | 40.70 (0.08) | **−6.62** [−6.75, −6.49] | −4.90 |
| P. hubeiensis | 843 | 321 | 5,936 / 7,472 | 40.47 (0.05) | 55.57 (0.15) | **−15.09** [−15.35, −14.83] | −11.30 |
| A. niger | 976 | 388 | 11,475 / 13,078 | 32.09 (0.07) | 50.26 (0.12) | **−18.17** [−18.39, −17.95] | −16.79 |
| N. crassa | 3,937 | 1,436 | 7,948 / 10,784 | 64.32 (0.06) | 63.91 (0.00) | +0.41 [+0.32, +0.50] | +0.05 |
| A. nidulans | 5,383 | 2,156 | 7,269 / 10,453 | 51.72 (0.06) | 51.09 (0.08) | +0.63 [+0.48, +0.78] | +0.58 |
| A. thermomutatus | 6,483 | 2,449 | 6,131 / 9,701 | 60.70 (0.02) | 59.57 (0.02) | +1.12 [+1.07, +1.18] | +2.85 |
| B. cinerea | 7,041 | 2,676 | 9,006 / 13,703 | 80.11 (0.09) | 77.44 (0.10) | +2.67 [+2.46, +2.88] | +2.35 |

- The 95% CI is from the repeat SDs (t with 4 degrees of freedom). It covers run-to-run noise only, not genome-to-genome variation.
- Exon and intron-chain levels give the same signs in all 8 genomes (`data/expC_analysis.txt`).
- **Repeat noise is small:** SD 0.00-0.09 for PASA training and 0.00-0.15 for BUSCO training. Differences of 0.3 points or more between arms are larger than run-to-run noise.

### 6.3 Conclusions

1. **Training on the full set did not rescue the low-count genomes.** *Measured, 4 genomes.* With 482-976 whole-genome complete models (266-388 final training models), PASA training still loses 6.6-18.2 locus F1 points. In experiment B, with half the models, the loss was 4.9-16.8. So a low complete-model count marks poor PASA data, not only a small training set (*inferred from these 4 genomes*).
2. **The production default of 500 was too low; the default is now 1,000 (user decision, 2026-09-29).** *Measured.* On whole-genome counts it stops only E. xenobiotica. S. commune (696), P. hubeiensis (843) and A. niger (976) pass it and lose 6.6-18.2 points.
3. **A threshold of about 1,000 whole-genome complete models separates the 8 genomes.** The losers have 482-976 and the winners 3,937-7,041. Experiment C has no genome between 977 and 3,936, so it does not place the threshold inside that range.
   - Experiment B has 16 genomes with 977-3,936 whole-genome models. Their PASA − BUSCO (training-chromosome training) ranges from −6.56 to +5.68, and 13 of 16 are within ±2.5 points. The losses among them are yeasts with few final training models on the training chromosomes: M. bicuspidata −6.56 (253 final), N. castellii −2.44 (175), S. cerevisiae −2.04 (392), K. capsulata −1.91 (247).
   - In the experiment B whole-genome sweep (section 5.2), the policy gain is +1.18 at 1,000 and 1,500 and +1.19 at 2,000 (vs +0.35 at 500).
4. **Experiment B's direction holds under whole-genome training.** *Measured.* PASA − BUSCO has the same sign in all 8 genomes, and the size is similar.

### 6.4 Limitations and a new finding about the train path

- 8 genomes; the threshold between 977 and about 2,000 is not pinned down.
- **The complete-model count depends strongly on how train was run.** All 9 rc.3 pilot genomes used here were trained through the BFD pipeline's shared-Trinity branch: `funannotate train --trinity <shared assembly>`, which runs PASA only, with a Trinity assembly made once per species (log line `PASA+PE ... using shared Trinity`, `pasa_tier=relaxed`). Which genome each shared assembly was built on was not checked. The gate wiring test (section 7) re-trained 5 of them with a full train on their own genome, using the same reads:

  | Genome | Complete models, pilot (shared Trinity) | Complete models, full own train |
  |---|---|---|
  | A. nidulans | 5,383 | 5,993 |
  | N. crassa | 3,937 | 3,715 |
  | S. commune | 696 | 522 |
  | P. antarcticum | 3,131 | 102 |
  | E. xenobiotica | 482 | 6,778 |

  - E. xenobiotica's poor PASA set (and its −14.4 loss) may therefore come from the shared-Trinity input, not from its RNA-seq. P. antarcticum's own-reads Trinity produced only 3,083 transcripts, so its full train is poor while its shared-Trinity set is good. *Measured counts. Section 8 measures the accuracy of the own-train PASA sets.*
  - Experiments B and C used the production (shared-Trinity) PASA sets, which is the right input for calibrating a production gate. But the gate result for a given genome can depend on which train path the pipeline chose.

## 7. Gate wiring test (2026-09-29)

- **Purpose:** a software check of funannotate `d18e67c` (identity measurement in train, `--min_rnaseq_identity` in predict; same code as `a9b814f`) against rc.3. It is not an accuracy test.
- **Image:** `funannotate-1.9.0-rc.3+d18e67c.sif`, built from the rc.3 image with only the funannotate Python package replaced (`pip --no-deps`), so all external tools are identical.
- **Runs (16 jobs, all exit 0):** full train with the BFD flags for A. nidulans and N. crassa (old and new image), and for S. commune, P. antarcticum, E. xenobiotica and C. siamense (new image). Predict ran in the production layout: a `training/` folder pruned with the pipeline's own rule, plus a copied `logfiles/` folder. Runs: old vs new predict for A. nidulans and N. crassa; `--min_rnaseq_identity 99` for S. commune and P. antarcticum; the default for E. xenobiotica. Four decision-only bracket runs on N. crassa set thresholds just below and just above the measured values.
- **Result: 41 checks, 33 PASS, 0 FAIL, 8 INFO** [`data/gate_wiring_checks.tsv`; `scripts/wiring_check_wiring.py`]:
  - train writes the read identity, and it matches experiment B exactly (A. nidulans 100.0%, N. crassa 96.67%, S. commune 94.04%, P. antarcticum 92.72%, E. xenobiotica 100.0%). The `training/` copy is written.
  - C. siamense stops at the map-rate gate (exit 3; 0.46% mapped) and still writes its report.
  - Old and new predict on the same training folder: identical training decisions (16 rows, excluding the new identity rows), identical `predict_training_gate.tsv` and identical `final_training_models.gff3`. The identity gate is recorded as disabled by default.
  - `--min_rnaseq_identity 99` switches S. commune (94.04%) and P. antarcticum (92.72%) to BUSCO training, and the PASA evidence stays in the run.
  - Brackets on N. crassa: identity 96.57 → PASA, 96.77 → BUSCO; PASA gate at the count (3,715) → PASA, count + 1 → BUSCO. Both gates trigger in the right direction at the exact boundary.
  - Predict found the identity report in the production layout, through `logfiles/`.
- **INFO items (not failures):**
  - The PASA GFF3 differs between old and new train. `trinity.fasta` already differs (A. nidulans: same 15,012 records; N. crassa: 28,467 vs 28,468). Trinity runs after the gate and the gate does not change its input, so this is consistent with Trinity run-to-run variation. It is not proven: a second rc.3 train was not run.
  - Gene models differ between old and new predict on the same training set. Whole-genome locus F1 differs by 0.1-0.2 points (A. nidulans 55.5/56.2 vs 55.5/56.1; N. crassa 68.4/75.3 vs 68.3/75.1), which is within the repeat noise.
- **Conclusion:** the new code does what it claims and does not change the training decisions or training sets at default settings.

## 8. Own-genome train vs shared-Trinity train (item 3, D121; analysed 2026-09-30)

### 8.1 Question and method

- **Question:** section 6.4 found that the complete-model count depends on the train path. Does a full train on the genome's own sequence (own genome-guided Trinity, then PASA) rescue PASA training for the genomes that lost?
- **Arms (3 repeats each):** the experiment C arms `pasa` and `busco` use the production shared-Trinity PASA set. The new arms `pasa_own` and `busco_own` use the PASA models, BAM and transcripts from a full own train with the same reads. `busco_own` trains from BUSCO, so it changes only the evidence.
- **Genomes (8):** the 4 experiment C losers, T. versicolor and N. castellii (short shared assemblies, above the 1,000 gate), and A. nidulans and N. crassa as controls. A. thermomutatus and B. cinerea have shared arms only.
- **Scoring:** as in experiment C. The mask is the union of the training models of all 4 arms, so the gene set differs slightly from section 6 (`data/expC_item3_scores.tsv`).
- **Repeat noise:** SD ≤ 0.19 locus F1 in every arm.

### 8.2 Results

**Masked locus F1, mean of 3 repeats** [`data/expC_item3_analysis.txt`]. Median lengths are Trinity transcript medians (shared from `wave1_triage.tsv`; own measured from `funannotate_train.trinity-GG.fasta`):

| Genome | Median bp, shared / own | Complete models, shared / own | pasa | busco | pasa_own | busco_own | Best own − best shared |
|---|---|---|---|---|---|---|---|
| E. xenobiotica | 435 / 1,854 | 482 / 6,778 | 52.71 | 67.57 | **76.98** | 74.58 | **+9.41** |
| S. commune | 642 / 656 | 696 / 522 | 33.52 | **40.20** | 34.46 | 40.15 | −0.06 |
| P. hubeiensis | 521 / 438 | 843 / 144 | 40.28 | **55.49** | 29.61 | 54.71 | −0.78 |
| A. niger | 549 / 2,077 | 976 / 7,147 | 29.90 | 47.73 | **54.31** | 51.72 | **+6.58** |
| N. castellii | 721 / 729 | 1,998 / 2,557 | 77.94 | 79.46 | 82.02 | **83.31** | **+3.85** |
| T. versicolor | 629 / 649 | 2,745 / 3,378 | 44.00 | 43.88 | **44.21** | 43.56 | +0.21 |
| N. crassa | 912 / 912 | 3,937 / 3,715 | 64.06 | 63.72 | **65.03** | 64.43 | +0.97 |
| A. nidulans | 1,535 / 1,724 | 5,383 / 5,993 | **51.20** | 50.50 | 50.67 | 50.31 | −0.53 |

- Exon and intron-chain levels give the same signs for the "best own − best shared" column in 7 of 8 genomes. The exception is S. commune at the exon level (0.00).
- N. crassa's shared and own assemblies are almost the same (28,467 and 28,468 transcripts, same median). So its difference comes from the PASA runs, not the assembly.

### 8.3 Conclusions

1. **Where the own assembly is long, the own train rescues PASA training.** *Measured, 2 genomes.* In E. xenobiotica and A. niger the own assembly median is 1,854 and 2,077 bp, against 435 and 549 bp shared. PASA training then beats BUSCO training by +2.40 and +2.59. The best arm is 9.41 and 6.58 points above the best production arm.
2. **Where the own assembly is also short, the own train does not help.** *Measured, 2 genomes.* In S. commune (656 bp) and P. hubeiensis (438 bp), BUSCO training stays best. P. hubeiensis own PASA training is worse (−25.10 against busco_own; 144 complete models). So for these two genomes the short assembly comes from the reads or the genome, not from the shared path. The cause is not known.
3. **Better PASA evidence helps even with BUSCO training.** *Measured.* busco_own − busco is +7.01 (E. xenobiotica), +3.99 (A. niger) and +3.85 (N. castellii). Only the evidence differs between these arms.
4. **The own-assembly median length separates the outcomes better than the shared one.** In the 6 genomes with a short shared assembly, the 2 with own median ≥ 800 bp gain. N. castellii (729 bp) gains from evidence only (BUSCO training still best). T. versicolor (649 bp) is flat. *Inferred from 8 genomes; no threshold is fitted.*
5. **Implication for BFD (not tested in the pipeline):** a genome whose shared assembly is short should get its own genome-guided Trinity assembly before falling back to BUSCO. If the own assembly is also short, the D120 rule (train from BUSCO, keep evidence) applies.

## 9. Error breakdown of the final gene models (2026-09-30)

### 9.1 Method

- **Runs:** replicate 1 of every arm in section 8 and experiment C: 36 predict runs, 10 genomes. The repeat SD is ≤ 0.19 F1, so one replicate per arm is enough.
- **Gene set:** the same mask and sequences as the masked scorer.
- **Classes** (`scripts/expC_error_breakdown.py`). Each scored RefSeq gene gets the first class that fits. Overlap means same-strand CDS overlap.
  - exact: a final model has the exact CDS chain of a RefSeq isoform.
  - missed: no final model overlaps the gene.
  - merged: one final model also overlaps another RefSeq gene.
  - split: two or more final models overlap the gene.
  - For 1:1 pairs:
    - ends_wrong: the same intron chain with a different start or stop;
    - exons_fewer and exons_more: fewer or more CDS segments;
    - splice_diff: the same number of segments with a boundary shifted.
- **Input check:** for each non-exact gene, did any EVM input (`gene_predictions.gff3`: Augustus, HiQ, GeneMark, SNAP, PASA) have the exact chain?
- **Per-source accuracy** (`scripts/expC_source_accuracy.py`): exact-chain sensitivity and precision of each EVM input alone, of `evm.round1.gff3`, and of the final models.
- **Checks:**
  - Every input source ends its CDS with a stop codon in 88-100% of models, the same convention as RefSeq. So the exact-chain comparison is not biased by the stop codon.
  - For E. xenobiotica busco_r1, exact + ends_wrong = 66.9% of genes, against gffcompare locus Sn 65.7.

### 9.2 Results

**Error classes, % of scored RefSeq genes, production (shared-path) BUSCO arm** [`data/expC_error_breakdown_summary.txt`; all arms in the same file]:

| Genome | exact | ends wrong | splice diff | exons fewer | exons more | split | merged | missed |
|---|---|---|---|---|---|---|---|---|
| A. niger | 38.1 | 4.8 | 6.3 | 12.7 | 7.2 | 2.6 | 0.8 | 27.5 |
| P. hubeiensis | 38.7 | 14.1 | 2.3 | 7.9 | 9.0 | 1.6 | 3.4 | 23.0 |
| S. commune | 30.6 | 4.8 | 6.8 | 11.8 | 9.8 | 2.6 | 7.8 | 26.0 |
| E. xenobiotica | 57.1 | 9.8 | 5.3 | 9.2 | 7.0 | 1.3 | 2.2 | 8.2 |
| A. thermomutatus | 58.6 | 6.8 | 6.2 | 8.5 | 10.8 | 6.3 | 0.2 | 2.7 |
| A. nidulans | 45.5 | 3.9 | 10.1 | 16.3 | 7.7 | 2.3 | 1.4 | 12.8 |
| B. cinerea | 72.1 | 6.7 | 2.9 | 5.2 | 6.1 | 0.5 | 0.3 | 6.2 |
| N. crassa | 53.9 | 6.8 | 4.3 | 8.5 | 8.4 | 1.8 | 0.4 | 15.9 |
| T. versicolor | 35.3 | 6.9 | 8.5 | 10.6 | 15.5 | 2.2 | 5.6 | 15.5 |
| N. castellii | 68.9 | 11.7 | 0.3 | 0.3 | 7.7 | 0.1 | 3.0 | 8.0 |

**GeneMark-ES alone against the final models, intron-chain F1** (ends ignored for multi-exon genes, as gffcompare does) [`data/expC_genemark_vs_final_intron_chain.txt`]:

| Genome | GeneMark-ES alone | Final, best shared-path arm | Final, best own-path arm |
|---|---|---|---|
| A. niger | **55.1** | 46.2 | 52.7 |
| P. hubeiensis | **50.9** | 44.5 | 44.1 |
| S. commune | **44.3** | 39.8 | 39.7 |
| E. xenobiotica | 74.1 | 63.5 | **74.6** |
| A. thermomutatus | **64.4** | 58.6 | not run |
| A. nidulans | 48.8 | **49.9** | 49.1 |
| B. cinerea | 74.7 | **78.4** | not run |
| N. crassa | **68.8** | 62.3 | 63.3 |
| T. versicolor | **48.2** | 43.5 | 43.7 |
| N. castellii | 74.4 | 71.2 | **77.6** |

**Where the missed genes are lost** (BUSCO arm; `data/expC_source_accuracy.tsv`, `data/expC_missed_repeat_check.txt`, `data/expC_noinput_contig_check.txt`):
- A large share of missed genes has no overlapping model from any predictor: A. niger 1,636 (19% of scored genes), P. hubeiensis 1,103 (19%), S. commune 2,312 (16%), N. crassa 655 (9%), A. nidulans 528 (8%).
  - Repeat masking does not explain them: ≤ 7% of these genes have ≥ 50% of their CDS in repeatmasker spans.
  - Contig length does not explain them either: most are on contigs ≥ 500 kb.
  - Most have no transcript or protein alignment.
- The rest of the missed genes are overlapped only by GeneMark or SNAP models (weight 1), and EVM does not call them. Examples: A. niger 650, S. commune 911, E. xenobiotica 335.
- The filter after EVM removes few scored genes, except in the two basidiomycetes: S. commune 386 and T. versicolor 385 (about 2.7% of scored genes). Almost all were removed with `remove_reason=repeat_match` (diamond hit to the repeat protein database).

### 9.3 Conclusions

1. **EVM output is worse than GeneMark-ES alone in 7 of 10 genomes** on the production arm, by 2.9 to 10.6 intron-chain F1 points. *Measured.* The final models are better only in A. nidulans (+1.1) and B. cinerea (+3.7). With a good own PASA set, E. xenobiotica (+0.5) and N. castellii (+3.2) also beat GeneMark.
2. **In 14-55% of non-exact genes, an EVM input already had the exact chain.** Most often that input was GeneMark (35,583 of 47,699 source hits over 36 runs). The combination step loses correct models that are present in its inputs. *Measured.* This is an upper bound: no real selector can pick the right input every time.
3. **Augustus from this training is weak on exact chains.** Augustus + HiQ exact sensitivity is 17-66% of GeneMark's in every genome. HiQ models are precise (48-83%) but few. *Measured.* Where the training data are good (B. cinerea, E. xenobiotica own), the final models beat GeneMark.
4. **Wrong start or stop affects 4-14% of genes.** gffcompare's locus, transcript and intron-chain levels do not count this error. So the F1 values in sections 1-8 do not measure it. It is largest in P. hubeiensis (14.1%) and N. castellii (11.7%).
5. **Genes with no prediction from any tool are the largest single loss** in A. niger, P. hubeiensis and S. commune. The cause is not known. One unverified possibility is that some of these RefSeq models are not real genes.
6. **Caveat on the GeneMark comparison.** Several of these RefSeq annotations come from JGI or Broad pipelines. The JGI pipeline uses GeneMark among its predictors. So agreement with GeneMark-ES may be inflated for those genomes. Which RefSeq sets were built with GeneMark was not checked.

## 10. EVM weight refit (2026-09-30)

### 10.1 Method

- **Inputs:** the saved EVM inputs of experiment B step-B runs. These predictors were trained on the training chromosomes and predicted the whole genome. Arms: `pasa.B` (PASA training) and `busco_r1.B` (BUSCO training), with the same fixed evidence.
- **Rerun:** only `funannotate-runEVM.py` is rerun (rc.3 image, Rust EVM, `-m 10 -i 1500`), with a new weights file (`scripts/evmrefit_evm_rerun.sh`). The post-EVM filters are not rerun. `evm.round1` and the final models differ by ≤ 1 F1 point (H99: 79.28 vs 80.01).
- **Scoring:** the holdout chromosomes against RefSeq, at the gene level (`scripts/evmrefit_evm_score.py`):
  - **ic:** intron-chain F1. It uses the gffcompare convention, so the ends are free for multi-exon genes.
  - **ex:** exact-CDS F1. The start and stop must also match.
- **Reproducibility:** the Rust EVM engine is not deterministic. Three runs with identical inputs differ by ±0.04 F1 (H99). So every comparison is against a rerun with the current weights ("base"), not against the stored output.
- **Fit set:** B. cinerea B05.10 and C. neoformans H99. Their references are curated (per the user), so they do not depend on one predictor.
  - Round 1: 29 sets (one-factor changes and a few combinations).
  - Round 2: 54 sets around the best, varying GeneMark 2-4, Augustus 1-2, HiQ 3-5 and PASA 3-6.
- **Validation set:** the other 38 experiment B genomes (both arms). They were not used to choose the weights. The finalists were 5 sets plus base.

### 10.2 Results

**Each source alone on the fit genomes (holdout, ic F1)** [`data/evmrefit_pilot_source_scores.tsv`]:

| | GeneMark-ES | Augustus + HiQ (PASA-trained) | PASA | SNAP | EVM, current weights |
|---|---|---|---|---|---|
| B. cinerea | 77.5 | 76.4 | 64.0 | 41.3 | 83.0 |
| C. neoformans H99 | 69.8 | 71.0 | 61.7 | 27.5 | 79.3 |

On these two curated references, EVM beats GeneMark alone by 5.5 and 9.5 points, unlike in section 9.

**Round 1, one-factor effects (mean over the 4 fit runs, ic F1 vs base)** [`data/evmrefit_pilot_analysis.txt`]:
- GeneMark 0 → −6.19; GeneMark 2 → +1.43; GeneMark 3 → +1.45; GeneMark 8 → +0.26.
- Augustus 2 → −2.36, but Augustus 2 + HiQ 4 + GeneMark 3 → +2.24. The ratio between the sources matters, not the absolute values.
- PASA 3 → +1.15; PASA 10 → −0.95; PASA 15 → −1.46.
- SNAP 0 → −1.63; SNAP 2 → −4.27.
- Protein and transcript weights: −0.34 to +0.06.

**Held-out validation, 38 genomes, delta vs current weights (mean [95% bootstrap CI over genomes]; genomes better / worse by > 0.1)** [`data/evmrefit_validation_summary.txt`]:

| Weights (Augustus/HiQ/GeneMark/SNAP/PASA/proteins/transcripts) | PASA arm, ic | PASA arm, ex | BUSCO arm, ic | BUSCO arm, ex |
|---|---|---|---|---|
| 2/5/4/1/6/1/1 | +4.82 [+3.53, +6.43]; 37 / 0 | +4.94; 38 / 0 | +4.13 [+3.38, +4.88]; 36 / 2 | +4.22; 36 / 2 |
| 1/3/2/1/3/1/1 | +4.74 [+3.54, +6.23]; 38 / 0 | +4.76; 38 / 0 | +4.06 [+3.33, +4.80]; 36 / 1 | +4.08; 36 / 0 |
| 2/3/3/1/4/1/1 | +4.73 [+3.52, +6.22]; 38 / 0 | +4.74; 38 / 0 | +4.02 [+3.31, +4.72]; 36 / 1 | +4.05; 36 / 1 |
| **1/3/2/1/4/1/1** | +4.44 [+3.21, +5.93]; 37 / 0 | +4.41; 38 / 0 | +3.78 [+3.12, +4.44]; **38 / 0** | +3.76; **38 / 0** |
| 2/3/2/1/4/1/1 | +2.31 [+1.99, +2.65]; 38 / 0 | +2.28; 38 / 0 | +2.32 [+1.94, +2.69]; 36 / 1 | +2.24; 37 / 1 |

- **The largest gains come where the PASA set is poor** (PASA arm): A. niger +20.7 to +22.6, E. xenobiotica +19.3 to +20.7, P. hubeiensis +11.3 to +11.8, S. commune +9.4 to +10.1. With PASA at weight 6, poor PASA models outweigh the ab initio models.
- **The smallest gains come on the references that are probably most curated:** S. cerevisiae S288C −0.31 to +2.24 and A. nidulans FGSC A4 +0.18 to +0.73. On the fit genomes the gain was +0.9 to +3.1 (in-sample).
- **Gene counts change little.** Where PASA is poor, predictions move towards the RefSeq count (A. niger PASA arm 1,320 → 2,013-2,018 against 2,216 RefSeq).
- **The only losses** are in the BUSCO arm: M. bicuspidata −0.75 and S. cerevisiae −0.31 with 2/5/4/1/6. The set 1/3/2/1/4 has no loss in any genome or arm.

### 10.3 Conclusions

1. **The current weights undervalue GeneMark and HiQ against PASA and plain Augustus.** *Measured, 40 genomes.* Four weight sets give +3.8 to +4.9 mean F1 on 38 held-out genomes, in both arms and at both levels. Their CIs overlap, so the data do not separate them.
2. **Recommended set: `augustus:1 hiq:3 genemark:2 snap:1 pasa:4`** (proteins 1, transcripts 1). It is the only finalist that is not worse by more than 0.1 in any held-out genome or arm (76 of 76). Its mean gain is +3.8 to +4.4, which is 0.3-0.4 below the best set. *This is a choice of safety over mean gain.*
3. **The truth-set caveat still applies to the size of the gain.** On the references that are probably most curated (S. cerevisiae, A. nidulans), the gain is small (about 0 to +2). Most RefSeq sets in the validation were not checked for provenance. The direction is the same on every reference.
4. **No code change is needed to use the new weights.** `funannotate predict -w augustus:1 hiq:3 genemark:2 pasa:4` sets them. The BFD pipeline passes `-w codingquarry:0 glimmerhmm:0 genemark:1` today. Changing funannotate's built-in defaults (`predict.py` StartWeights) is a separate decision.
5. **Not tested:** GlimmerHMM and CodingQuarry weights (both 0 in BFD); `--repeats2evm`; experiment C whole-genome training; own-train PASA sets.

### 10.4 Genomes without RNA-seq, and full-predict confirmation (D127)

- **Without RNA-seq** (experiment B `norna` arm, 15 genomes, EVM-only, intron-chain F1 vs current weights) [evm_refit/runs/*/norna.B]: `augustus:1 hiq:3 genemark:2 snap:1` +8.59 (15 of 15 better); `augustus:2 hiq:5 genemark:4 snap:1` +9.52 (15 of 15). The choice between the two used the same 15 genomes.
- **Full predict, post-EVM filters included** (experiment B design, PASA arm, 2 repeats per weight set; `data/evmrefit_confirm_summary.txt`). Current BFD weights vs the refit set:

  | Genome (RefSeq) | Intron-chain F1 | Exact-CDS F1 | gffcompare locus F1 | Genes | Repeat range |
  |---|---|---|---|---|---|
  | B. cinerea B05.10 | 83.17 → 84.34 (+1.16) | +1.20 | +1.31 | +125 | ≤ 0.13 |
  | C. neoformans H99 | 80.07 → 83.16 (+3.09) | +3.00 | +3.10 | −1 | ≤ 0.04 |
  | P. ostreatus PC9 | 47.00 → 51.60 (+4.60) | +4.55 | +4.76 | −180 | 0.10 |
  | P. blakesleeanus NRRL 1555 | 38.62 → 42.50 (+3.88) | +3.83 | +3.99 | +72 | ≤ 0.09 |

  - PC9 and NRRL 1555 were trained here with the production shared-Trinity path; their RefSeq GFFs were downloaded from NCBI. NRRL 1555's RefSeq set is a JGI annotation. Both have low absolute F1; the reason was not checked.
- **Adopted (user, 2026-09-30):** BFD `predict_evm_weights` = `augustus:1 hiq:3 genemark:2 snap:1 pasa:4` (genomes with a PASA set) and `predict_evm_weights_norna` = `augustus:2 hiq:5 genemark:4 snap:1`, in Fungi_BFD and nf_funannotate1. funannotate's built-in StartWeights are unchanged.

## 11. Open items

- **Decided (user, 2026-09-29): the `--min_pasa_complete_models` default is now 1,000** ("sure set --min_pasa_complete_models to 1000 if that is justified, may be a high bar but is show the accuracy need"). 1,000 may be conservative, because no tested genome had 977-3,936 complete models.
- Genomes with 1,000-3,900 whole-genome complete models under whole-genome training, to place the threshold more exactly.
- **Done (section 8):** shared-Trinity vs own-train PASA sets. Open: why S. commune and P. hubeiensis give short assemblies on both paths; whether the BFD pipeline should build an own genome-guided assembly when the shared one is short.
- **EVM weights (section 10):** a pending user decision: adopt `-w augustus:1 hiq:3 genemark:2 pasa:4` in BFD and/or as funannotate defaults. Before that, confirm with full predict runs (post-EVM filters included) on a few genomes. Not yet tested: GlimmerHMM/CodingQuarry weights and `--repeats2evm`.
- **Truth-set check:** find which RefSeq annotations were built with GeneMark (JGI/Broad sources), and repeat the section 9 comparison without them.
- **Genes with no prediction from any tool** (16-19% of scored genes in A. niger, P. hubeiensis, S. commune): check whether they have homology or expression support, to separate real misses from doubtful RefSeq models.
- **Start/stop accuracy:** add an exact-CDS level to the scorer. gffcompare's intron-chain level does not count wrong starts or stops (4-14% of genes).
- **repeat_match filter in basidiomycetes:** 2.7% of scored RefSeq genes in S. commune and T. versicolor are removed after EVM. Check whether they are TE-derived.
- Why selection keeps so few training models in some yeasts, and whether a yeast-specific rule helps (section 5.2).
- Low map rate: the map-rate gate still stops train and removes all RNA-seq evidence. Test whether the reads that map help as evidence (C. siamense Cg363, D. hansenii CBS767 with `--min_rnaseq_map_rate 0`).
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
  - `expC_scores.tsv`, `expC_summary.tsv`, `expC_analysis.txt`: experiment C (section 6). `gate_wiring_checks.tsv`: the gate wiring test (section 7).
  - `expB_complete_models.tsv`, `expB_complete_threshold_by_genome.tsv`, `expB_complete_threshold_analysis.txt`, `expB_final_models_sweep.txt`, `expB_combined_gate_sweep.txt`: experiment B complete-model gate (section 5).
- Sections 8-9 (2026-09-30): `expC_item3_scores.tsv`, `expC_item3_summary.tsv`, `expC_item3_analysis.txt` (item 3); `expC_error_breakdown.tsv`, `expC_error_breakdown_summary.txt`, `expC_source_accuracy.tsv`, `expC_genemark_vs_final_intron_chain.txt`, `expC_missed_repeat_check.txt`, `expC_noinput_contig_check.txt` (error breakdown). Per-gene classes: `experiment_C/error_breakdown_genes/` in the working data.
- Section 10: `evmrefit_*` data (weight sets, pilot analysis, source scores, validation deltas and summary) and `scripts/evmrefit_*` (they run from the working-data `evm_refit/` folder under their names without the prefix).
- `figures/`: PNG (embedded above) and PDF versions of Figures 1-9. `figures/make_figures.py` regenerates all of them from `data/` (`/usr/bin/python3.12 docs/assessment_pasa2.6_fun1.9/figures/make_figures.py`; needs matplotlib).
- `scripts/`: the analysis scripts that produced the tables (`training_set_vs_refseq.py`, `titration_analysis.py`, `intron_discordance.py`, `production_f1_scan.py`, `production_identity.py`, `predict_scorer.py`, `diversity.py`, the R13 scripts, and `expB_run_identity.sh`, `expB_measure_identity.py`, `expB_analyze_identity.py`, `expB_count_complete.py`, `expB_analyze_complete_threshold.py`, `expB_gate_sweeps.py`, `expC_*` and `wiring_*`). The `expC_*` scripts run from the working-data `experiment_C/` folder under their names without the `expC_` prefix, because they import each other (`masked_score`, `error_breakdown`, `source_accuracy`).
- Full code review of the PASA fork: hyphaltip/PASApipeline, `CODE_REVIEW_20260925.md` on branch `rust_optimize`.
- Working data (UCR HPCC): `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/` (renamed 2026-09-26 from `do_pasa_rust_vs_perl/`, which is now a symlink).
