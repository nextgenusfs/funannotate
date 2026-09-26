# Selection and validation of ab initio training data in funannotate (draft methods text)

Status: draft for the funannotate paper, 2026-09-26. Data location: the result files named below are in the BFD project directory `/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/do_pasa_rust_vs_perl/` (UCR HPCC); they are not part of this repository. Every number comes from the result files named in each table caption. The tables are produced by `methods_tables.py`, not typed by hand. Decision history: `DECISIONS.md` (entry numbers D07-D82). Code: funannotate-live working tree based on commit `41a2fd7` (uncommitted at the time of writing).

## 1. What problem this solves

funannotate trains Augustus and SNAP from gene models that PASA builds from RNA-seq transcripts. When those models are wrong, the ab initio predictors learn wrong gene structure, and the final gene set loses genes. In an 18-genome pilot, retraining on PASA models reduced BUSCO completeness of the predicted proteome by up to 60 points (for example *Drepanopeziza brunnea* 97.1% → 46.8%, *Ophiostoma quercus* 98.8% → 37.1%). We traced this to three causes:

1. Many PASA training models were partial open reading frames (ORFs). Only 5-26% were complete ORFs in filamentous fungi, against about 75% in yeasts.
2. Some RNA-seq did not come from the annotated genome (infected host tissue, a different species, or an out-of-date read file).
3. A PASApipeline defect (duplicate GFF3 rows, "F1") corrupted the genome coordinates of PASA training models between 2026-07-06 and 2026-09-24.

We then changed how funannotate decides which data trains the ab initio predictors, and measured each change against RefSeq annotations.

## 2. How funannotate now chooses its training data

The decisions happen in a fixed order. Each decision is written to `logfiles/training_decisions.tsv` (columns: command, stage, decision, value, threshold, outcome, reason, timestamp). It is also echoed to the run log as a line that starts with `TRAINING-DECISION`, and predict prints a summary table before Augustus training.

**funannotate train**
1. **RNA-seq concordance gate** (`--min_rnaseq_map_rate`, default 10%). The first 200,000 reads (`--rnaseq_gate_reads`) are mapped to the genome with minimap2 (`-x splice:sr`). If fewer than 10% map with MAPQ ≥ 1, train stops before Trinity and PASA with exit code 3. A workflow can then route the genome to prediction without RNA-seq. Stage name: `rnaseq_gate`. Report: `logfiles/train_rnaseq_gate.tsv`.
2. **PASA options.** The opt-in PASA flags `--pasa_unspliced_join_spliced` and `--pasa_one_alignment_per_cdna` are passed only if the installed PASA launcher supports them (PASApipeline ≥ v2.6.1-rc.2). Stage name: `pasa_options`.
3. **One PASA model per locus** (`getBestModel`). PASA and TransDecoder produce several models per locus. Models on the same strand whose transcript spans overlap by at least `--pasa_alignment_overlap` (30%) of the shorter model are joined into transitive clusters. One model is kept per cluster, ranked by: complete ORF, number of CDS exons, CDS length, transcript abundance (TPM), then gene ID.

**funannotate predict**
4. **Initial training source per predictor** (`training_mode_initial`): pre-trained parameters, PASA models, or BUSCO models.
5. **PASA training-set gate** (`--min_pasa_complete_models`, default 500). The number of complete-ORF models in the PASA GFF3 is counted: ATG start, one stop codon at the end, and a CDS length divisible by 3. If it is below the threshold, the predictors that would train from PASA train from BUSCO instead. The PASA models are still used as EVM evidence. The gate is skipped when no predictor trains from PASA, or when Augustus output already exists. Stage name: `pasa_gate`. Report: `logfiles/predict_training_gate.tsv`.
6. **Training-set selection** (`selectTrainingModels`), one record per step:
   - `select_complete_orf`: partial ORFs are removed;
   - `select_keeper_filter`: models whose introns match RNA-seq or protein hints ("keepers", from filterGenemark.pl) are preferred when at least 200 exist;
   - `select_multi_cds`: multi-exon models are required when at least 200 exist;
   - `select_single_exon_support` / `select_single_exon_admit`: see step 7;
   - `select_redundancy`: redundant models are removed (DIAMOND, ≥80% identity and coverage);
   - `select_overlap`: one model is kept per transitive overlap cluster on either strand;
   - `select_final`: the number of training models written.
7. **Single-exon training genes** (default on; `--no_training_single_exon` turns it off).
   - Without this step, the multi-exon requirement removes every single-exon gene, although about 14-22% of fungal RefSeq genes are single-exon.
   - Complete single-exon PASA models whose CDS is at least 80% covered by a same-strand protein2genome alignment are admitted.
   - They are capped at share/(1 − share) × the number of multi-exon training models. The share is the fraction of single-CDS GeneMark-ES gene models (fallback: near-full-length protein alignments).
   - Stage names: `single_exon_share`, `select_single_exon_admit`.
8. **Minimum training models** (`--min_training_models`). If too few models remain, training falls back to BUSCO. Stage name: `min_training_models`.
9. **Final training source** per predictor and the number of models (`final_training_source`).

## 3. How we evaluated the changes

- **Reference genomes.** Five fungal genomes with RefSeq annotation: *Neurospora crassa* OR74A, *Aspergillus nidulans* FGSC A4, *Botrytis cinerea* B05.10, *Cryptococcus neoformans* H99 and *Schizophyllum commune* H4-8.
- **RNA-seq identity.** Measured as the median NM-based identity of 50,000 mapped reads:
  - same strain (≥ 99.5%): *A. nidulans*, *B. cinerea*, *C. neoformans*;
  - divergent: *N. crassa* (reads from strain HJDF, 96.7%) and *S. commune* (94.0%).
- **Held-out chromosomes.** Chromosomes ≥ 1 Mb were sorted by length and assigned alternately to a training set and a holdout set; the mitochondrion was excluded.
  - Step A trains Augustus and SNAP on the training chromosomes only.
  - Step B predicts the whole genome with the trained parameters (`-p`).
  - Only holdout chromosomes are scored, so training genes are never scored.
- **Two evidence settings in step B:**
  - "fixed" gives every variant the same EVM evidence, which isolates the training effect;
  - "own" gives each variant its own PASA models as EVM evidence, which is closest to production.
- **Scores:**
  - gffcompare (v0.12.10) at the CDS level: locus, exon and intron-chain sensitivity (Sn) and precision (Pr);
  - an exact CDS-chain match to RefSeq protein-coding mRNAs, split into single-exon and multi-exon genes;
  - BUSCO (v5.8.0, fungi_odb10, protein mode) of the predicted proteome.
- **Fixed software.** All variants ran in the same container (funannotate-1.9.0-rc.1, PASApipeline with the F1 fix). Only the funannotate Python code differed, bind-mounted from frozen snapshots, so the ab initio and EVM binaries were identical.
- **PASA inputs.** "rc1" is the standard rc.1 PASA output. "R1" adds a corrected minimap2 alignment parser. "R2" lets unspliced alignments join spliced clusters.

## 4. Results

- **Locus rule (Table M1).** Clustering on transcript spans (tx_strand) was at least as good as clustering on CDS spans on all five genomes. The CDS-span rules added 33-100 predicted genes, changed locus Sn by −0.7 to +0.3 points, and lowered locus Pr by 0.4-2.6 points.
- **Ranking (Tables M3, M8).**
  - Guarded ranking gave 117-235 more exact RefSeq intron chains per genome in the PASA models passed to EVM (Table M8).
  - It did not change the final predictions by more than 1 point in either direction (Table M3), and it gave fewer training models.
  - We kept complete-first as the default. Guarded ranking is available as a parameter.
- **Training-set construction (Tables M2, M7).**
  - Against RefSeq, the new selection raised the fraction of training models that exactly match a RefSeq gene: *N. crassa* 38.9% → 78.7%, *A. nidulans* 50.1% → 62.1%, *B. cinerea* 75.5% → 89.4% (Table M7). This benchmark had no keeper set, so it is an upper bound.
  - In real runs with hints, the keeper filter already removed most partial models when ≥ 200 keepers existed. The new selection then changed final accuracy by ≤ 1 point (Table M2).
  - On *S. commune* (30 keepers, 93 complete models), removing partial models without the gate reduced locus Sn from 29.2 to 24.9 (Table M2).
- **PASA gate (Table M6).** On *S. commune*, BUSCO training (what the gate chooses at 93 < 500 complete models) gave locus 37.4 / 49.2 and proteome BUSCO 97.5%, against 29.2 / 44.2 and 86.9% for the previous production path. Removing partial models and the gate must therefore be used together.
- **Single-exon training genes (Table M4).** Admitting protein-supported single-exon genes raised single-exon Sn by 5.2-7.6 points on four genomes. Multi-exon Pr rose by 0.4-2.4 points. Multi-exon Sn changed by 0.0 to −0.9. Single-exon Pr fell by 0.8-5.3. On *S. commune* no single-exon gene had protein support, so nothing changed.
- **Divergent RNA-seq (Table M5).** On *N. crassa*:
  - PASA training without single-exon genes stayed below BUSCO training, even with the corrected parser (fixed evidence 64.6 / 73.0 vs 66.0 / 74.1).
  - Adding single-exon training genes to the corrected PASA models gave the best result in both evidence settings (fixed 66.4 / 74.9; own 67.1 / 75.0).
  - Relaxing PASA validation (identity 90%, no splice-site requirement) added training models but did not improve predictions.
- **Same-strain RNA-seq.** BUSCO training was worse than PASA training on *C. neoformans* (78.1 / 80.6 vs 79.9 / 80.8) and *B. cinerea* (82.0 / 82.0 vs 83.3 / 82.2), and about equal on *A. nidulans* (Table M2).
- **Choice of the single-exon share estimator.** Against RefSeq single-CDS genes on the training chromosomes (*N. crassa*, *A. nidulans*, *B. cinerea*: 22.0%, 14.0%, 21.5%):
  - GeneMark-ES gave 25.5%, 20.5% and 22.8% (mean absolute error 3.8 points);
  - near-full-length protein alignments gave 11.8%, 11.5% and 7.9% (error 8.8 points).
  - So GeneMark-ES is the default source.

## 5. Scale in the BFD production set

- **Genomes checked:** 21,762 production genomes with a predict log.
  - 8,010 trained Augustus on PASA models;
  - 1,595 of these had fewer than 500 complete PASA models.
- **Read concordance.** 55 of 1,439 species mapped fewer than 10% of reads to their own genome.
- **F1 defect.** 5,599 PASA training sets dated 2026-06-29 to 2026-09-24 carried CDS lengths not divisible by 3. Sets dated before 2026-06-29 had none (1,116 sets, maximum 0.000). Another 3,052 genomes reused ab initio parameters from an affected representative genome.
- **Low-keeper genomes.** Fewer than 200 keepers occurred in 30% of F1-affected genomes and in 10% of clean genomes, in a random sample of 400.

## 6. Limitations

- Each comparison uses one genome and one run, with no replicates. Differences of 1-2 points show direction, not significance.
- The 500-model gate threshold was chosen before the F1 fix. It is being calibrated in a titration experiment:
  - 4 genomes; N = 50-2,000 complete training models; 3-5 random draws per N;
  - rule: the smallest N where the lower 95% bound of (PASA-trained − BUSCO-trained) holdout locus F1 is ≥ 0.
  - A cross-genome validation on about 40 RefSeq genomes follows.
- The Swiss-Prot protein evidence contains curated proteins of these model species, so protein support for single-exon genes is strongest for them.
- RefSeq annotations are the reference, but they are not error-free.

## 7. Tables

Generated by `methods_tables.py` from `predict_arms/scorecard.tsv`, `predict_arms/single_exon_scores.tsv`, `refseq_benchmark/benchmark.tsv` and `refseq_benchmark/rank_benchmark.tsv` on 2026-09-26.

### Table M1. Locus rule for one PASA model per locus (step B with each rule's own PASA evidence)

Holdout-chromosome gffcompare, CDS level. Values are locus Sn / Pr (%), predicted genes.

| Genome | tx_strand | cds_strand | cds_blind |
|---|---|---|---|
| N. crassa OR74A | 63.8 / 71.8 (3447) | 63.5 / 70.9 (3480) | 63.5 / 70.8 (3478) |
| A. nidulans FGSC A4 | 54.7 / 55.7 (4911) | 54.7 / 55.1 (4975) | 54.7 / 55.1 (4975) |
| B. cinerea B05.10 | 82.9 / 81.4 (5711) | 83.2 / 80.3 (5811) | 83.2 / 80.2 (5811) |
| C. neoformans H99 | 80.5 / 81.4 (2691) | 79.8 / 78.8 (2754) | 79.9 / 79.0 (2753) |
| S. commune H4-8 | 25.0 / 36.0 (4984) | 25.1 / 35.6 (5060) | 25.1 / 35.6 (5059) |

### Table M2. Training-set construction, same fixed EVM evidence (training effect only)

Locus Sn / Pr (%); training models in parentheses. old = funannotate 41a2fd7; oldnew = old getBestModel + new selection; tx = new getBestModel (tx_strand, complete-first) + new selection; busco = BUSCO-forced training.

| Genome | old | oldnew | tx | busco |
|---|---|---|---|---|
| N. crassa OR74A | 63.3 / 71.6 (418) | 63.5 / 71.8 (418) | 64.3 / 72.6 (465) | 66.0 / 74.1 |
| A. nidulans FGSC A4 | 54.6 / 55.7 (1103) | 54.7 / 55.8 (1103) | 54.6 / 55.8 (1122) | 54.8 / 56.8 |
| B. cinerea B05.10 | 83.3 / 82.2 (1229) | 83.1 / 82.0 (1229) | 83.3 / 82.2 (1284) | 82.0 / 82.0 |
| C. neoformans H99 | 80.1 / 80.9 (814) | 80.1 / 81.0 (814) | 79.9 / 80.8 (832) | 78.1 / 80.6 |
| S. commune H4-8 | 29.2 / 44.2 (552) | 24.5 / 35.5 (87) | 24.9 / 36.0 (92) | 37.4 / 49.2 |

### Table M3. Ranking inside a locus: complete-first vs guarded (complete only if CDS ≥ 80% of the locus's longest)

| Genome | Evidence | complete-first locus Sn / Pr (train models) | guarded locus Sn / Pr (train models) |
|---|---|---|---|
| N. crassa OR74A | fixed | 64.3 / 72.6 (465) | 63.9 / 71.9 (441) |
| N. crassa OR74A | own | 63.8 / 71.8 (465) | 63.6 / 71.6 (441) |
| A. nidulans FGSC A4 | fixed | 54.6 / 55.8 (1122) | 54.6 / 55.8 (1111) |
| A. nidulans FGSC A4 | own | 54.7 / 55.7 (1122) | 54.7 / 55.8 (1111) |
| B. cinerea B05.10 | fixed | 83.3 / 82.2 (1284) | 83.3 / 82.1 (1251) |
| B. cinerea B05.10 | own | 82.9 / 81.4 (1284) | 83.3 / 82.0 (1251) |

### Table M4. Single-exon training genes (R6 option b): se (on) vs tx2 (off), same code, fixed evidence

Exact CDS-chain match to RefSeq protein-coding mRNAs on holdout chromosomes. Sn per RefSeq gene, Pr per predicted model.

| Genome | RefSeq single / multi genes | Single Sn | Single Pr | Multi Sn | Multi Pr |
|---|---|---|---|---|---|
| N. crassa OR74A | 877 / 3035 | 49.4 → 54.7 | 65.9 → 60.6 | 59.9 → 59.0 | 65.4 → 67.8 |
| A. nidulans FGSC A4 | 661 / 4353 | 71.3 → 78.5 | 58.6 → 55.6 | 47.7 → 47.4 | 50.6 → 51.8 |
| B. cinerea B05.10 | 1236 / 4484 | 64.2 → 71.8 | 74.0 → 73.2 | 79.5 → 79.2 | 77.2 → 79.0 |
| C. neoformans H99 | 77 / 2649 | 50.6 → 55.8 | 34.5 → 32.8 | 74.2 → 74.2 | 76.2 → 76.6 |
| S. commune H4-8 | 1237 / 5988 | 7.6 → 7.3 | 29.0 → 28.3 | 23.7 → 23.8 | 30.5 → 30.6 |

### Table M5. PASA training vs BUSCO training, divergent reads (N. crassa OR74A; RNA-seq from strain HJDF, median read identity 96.7%)

| Arm | PASA input | Single-exon training | Training | Fixed evidence locus Sn / Pr | Own evidence locus Sn / Pr |
|---|---|---|---|---|---|
| tx | rc1 (old minimap2 parser) | no | PASA | 64.3 / 72.6 | 63.8 / 71.8 |
| se | rc1 | yes | PASA | 65.1 / 73.6 | n/a |
| txR1 | R1 | no | PASA | 64.4 / 72.7 | 65.4 / 73.1 |
| txR1R2 | R1 + R2 | no | PASA | 64.6 / 73.0 | 65.8 / 73.6 |
| txID90 | R1 + gmap, 90% identity | no | PASA | 64.5 / 72.6 | 64.6 / 71.9 |
| txRel | R1, relaxed validation | no | PASA | 64.5 / 72.6 | 65.8 / 72.9 |
| seR1R2 | R1 + R2 | yes | PASA | 66.4 / 74.9 | 67.1 / 75.0 |
| busco | rc1 (evidence only) | n/a | BUSCO | 66.0 / 74.1 | n/a |
| buscoR1R2 | R1 + R2 (evidence only) | n/a | BUSCO | 66.1 / 74.2 | 66.8 / 74.4 |

### Table M6. Genome with few complete PASA models (S. commune H4-8: 93 complete of 822 PASA models, 30 keepers; divergent reads)

| Arm | Locus Sn / Pr | Intron-chain Sn / Pr | Exon Sn / Pr | Proteome BUSCO C (%) |
|---|---|---|---|---|
| busco | 37.4 / 49.2 | 42.3 / 48.6 | 69.8 / 79.6 | 97.5 |
| old | 29.2 / 44.2 | 34.0 / 43.8 | 56.4 / 78.1 | 86.9 |
| txR1R2 | 28.8 / 39.2 | 32.5 / 39.0 | 54.2 / 75.5 | 85.0 |
| tx | 24.9 / 36.0 | 28.1 / 36.2 | 47.0 / 73.8 | 79.4 |

### Table M7. Training set against RefSeq: old pipeline vs new selection (no filterGeneMark keeper set; upper bound, see text)

| Genome | Old: models, exact % | New: models, exact % | Exact RefSeq genes old → new | Redundant old → new |
|---|---|---|---|---|
| N. crassa OR74A | 2889, 38.9 | 1515, 78.7 | 1124 → 1192 | 46 → 0 |
| A. nidulans FGSC A4 | 5007, 50.1 | 4123, 62.1 | 2508 → 2559 | 37 → 4 |
| B. cinerea B05.10 | 5637, 75.5 | 4866, 89.4 | 4255 → 4348 | 59 → 0 |

### Table M8. getBestModel ranking against RefSeq (PASA models passed to EVM): exact CDS / exact intron chain

| PASA source | Genome | complete-first | structure-first | guarded |
|---|---|---|---|---|
| rc1 | N. crassa OR74A | 1549 / 2060 | 1530 / 2225 | 1536 / 2183 |
| rc1 | A. nidulans FGSC A4 | 3224 / 3210 | 3202 / 3221 | 3222 / 3221 |
| rc1 | B. cinerea B05.10 | 5455 / 5030 | 5419 / 5147 | 5447 / 5113 |
| R1 | N. crassa OR74A | 2207 / 3020 | 2174 / 3255 | 2187 / 3206 |
| R1 | A. nidulans FGSC A4 | 3258 / 3234 | 3238 / 3248 | 3258 / 3249 |
| R1 | B. cinerea B05.10 | 5767 / 5398 | 5726 / 5515 | 5761 / 5475 |
| R1R2 | N. crassa OR74A | 2303 / 3060 | 2275 / 3259 | 2290 / 3229 |
| R1R2 | A. nidulans FGSC A4 | 3259 / 3235 | 3239 / 3248 | 3259 / 3249 |
