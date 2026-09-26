# Transcript evidence, alignment and the value of RNA-seq in funannotate (draft methods and results text)

Status: draft for the funannotate paper, 2026-09-26. It is the companion to `training_data_selection_methods.md`, which covers training-data selection (gates, R3/R5, single-exon training). Owner: the REVIEW session. Every number comes from the files or DECISIONS entries named in brackets. This file records measured results and marks each conclusion as measured, inferred or not tested.

## 1. Summary of conclusions

1. **Training Augustus/SNAP from RNA-seq-derived PASA ORFs is rarely better than training from BUSCO genes by more than about 1 point, and it is clearly worse when there are few complete models.** *Measured; interim (3 genomes). Experiment A, D95/D97.*
   - Below about 300 complete PASA training models, BUSCO training gives 1-5 points higher holdout locus F1 in all three genomes tested.
   - Above about 500 models the two sources are within ±1 point. Botrytis is the exception, where PASA training is ahead by +0.7 to +1.7.
   - On a genome with 93 complete models (S. commune), BUSCO training beat the previous production path by +8.2 locus Sn / +5.0 Pr, and PASA with R1 + R2 by +8.6 / +10.0 [D79, Table M6].
   - Most of BUSCO's edge comes from single-exon genes, which PASA-derived training sets lacked until single-exon training (R6 b) was added [D73, D80].
2. **The main value of RNA-seq is as evidence for the final gene models, mainly through PASA models in EVM.** *Measured.*
   - On N. crassa, the alignment fixes R1 + R2 changed holdout locus Sn/Pr by +2.0/+1.8 when both the training set and the EVM evidence changed (txR1R2 own minus tx own).
     - Training alone: +0.3/+0.4 (txR1R2 fixed minus tx fixed).
     - Evidence alone, within txR1R2: +1.2/+0.6 (own minus fixed); by difference about +1.7/+1.4, if the effects add [scorecard.tsv, Table M5, D73].
   - Botrytis (same-strain reads) shows the same pattern: +0.7/+0.5 combined (txR1 own minus tx own); training alone −0.2/−0.3; evidence alone +0.5/0.0.
   - EVM weights explain why: PASA models carry weight 6, against 1 for Augustus, GeneMark, SNAP, proteins and transcripts, and 2 for Augustus HiQ [funannotate-predict.log, EVM Weights line].
3. **RNA-seq alignments also feed Augustus hints** (exon and intron hints from `bam2ExonsHints`). *Not measured separately.* No hints-on vs hints-off arm was run, so the share of the gain that comes through hints is unknown.
4. **GeneMark did not use RNA-seq in these experiments.** *Measured (configuration).* The arms ran GeneMark-ES (ab initio) from a precomputed `genemark.genome.gtf`. GeneMark-ET/EP modes with RNA-seq or protein hints were not tested.
5. **Practical conclusion:**
   - Keep RNA-seq in every run, for evidence.
   - Train from PASA only when there are enough complete models, which is what the gate controls. BUSCO training is a safe default below that threshold. The threshold is being calibrated (Section 6).

## 2. A long-standing error in converting minimap2 transcript alignments

- **What the code did** (funannotate `library.bam2gff3` / `bam2ExonsHints`, logic from 2018, commits `a9ed8f2`, `bc661b4`, `b77692b`; `ce48682` in 2023 only reformatted it). The code walked the minimap2 `cs` tag and:
  - did not advance the genome position on substitutions or deletions, so exon boundaries after the first mismatch were shifted;
  - counted soft-clipped bases as aligned;
  - checked intron motifs against the SAM flag instead of either orientation, which dropped reverse-oriented spliced contigs from unstranded libraries.
- **Effect in PASA validation** (spliced minimap2 alignments that pass) [D13, D24, D44]:

  | Genome (RNA-seq) | Before the fix | After the fix |
  |---|---|---|
  | N. crassa OR74A (divergent reads, about 96-97% identity) | 55 / 6,973 (0.8%) | 9,045 / 14,215 (64%) |
  | Botrytis cinerea B05.10 (same strain) | 999 / 8,076 (12%) | 13,406 / 15,192 (88%) |
  | A. nidulans FGSC A4 (same strain) | 2,910 / 3,869 (75%) | 6,855 / 7,247 (95%) |

  - Before the fix, about 70% of the alignments PASA accepted in N. crassa were single-exon. The converter fed PASA almost only single-exon alignments.
- **The fix:** coordinates come from the CIGAR string; clips are excluded from the Target range; consensus intron motifs are accepted in either orientation; the strand is the alignment strand for PASA import and the motif strand for EVM evidence and Augustus hints [PR nextgenusfs/funannotate#1210, merged; v1.9.0-rc.2/rc.3].
- **Exact RefSeq CDS chains in the PASA training models, before → after:** N. crassa 1,481 → 2,096 (+42%); Botrytis 5,353 → 5,644 (+5.4%); A. nidulans 3,171 → 3,220 (+1.5%) [D24, D44].

## 3. Splice-site accuracy by aligner (divergent reads)

- **Introns in spliced alignments that exactly match a RefSeq intron**, N. crassa [D36, D37; `code_review_20260925/intron_discordance.py`]:

  | Aligner | Exact | At 95-97% alignment identity | Main error type |
  |---|---|---|---|
  | minimap2 (fixed conversion) | 89.5% | 88% | 1.2% of introns shifted or non-canonical |
  | gmap | 61.0% | 59% | junctions slid within repeats (4,910) |
  | blat | 60.8% | 57% | non-canonical junctions (10,856) |
- **Non-matching minimap2 introns** are mostly inside RefSeq genes but different from the annotated intron (6.3% of all introns), or outside annotated genes (3.0%). Those are likely alternative isoforms or unannotated loci, not alignment errors.
- **Aligner choice in PASA:** gmap as PASA's own aligner added nothing over blat and cost 67% more runtime. gmap alone gave the weakest training set (1,172 exact chains, against 2,064 for fixed minimap2 alone) [D51].
- **gmap crash:** gmap (2021-12-17, 2023-04-28 and 2025-07-31) segfaults on some tandem-repeat contigs. PASA now retries in chunks and skips the offending transcripts [PASA v2.6.1-rc.2]. The root cause and a patch are documented [D56].

## 4. Relaxed PASA validation does not improve predictions

- Relaxing PASA's splice-boundary rule and identity threshold (`--pasa_num_bp_splice 0 --pasa_min_avg_per_id 90`) raised valid spliced alignments on N. crassa from 64% to 88%, and exact chains in the training models by 26% [D65].
- Holdout predictions did not clearly improve: txRel own 65.8/72.9 against txR1 own 65.4/73.1, that is +0.4 Sn / −0.2 Pr [D72, Table M5]. (Corrected: an earlier version compared against txR1R2, which also includes R2.)
- **Lesson:** gains measured on training models do not reliably carry through to prediction accuracy. Defaults were chosen from prediction-level scores.

## 5. A PASA output bug in production (F1)

- **The bug:** a fork commit (`bce776a`, 2026-06-29) printed every assembly GFF3 row twice. TransDecoder's genome mapping then gave frame-broken CDSs.
- **Production scan** of 8,007 PASA-trained BFD genomes [D48; `production_f1_scan.tsv`]:
  - the share of models whose CDS length is not divisible by 3 had a median of about 0.45 from 2026-07-06 to 09-10;
  - 5,599 of the 6,891 genomes in the window were affected; 0 of 1,116 genomes before it were.
- **Examples:** A. niger CBS 101883 and C. neoformans H99 had same-strain reads, yet only 74 and 160 complete PASA models. After the fix, H99 had 2,002 complete models of 2,704.
- **Fixed in:** PASA v2.6.1-rc.1 and later.

## 6. Calibrating the complete-model gate (in progress)

- **Design:** experiment A titrates complete PASA training models (N = 50-2,000; 5 draws for N ≤ 500, 3 above) within 4 RefSeq genomes, against BUSCO training on the same fixed evidence. Experiment B tests the crossover across about 40 RefSeq genomes [D81].
- **Decision rule (conservative):** the threshold is the smallest N at which the lower 95% bound of (PASA-trained − BUSCO-trained) holdout locus F1 is at least 0.
- **Interim, 3 genomes:** BUSCO training wins below about 300 models, and the difference is under 1 point above 500. Under the conservative rule, N* is 1,000 for Botrytis and not reached for A. nidulans or N. crassa [D95].
  - Caveat: the interim comparator used an older code snapshot [D97]. Repeats on matching code, H99 and experiment B are pending.
- **Gate variable:** the count of complete models tracks accuracy slightly better than the count of filterGeneMark keepers (Spearman 0.84-0.95 against 0.80-0.88).

## 7. Read identity of production RNA-seq

- **Median transcript-to-genome identity** across 8,007 PASA-trained BFD genomes: 69% at ≥ 99%, 31% at 90-99%, 0.6% below 90% [D52; `production_identity.tsv`].
- Transcript identity reads 0.4-0.9 points below read identity, because of Trinity assembly errors (Botrytis 99.58% against 100%). So the true divergent share is about 20-31%.
- Species-level RNA-seq selection often pairs reads from another strain with the reference genome. N. crassa's RNA-seq came from wild isolate HJDF, not OR74A [D19].

## 8. Not tested / open

- The contribution of Augustus hints (no hints-on vs hints-off arm).
- GeneMark-ET/EP with RNA-seq or protein hints.
- Whether EVM weights (PASA = 6) are optimal after the fixes. Weights were not refit on held-out data.
- The gate threshold is pending experiment A (final) and experiment B.
