# Decision log: PASA / funannotate training-set work (from 2026-09-25)

This log is shared by two Claude sessions working for jstajich:
- **REVIEW**: the PASA-side review session. It owns PASApipeline, `bam2gff3`, R13 and the scorers. Its report is `PASApipeline/CODE_REVIEW_20260925.md`.
- **SELECT**: the selection session. It owns the funannotate-live selection code (R3/R5, the locus rule) and `refseq_benchmark/benchmark.py`.

**Rules:**
- Add entries only. Do not rewrite old entries.
- To change a decision, add a new entry that supersedes the old one, and set the old one's status to `superseded by Dnn`.
- Each entry states its evidence as measured numbers with file paths, or says "no data".

**Status values:**
- `proposed`: not yet approved by the user.
- `accepted`: approved by the user.
- `done`: implemented, with its branch or commit named.
- `rejected`.
- `superseded`.

**Code locations:**

| Repo | Worktree / branch | Owner |
|---|---|---|
| funannotate-live | main working tree (uncommitted), base `41a2fd7` | SELECT |
| funannotate-live | `~/projects/funannotate/funannotate-live-bam2gff3`, branch `fix/bam2gff3-cigar` | REVIEW |
| PASApipeline | `rust_optimize` working tree (review files only) | REVIEW |
| PASApipeline | worktree for R2, branch `fix/unspliced-orient-clustering` (planned) | REVIEW |
| `refseq_benchmark/` | `benchmark.py`, `run.sh`, `library_old.py`, `train_old.py` | SELECT |
| `refseq_benchmark/` | `diversity.py`, `predict_scorer.py` | REVIEW |

**Merge order (proposed):** SELECT's R3/R5 first, then REVIEW's R1 rebased on top. Both change `library.py`, but in different functions. Before any merge, a Fable 5.1 subagent reviews the diff adversarially.

---

### D01: Rewrite minimap2 SAM→GFF3 conversion from the CIGAR string (R1)
- **Owner / status:** REVIEW / done on branch `fix/bam2gff3-cigar` (uncommitted); **needs user acceptance to merge**.
- **Decision:** replace the cs-string walk in `bam2gff3` and `bam2ExonsHints` with one parser, `parse_minimap2_splice_record()`. It works as follows:
  - Coordinates come from the CIGAR: M/=/X consume both sequences, D the genome only, N is an intron, I and S the transcript only.
  - The Target range excludes soft clips.
  - Every intron must be a consensus motif, and all introns must have the same orientation.
  - For the PASA GFF3, column 7 is the SAM strand. For the Augustus hints, a spliced alignment takes its strand from the intron motif.
- **Evidence:**
  - Tests: `tests/test_bam2gff3.py` has 18 tests; 11 failed on the old code and all 18 pass on the new code.
  - Real minimap2 check: the old code emitted 2 of 4 transcripts with shifted exons. The new code emits 4 of 4 with exact exons.
  - Real PASA data (Cordyceps, old run), from `alignment.validations.output`:
    - Spliced custom alignments: 7,409 of 7,746 (96%) failed validation, mostly "Splice site validations failed".
    - Single-exon custom alignments: 7,947 of 8,617 passed.
    - 4,627 transcripts failed as custom alignments but passed as PASA's own minimap2 and blat alignments.
- **Alternatives considered:**
  - Patch the cs walk. Rejected: the CIGAR is the authoritative source of coordinates.
  - Use `ts:A` for column 7. Rejected: PASA needs the alignment strand there. Fable tested this, and the ts:A version failed the emulated PASA check.
- **Code:** `funannotate/library.py`: `parse_minimap2_splice_record`, `bam2gff3`, `bam2ExonsHints`; `tests/test_bam2gff3.py`.

### D02: Mark the duplicate-GFF3-row bug (F1) as resolved
- **Owner / status:** REVIEW / done (user confirmed 2026-09-25).
- **Evidence:** rc.1 runs have 0 duplicate rows in `pasa_assemblies.gff3`. The bug was fixed in PASApipeline `4376a22`.

### D03: The rc.1 rust/norust comparison does not test the Rust assembler (F5)
- **Owner / status:** REVIEW / done. The user decided not to test F5, and `pasa_rust` will be removed from `bin/`.
- **Evidence:**
  - The rc.1 image does not set `PASA_ASSEMBLER`, and `PASA_alignment_assembler.pm:84-86` picks C++ `pasa` first.
  - With IDs removed, the assembly GFF3 coordinates are identical between the arms for N. crassa and A. nidulans.

### D04: Do not try to rescue incomplete models by short genome extension
- **Owner / status:** REVIEW / proposed.
- **Decision:** a 6-9 codon genome extension does not rescue incomplete models, so R3 drops incomplete models from training.
- **Evidence (N. crassa rc.1 against RefSeq):**
  - Only 2-4% of open ends are within 9 codons of the reference start or stop.
  - The median gap is 146-435 codons.
  - An in-frame stop occurs within 9 codons downstream in only 4-5% of cases, against 26-29% in the shifted-frame baseline.
  - Scripts: `r13_genome_extension.py`, `r13_reference_compare.py`.
- **Alternative still open:** recover the missing parts upstream in PASA (see D06).

### D05: Shared evaluation protocol
- **Owner / status:** REVIEW + SELECT / accepted by both sessions; user approved the collaboration 2026-09-25.
- **Training-set level:** `benchmark.score()` (SELECT) plus `diversity.diversity()` (REVIEW), against RefSeq.
- **Prediction level:** `predict_scorer.py` (REVIEW).
  - It keeps protein-coding CDS only, uses gffcompare 0.12.10, and scores Sn/Pr at base, exon, intron, intron-chain, transcript and locus level.
  - It uses a chromosome holdout: chromosomes of at least 1 Mb, sorted by length, alternate between training and holdout. Training sets are filtered to the training chromosomes before Augustus/SNAP training.
- **Genomes:** N. crassa, A. nidulans, Botrytis B05.10.
- **Selection rule:** choose from the scorecard, not by argument. The SELECT session confirms compute with the user before launching the predict arms.

### D06: Investigate evidence lost inside PASA before accepting a single-exon cap
- **Owner / status:** REVIEW / in progress (R13 rerun job 29106979).
- **Evidence (N. crassa rc.1):** in 30.8% of incomplete models, the minimap2 Trinity alignments together cover at least 90% of the RefSeq gene, but the PASA assemblies together do not. The combined Cordyceps evidence in D01 suggests that the custom-alignment validation failures are a major cause.
- **Next step:** in the rerun, cross-tabulate each lost gene's transcripts against valid/failed status per aligner.

### D07: Evidence-quality gates in train and predict
- **Owner / status:** SELECT / done in the funannotate-live main working tree (uncommitted); Nextflow wiring in `Fungi_BFD/nextflow` (uncommitted). Needs user acceptance to commit.
- **Decision:**
  - `funannotate train --min_rnaseq_map_rate` (default 10%): maps the first 200k reads with minimap2 splice:sr; if too few map, it exits 3 before Trinity/PASA.
  - `funannotate predict --min_pasa_complete_models` (default 500): counts complete ORFs in the PASA GFF3; if too few, Augustus/SNAP train from BUSCO.
  - Nextflow passes these flags only when a param is set (older images reject unknown flags). Exit 3 is routed to the existing `.pasa_train_failed` degrade path, and train logs are copied to `genome_annotation/<out>/logfiles/`.
- **Evidence:** 18-genome pilot (`do_annotation_triagePASArerun/qc_busco_compare/`).
  - Failing genomes mapped 0.1-2% of reads and had 4-281 complete PASA models; Meyerozyma mapped 97.3% with 1,258.
  - Production sweep (`qc_training_sweep_20260925/`): 2,572 rerun candidates; 55 of 1,439 species map below 10%.
- **Alternatives:** a Trinity transcript-alignment-rate gate (not needed; the read rate is more direct).
- **Code:** `library.py` `rnaseq_concordance_gate`, `sample_read_map_rate`, `run_rnaseq_concordance_gate`, `count_complete_orf_models`, `pasa_training_gate`, `run_pasa_training_gate`; `train.py` main (after read normalization); `predict.py` training-source block. Tests: `tests/test_training_gates.py`.
- **Caveat:** both thresholds come from about 12 comparable genomes. They should be re-derived from the predict-arm scorecard (D11).

### D08: R3, train only on complete ORFs, judged from the protein
- **Owner / status:** SELECT / done (uncommitted); Fable 5.1 review pending.
- **Decision:** `selectTrainingModels` drops every model that is not a complete ORF before any other filter.
  - Completeness is tested by `library.is_complete_model()`: M start, `*` at the end, no internal `*`, codon_start 1.
- **Evidence:**
  - First version tested `gff2dict` `cds_transcript`, which is in genome orientation for minus-strand genes. It rejected 848 of 848 minus-strand Meyerozyma models; the true number of complete ones is 628.
  - The protein test finds 1,235 complete models (607 plus, 628 minus), which equals the coordinate-based `count_complete_orf_models`.
  - REVIEW (D04): 96% of incomplete N. crassa models are in-frame fragments of real genes; extending their ends rescues 2-4%.
- **Alternatives:** use TransDecoder's `complete` label (equivalent on N. crassa: 2,146 = 2,146; the protein test does not need a second file).
- **Code:** `library.is_complete_model`, `library.selectTrainingModels`. Tests: `tests/test_training_selection.py`, including a minus-strand regression test.

### D09: R5, one model per locus in getBestModel, ranked by structure before TPM
- **Owner / status:** SELECT / done (uncommitted); the locus rule is pending D10.
- **Decision:** cluster models transitively (union-find). Keep one per cluster, ranked by complete ORF, then CDS exon count, then CDS length, then TPM, then ID.
  - Replaces a one-sided, strand-blind overlap test that kept every model winning its own window, with TPM ranked first.
- **Evidence (Meyerozyma rc.1):**
  - Old code: 1,640 models; new (tx_strand): 2,503.
  - 920 of the new models are at loci where the old code kept no model (CDS disjoint from every old model), and 839 of them are complete.
  - Only 2 of the 2,503 overlap an opposite-strand CDS, so there is no antisense flood.
- **Code:** `train.getBestModel`, `library.cluster_overlapping`.

### D10: Locus rule for getBestModel, chosen by RefSeq benchmark
- **Owner / status:** SELECT / proposed. The benchmark is job 29107015 → `refseq_benchmark/benchmark.tsv`.
- **Options:**
  - **tx_strand** (current): transcript span, same strand, ≥30% of the shorter model. In Meyerozyma it still merges 1,283 CDS-disjoint loci and loses 2,894 CDS loci.
  - **cds_strand:** CDS span, same strand.
  - **cds_blind:** CDS span, either strand.
- **Decide on:** exact_cds, refseq_genes_exact against the step1 ceiling (N. crassa 1,584, per REVIEW), redundant_models (target 0), no_refseq_overlap, and then the D11 prediction scores.

### D11: Predict-arm experiment (training-set variant → prediction accuracy)
- **Owner / status:** SELECT / accepted by the user 2026-09-25 ("use slurm as much as you need"). Running.
- **Design:** see `REPORT_training_selection.md` section D3.
  - rc.1 image for every arm; new code bind-mounted from the frozen snapshot `predict_arms/code_new/` (diff md5 75bec9b5…).
  - Step A trains on train chromosomes only (REVIEW's `predict_scorer.py split`).
  - Step B predicts the full genome with `-p parameters.json` and evidence identical across arms.
  - Scoring: `predict_scorer.py score` on holdout chromosomes, plus BUSCO.
- **Arms:** old; oldPASA+newSelect; tx_strand; cds_strand; cds_blind; busco (forced BUSCO training). Plus production-like step B runs for tx/cds with their own PASA evidence.
- **Genomes:** N. crassa OR74A, A. nidulans FGSC A4, Botrytis cinerea B05.10.
- **Jobs:** prep 29107057 (N. crassa, A. nidulans) + 29107060 (Botrytis, after PASA job 29107008).

### D12: A/B test of R1 (bam2gff3 fix) inside the unchanged rc.1 image
- **Owner / status:** REVIEW / running. Started 2026-09-25 as jobs 29107189 (R1 fix, N. crassa and A. nidulans) and 29107193 (A. nidulans baseline with KEEP_R13).
- **Decision:** test R1 by binding a patched `library.py` over the image copy, instead of rebuilding the image.
  - `run_train_compare.sh` has a new opt-in variable `LIBPY_OVERRIDE`. The variable `RESROOT` redirects output (`r1_fix/`, `r13_rerun/`). `KEEP_R13=1` keeps the PASA intermediate files. With none of these set, the script behaves exactly as before.
  - The run records the md5 of the override file in `inputs.tsv`.
- **Why this is valid:** the image's `site-packages/funannotate/library.py` has the same md5 (574dc94b…) as funannotate-live `41a2fd7:funannotate/library.py`, which is the base of branch `fix/bam2gff3-cigar`. So the only difference between the arms is the R1 change.
- **Frozen input:** `do_pasa_rust_vs_perl/r1_fix_library.py`, a copy of the worktree file taken at submission.
- **Comparison:**
  - Baselines: `r13_rerun/<genome>.rust` against `r1_fix/<genome>.rust`.
  - Metrics:
    - validation outcome by aligner and segment count (`r13_lost_evidence.py`);
    - single-exon share of valid alignments and of assemblies (M2);
    - `benchmark.score` + `diversity` of `funannotate_train.pasa.gff3` against RefSeq.

### D13: Finding: PASA's valid alignments are about 70% single-exon; the imported minimap2 channel gives almost no spliced evidence
- **Owner / status:** REVIEW / measured; this is not a code decision.
- **Evidence** (N. crassa rc.1, `r13_rerun/Neurospora_crassa_OR74A.rust/r13/alignment.validations.output.gz`, summarized in `r13_lost_evidence.py` output):

  | Aligner | Spliced alignments valid | Single-exon alignments valid |
  |---|---|---|
  | custom (minimap2 via old bam2gff3) | 55 of 6,973 (0.8%) | 8,432 of 13,713 |
  | blat | 6,193 of 18,320 (34%) | 5,910 of 10,119 |

  - Top custom failure: "Splice site validations failed" (6,725).
  - Top blat failure: "Incontiguous alignment + splice site failed" (7,100).
- **Genes whose Trinity alignments cover ≥90% of the RefSeq gene but whose PASA assemblies do not (947):**
  - 65% lose coverage at validation.
  - 35% (329) have valid alignments that cover the gene; the coverage is lost later, in clustering or assembly. R2 and F4 target this group.
- **R13 member test (training models, current code):**
  - 88% of incomplete models come from fragment transcripts.
  - 8.7% have a member transcript with a complete ORF. But the median ORF length difference is -36 aa, so that is usually a different ORF.
  - Alignment truncation is therefore a minor cause of missing start/stop codons.
- **Consequence:** R1 is expected to change the spliced/unspliced balance of valid evidence, so it has the highest priority on the PASA side. D12 will measure it.

### D14: Experiment: gmap instead of blat as PASA's own aligner, together with R1
- **Owner / status:** REVIEW / running. Job submitted 2026-09-25. Output: `r1_fix_gmap/Neurospora_crassa_OR74A.rust`.
- **Why:** D13 shows blat validates only 34% of spliced Trinity alignments. The top failure is "Incontiguous alignment + splice site failed".
- **Arm:** R1 `library.py` override plus `EXTRA_TRAIN_ARGS="--aligners minimap2 gmap"`. funannotate removes minimap2 from PASA's list, so PASA runs gmap and imports the fixed minimap2 GFF3 as custom alignments.
- **Script change:** `run_train_compare.sh` has a new opt-in variable `EXTRA_TRAIN_ARGS`, appended to the `funannotate train` command. It is empty by default.
- **Compare against:** `r1_fix/Neurospora_crassa_OR74A.rust` (R1 with blat) and `r13_rerun/...` (baseline). Same metrics as D12.

### D15: User decisions, 2026-09-25
- **D04 → accepted.** Drop incomplete models from training (R3). There will be no genome-extension or rescue step.
- **R2 and F4 → accepted as opt-in flags.**
  - New `Launch_PASA_pipeline.pl` options, off by default, so default PASA output stays upstream-compatible. funannotate turns them on.
  - R2 lets `?` single-exon alignments join the spliced cluster or subcluster that contains them.
  - F4 keeps one valid alignment per transcript per cluster when more than one aligner is used.
  - Code: PASApipeline worktree `../PASApipeline-r2`, branch `fix/unspliced-orient-clustering`. A Fable review is required before merge.
- **Aligner default → accepted, conditional on data.** Propose changing funannotate train's default `--aligners` from `minimap2 blat` to `minimap2 gmap` only if the D14 scorecard shows gmap is better on RefSeq metrics. Report the runtime cost next to the accuracy gain.
- **Still waiting for data before the user decides:**
  - D01: merging R1 (needs the D12 A/B results and the Fable review).
  - D10: the locus rule (SELECT's benchmark).
  - D14: the aligner result.

### D17: D08 (R3, complete-ORF filter, no rescue) accepted
- **Owner / status:** SELECT / accepted (user decision recorded in D15, item 1). Code is done in the funannotate-live main working tree and not yet committed. It still needs the Fable 5.1 adversarial review (protocol) before merge.
- **Scope:** only the training set (`selectTrainingModels`) drops incomplete models. The PASA GFF3 passed to EVM as evidence is unchanged by R3.

### D16: R1 changes after the Fable 5.1 pre-merge review
- **Owner / status:** REVIEW / done on branch `fix/bam2gff3-cigar` (uncommitted). The D12/D14 jobs were cancelled and resubmitted with this version.
- **Review verdict:** APPROVE WITH CHANGES. The review ran 18+15 tests and an old-vs-new comparison on a real 34,362-record BAM.
  - Kept records: old 16,150, new 21,845. Runtime: 10.3 s old, 2.9 s new.
  - The motif-derived strand agrees with minimap2's ts:A tag for 13,370 of 13,379 records.
- **Changes made:**
  1. **HIGH:** `bam2gff3` also feeds predict's `harmonize_transcripts`, which produces EVM evidence and Augustus hints. So there is a new parameter, `strand="align"` (the default, for PASA in train.py and update.py) or `strand="splice"` (motif strand, now passed by `harmonize_transcripts`). Without this, 27.6% of kept records (antisense-oriented spliced contigs) would have given evidence on the wrong strand.
  2. **MEDIUM:** identity is now gap-compressed: substitutions + indel events. This is closer to the old values and to blat/psl per-id, which PASA compares with MIN_AVG_PER_ID. Counting indel bases would have lowered the identity of 71% of records.
  3. **LOW:** the flag filter is now `0xF04`, which also drops QC-fail (0x200) and duplicate (0x400) records.
  4. **LOW:** a record is skipped when its cs tag and CIGAR disagree on the number of introns.
- **Accepted as intended behavior changes:**
  - The Target range excludes soft clips, so PASA's MIN_PERCENT_ALIGNED now rejects heavily clipped transcripts.
  - Any non-consensus intron rejects the record. The old code checked only the last intron; 3.4% of old-kept records drop, and PASA would reject them anyway.
- **Tests:** `tests/test_bam2gff3.py` now has 23 tests, all passing. Existing test files pass. The real minimap2 check gives 4 of 4 transcripts with exact coordinates, and splice mode gives the correct gene strand.
- **Frozen A/B input:** `r1_fix_library.py`, md5 dc7ab460af04….

### D18: R2 + F4 opt-in PASA flags: implementation and end-to-end test
- **Owner / status:** REVIEW / implemented on PASApipeline branch `fix/unspliced-orient-clustering` (worktree `../PASApipeline-r2`, uncommitted). Test jobs submitted 2026-09-25: 29107369 (R2 only) and 29107370 (R1 + R2), N. crassa rust arm.
- **Code:**
  - `PerlLib/Unspliced_orient_join.pm` (new). A `?` alignment gets an orientation only when all spliced alignments covering at least L% of its span share that orientation. It stays `?` if both orientations or none cover it.
  - `PerlLib/t/unspliced_orient_join.t`: 14 tests, all passing.
  - `assign_clusters_by_stringent_alignment_overlap.dbi -U`: uses the stringent threshold L.
  - `subcluster_builder.dbi -U`: uses its -m threshold, default 50.
  - `ensure_single_valid_alignment_per_cdna_per_cluster.pl -S`: prefers a spliced alignment, then the higher score.
  - `Launch_PASA_pipeline.pl`: new flags `--UNSPLICED_JOIN_SPLICED` (adds -U to stringent clustering and subclustering) and `--ONE_ALIGNMENT_PER_CDNA` (runs the one-per-transcript step with -S when more than one aligner is used, including custom). Both are off by default, per D15.
  - Not covered: `--gene_overlap` clustering (only a notice is printed; subclustering still applies) and `sharding/shard_subcluster_builder.sh` (does not pass -U yet).
  - Perl syntax check passes in the rc.1 image.
- **Why the test is valid:** the image's `Launch_PASA_pipeline.pl`, the three changed scripts and all 96 `PerlLib` files are byte-identical to PASApipeline `23d67c0`, the base of the branch. The worktree files are frozen in `r2_pasa/` and bound over the image copies.
- **Test harness only:** funannotate has no way to pass extra PASA flags, so `r2_test_train.py` is the image's `train.py` (identical to funannotate `41a2fd7`) plus two lines that append `$PASA_EXTRA_ARGS` to the Launch command. It belongs to no branch. A real funannotate change to pass the flags is still to be designed.
- **Script change:** `run_train_compare.sh` has a new opt-in variable `EXTRA_BINDS`, which is recorded in `inputs.tsv` together with `PASA_EXTRA_ARGS`.
- **Compare:** `r13_rerun` (baseline) against `r1_fix` against `r2_only` against `r1_r2`, with the D12 metrics plus assembly-level single-exon share and overlap with spliced assemblies (M2).
- **Numbering note:** entries are not in numeric order in this file. Always take the highest D number: `grep -o '^### D[0-9]*' DECISIONS.md | sed 's/### D//' | sort -n | tail -1`.

### D19: Finding: the N. crassa RNA-seq is from a divergent wild isolate, so the benchmark must be read per genome
- **Owner / status:** REVIEW / measured. Proposes a user decision on a read-identity gate (see below).
- **Evidence:**
  - The reads come from SRA SRR33994706, SRR33994703 and SRR33994701 (`get_sra/rnaseq_fix_reps_20260830.samples.rnaseq_sra.csv`). They were selected by species taxid 5141.
  - The BioSample for SRR33994706 says: strain "Neurospora crassa HJDF", isolation source "red mold tofu" (Fujian). This is not OR74A.
  - Identity of blat alignments to the OR74A genome, from `alignment.validations.output`:

    | Genome | Share of alignments by identity |
    |---|---|
    | N. crassa | 95-98%: 73.9%; 90-95%: 18.7%; ≥99.5%: 0.6% |
    | A. nidulans | ≥99.5%: 98.8% |
  - Spliced alignments that pass validation: custom 0.8% and blat 34% for N. crassa, against 75% and 89% for A. nidulans.
- **Interpretation:**
  - With 2-5% divergence, the old `bam2gff3` shift bug (D01) breaks almost every spliced minimap2 alignment.
  - PASA's MIN_AVG_PER_ID 95 removes about 19% of alignments.
  - So N. crassa tests divergent reads, a realistic case for species-level SRA selection in BFD. A. nidulans tests same-strain reads.
  - The existing RNA-seq gate measures mapping rate only: 95.8% for N. crassa, which hides the divergence.
- **Proposed follow-up (needs user decision):**
  1. Report the scorecard per genome, and label N. crassa as "divergent RNA-seq (~96-97% identity)".
  2. Add read or transcript identity to the RNA-seq gate: for example the median identity of the sampled read alignments (NM/aligned length) or of the transcript alignments. Warn or branch when it is below about 99%. No threshold has been tested.
  3. For divergent data, consider a lower `--pasa_min_avg_per_id` (for example 90). Untested; it would need its own A/B run.

### D20: Botrytis B05.10 reads are same-strain; its read file is stale relative to its query cache
- **Owner / status:** SELECT / done (measurement). Follows REVIEW's D19 (N. crassa reads are strain HJDF, not OR74A).
- **Evidence:**
  - `rnaseq_reads/Botrytis_cinerea_norm_R1.fastq.gz` contains SRR31907356 (1,430,664 reads; a dsRNA nano-carrier study, strain not given) and SRR27676449 (1,642,487 reads).
  - Neither run is among the 5 in `sra_query/Botrytis_cinerea.sra_query.csv`, so the read file predates the current cache (the same pattern as Ramularia ERR1530623).
  - minimap2 splice:sr against B05.10, 50k reads per run: median NM-based identity is 100.0% for both runs (10th percentile 98.7%).
- **Benchmark labels:** A. nidulans FGSC A4 = clean same-strain; Botrytis B05.10 = clean (identity); N. crassa OR74A = divergent reads (strain HJDF, about 95-98% identity). Results are reported per genome, never averaged across genomes.

### D21: Proposed: add a read-identity check to the RNA-seq concordance gate
- **Owner / status:** SELECT / proposed (REVIEW suggestion). **Needs user decision.** Not implemented.
- **Proposal:** `sample_read_map_rate` already aligns 200k reads. It could also report the median NM-based identity of the mapped reads, and train would warn (not fail) below about 99%.
- **Evidence:** the mapping-rate gate passed N. crassa at 95.8% while its reads come from another strain (about 95-98% identity). Botrytis and A. nidulans have about 99.5-100% identity.
- **Open question:** whether to warn only, or to also select relaxed PASA thresholds automatically. There is no data yet on which is better; the N. crassa predict arms (D11) will show how much divergent reads cost.

### D22: Cross-strain aligner comparison (N. crassa = divergent reads) and a crash-safe gmap wrapper
- **Owner / status:** REVIEW / running. Jobs submitted 2026-09-25: 29107491-29107496. Output: `do_pasa_rust_vs_perl/xs/<arm>/<genome>.rust`.
- **User direction (2026-09-25):**
  - Cross-strain RNA-seq should still work, so N. crassa becomes a test category.
  - When reads are about 95% identical, try gmap instead of minimap2.
  - Keep fixing the minimap2 parsing (R1), and compare the two.
- **Arms** (N. crassa; A. nidulans as the clean control where marked `+An`):

  | Arm | PASA alignment evidence |
  |---|---|
  | baseline (`r13_rerun`) | old minimap2 custom + blat |
  | `r1_fix` | fixed minimap2 + blat |
  | `xs/mm2fix_only` +An | fixed minimap2 only (PASA's own aligner skipped) |
  | `xs/gmap_only` +An | gmap only (custom import skipped) |
  | `xs/blat_only` +An | blat only |
  | `xs/mm2fix_gmap` +An | fixed minimap2 + gmap |
  | `xs/mm2old_gmap` | old minimap2 + gmap (isolates the aligner effect without R1) |
  | `xs/mm2fix_gmap_id90` | fixed minimap2 + gmap, `--pasa_min_avg_per_id 90` |
- **gmap crash found:**
  - The first gmap attempt (29107331) died. gmap segfaults (SIGSEGV) on transcript `Trinity_GG_5605_c0_g1_i2`, a 290-bp contig with a tandem repeat.
  - Reproduced with gmap 2025-07-31 (the image version), 2023-04-28 and 2021-12-17, under every option set tried (-B 2..5, with and without -x 50, -n 0/1, --no-chimeras).
  - So one bad contig kills PASA's whole gmap step.
- **Fix:** `scripts/process_GMAP_alignments_gff3_chimeras_ok.pl` on branch `fix/unspliced-orient-clustering`.
  - gmap now writes to a temp file. On success the file is copied to STDOUT, so the content is the same; line order already varies between gmap -t runs.
  - Only on failure, and only for GFF3 output, it reruns in 500-sequence chunks. Failing chunks are split down to single transcripts; crashing transcripts are skipped and listed in `gmap.failed_transcripts.txt`. SAM mode still dies on error.
  - Tested in the image on 300 real contigs plus the crashing one: exit 0, 297 aligned (the same set as a gmap run without the bad contig), and only the bad contig skipped.
  - This is a crash fix and is active by default: it changes nothing unless gmap fails. **Needs user acknowledgement**, because D15 asked for opt-in PASA changes.
- **Test harness only:** `r3_test_train.py` is the image's `train.py` plus `PASA_EXTRA_ARGS`, `SKIP_CUSTOM_ALIGNMENTS=1` and `SKIP_PASA_ALIGNERS=1` switches, which edit the Launch command. It belongs to no branch.
- **Metrics:**
  - Validation by aligner and segment count.
  - Single-exon share of valid alignments and of assemblies.
  - Training set: `benchmark.score` + `diversity` against RefSeq.
  - Runtime (`run_status.tsv`), as D15 requires for the aligner decision.

### D23: RefSeq benchmark result; D10 locus rule stays tx_strand pending the predict arms
- **Owner / status:** SELECT / done (benchmark); the rule choice is **proposed**, pending D11 own-evidence arms. Supersedes nothing; updates D10.
- **Data:** `refseq_benchmark/benchmark.tsv` (job 29107015). Reads: N. crassa divergent (D19); A. nidulans and Botrytis clean (D20 above).
- **Training set (new select, any locus rule) vs old pipeline (old getBestModel + old select):**

  | Genome | Old: models, exact % | New: models, exact % | Exact RefSeq genes, old → new |
  |---|---|---|---|
  | N. crassa | 2,889, 38.9% | 1,515, 78.7% | 1,124 → 1,192 |
  | A. nidulans | 5,007, 50.1% | 4,123, 62.1% | 2,508 → 2,559 |
  | Botrytis | 5,637, 75.5% | 4,866, 89.4% | 4,255 → 4,348 |

  Redundant models: 46/37/59 → 0/4/0. The three locus rules give the same training set, because selectTrainingModels then clusters on any overlap. So the locus rule affects only the EVM evidence.
- **getBestModel output (EVM evidence), Botrytis as the example:**

  | Rule | Models | Exact | Exact % | No RefSeq overlap | Redundant | Antisense-only |
  |---|---|---|---|---|---|---|
  | old | 8,290 | 5,348 | 64.5 | 444 | 399 | 44 |
  | tx_strand | 8,423 | 5,455 | 64.8 | 501 | 369 | 76 |
  | cds_strand | 9,180 | 5,723 | 62.3 | 764 | 469 | 90 |
  | cds_blind | 9,151 | 5,723 | 62.5 | 764 | 469 | 63 |

  A. nidulans has the same pattern (tx 3,224 exact, 46.4%; cds_blind 3,295, 44.6%). On N. crassa all three new rules tie (about 1,550 exact).
- **Conclusions:**
  - cds_strand is dominated by cds_blind (same exact matches, fewer antisense), so it is dropped.
  - tx_strand against cds_blind is a recall/precision trade-off at the evidence level. The D11 arms B.own (tx vs cdsB) decide.
  - Until then tx_strand stays as the working-tree default.
- **Side finding (for R6):** every training set has 0% single-exon genes, in the old code and the new. The multi-CDS requirement removes them once at least 200 multi-exon models exist. RefSeq single-CDS mRNAs: N. crassa 21.4%. Augustus/SNAP never learn single-exon genes. **Needs a user decision (R6), no data yet on the effect.**

### D24: Result of the R1 A/B test (D12, D16): large gain on divergent reads, small gain on same-strain reads
- **Owner / status:** REVIEW / measured. This supports merging R1 (D01); merge still needs user approval.
- **Runs:** baseline `r13_rerun/<g>.rust` against `r1_fix/<g>.rust`. Same rc.1 image; only `library.py` differs (md5 dc7ab460…).
- **Validation of the custom (minimap2) alignments,** from `alignment.validations.output`:

  | Genome | Spliced valid | Spliced failed | Single-exon valid |
  |---|---|---|---|
  | N. crassa: baseline → R1 | 55 → 9,045 | 6,918 → 5,170 | 8,432 → 8,253 |
  | A. nidulans: baseline → R1 | 2,910 → 6,855 | 959 → 392 | 3,193 → 3,178 |

  - The blat numbers are identical between arms, as expected.
- **Training models against RefSeq** (`funannotate_train.pasa.gff3`, image's old selection, `code_review_20260925/training_set_vs_refseq.py`):

  | Genome | Models | Single-CDS | Complete | Exact chains | RefSeq genes hit | Exon-count KS |
  |---|---|---|---|---|---|---|
  | N. crassa: baseline | 6,817 | 55.5% | 31.6% | 1,481 | 5,086 | 0.341 |
  | N. crassa: R1 | 8,024 | 46.3% | 35.1% | **2,096 (+42%)** | 5,941 | 0.250 |
  | A. nidulans: baseline | 6,903 | 24.1% | 75.0% | 3,171 | 6,241 | 0.103 |
  | A. nidulans: R1 | 6,948 | 23.9% | 76.0% | 3,220 (+1.5%) | 6,282 | 0.101 |
- **Runtime:** N. crassa 1,266 → 1,520 s; A. nidulans 905 → 1,015 s. PASA has more valid alignments to process.
- **For SELECT:** `r1_fix/<g>.rust/pasa.step1.gff3` and `kallisto.tsv` can feed the new selection, to cross selection rule × R1.

### D25: Selection × R1 cross, evidence level
- **Owner / status:** SELECT / running (job 29107591).
- **Design:** the same `refseq_benchmark/benchmark.py` (old / tx_strand / cds_strand / cds_blind getBestModel, then new or old selection), run on REVIEW's R1-fixed PASA output (`r1_fix/<genome>.rust/pasa.step1.gff3` + `kallisto.tsv`) via the env override `BENCH_RES`.
- **Output:** `refseq_benchmark_r1/benchmark.tsv`, comparable row for row with `refseq_benchmark/benchmark.tsv`.
- **Genomes:** N. crassa and A. nidulans (job 29107591); all three genomes rerun after REVIEW Botrytis R1 job 29107590 (dependent job below).
- **Next:** after this cross, add predict arms (D11 design) for the best selection variant on R1 PASA input.

### D26: R2 fix: subcluster_builder -U must store the assigned orientation; R2 arms rerun
- **Owner / status:** REVIEW / fixed on branch `fix/unspliced-orient-clustering`. Arms resubmitted as 29107594 (R2 only) and 29107595 (R1 + R2), now on N. crassa and A. nidulans.
- **Failure in the first R2 runs (29107369/70):**
  - `classify_alt_splice_isoforms_per_subcluster.dbi` (the `--ALT_SPLICE` step) takes the orientation of a subcluster member.
  - With `-U`, a subcluster could hold a `?` assembly next to `+` assemblies, so it tried to insert `orient='?'` into `splice_variation`, which fails the CHECK constraint (`orient IN ('', '+','-')`). PASA died.
- **Fix:** when `subcluster_builder.dbi -U` assigns a `?` assembly an orientation, it now writes that orientation to `align_link.spliced_orient` for the assembly and sets it on the object, and it logs the change.
  - This is safe for the loader: `Ath1_cdnas.pm` raises an error only when the computed orientation is `+` or `-` and differs from the stored one. A single-exon assembly computes `?`, so the stored orientation is accepted.
- **Frozen file refreshed:** `r2_pasa/scripts/subcluster_builder.dbi`. The other R2 files are unchanged.
- **Lesson:** Perl syntax checks and unit tests did not catch this. The end-to-end run did. Every PASA-side change needs a full run before merge.

### D27: User accepts the gmap crash fallback as default behavior (not opt-in)
- **Owner / status:** REVIEW (code owner) / accepted by the user 2026-09-25, relayed through SELECT. User's words: "I think this change is okay to be default rather than opt-in since it solves crash problems."
- **Scope:** the chunked, bisecting rerun in `scripts/process_GMAP_alignments_gff3_chimeras_ok.pl` (PASApipeline worktree `~/projects/funannotate/PASApipeline-r2`, branch `fix/unspliced-orient-clustering`, uncommitted). It stays on by default. It changes nothing unless gmap fails. This is an exception to D15's opt-in rule, which still applies to R2 and F4 (`--UNSPLICED_JOIN_SPLICED`, `--ONE_ALIGNMENT_PER_CDNA`).

### D28: Status update for the gmap crash fallback (D22)
- **Owner / status:** REVIEW / accepted as default behavior, following D27. The approval was relayed by the SELECT session, quoting the user; REVIEW asked the user to confirm directly.
- **What this means:** the chunked/bisection fallback in `process_GMAP_alignments_gff3_chimeras_ok.pl` stays on by default, as an exception to the D15 opt-in rule. It acts only when gmap fails.
- R2 (`--UNSPLICED_JOIN_SPLICED`) and F4 (`--ONE_ALIGNMENT_PER_CDNA`) stay opt-in.

### D29: Correction to D23: in a real predict run, R3 changes the training set only when filterGeneMark keeps fewer than 200 models
- **Owner / status:** SELECT / done (measurement). Qualifies D23 and D08 (R3 is still accepted, but it matters in fewer genomes than D23 implied).
- **Why D23 overstated it:** the benchmark ran selectTrainingModels with an EMPTY keeper GTF (no predict hints at that stage), so keeperCheck was False in every variant. A real predict run has hints.ALL.gff, and filterGenemark.pl keeps only models whose introns match RNA-seq/protein hints. That already removes most incomplete models.
- **Evidence (predict arms, step A, real hints):**
  - A. nidulans old.A vs oldnew.A: 3,668 PASA genes → 1,377 filterGeneMark keepers (1,107 multi-CDS) → 1,103 training models in BOTH arms. Only 12 CDS lines differ between the two `final_training_models.gff3` files.
  - N. crassa old.A: 657 keepers (418 multi-CDS) → 418 training models.
- **Where R3 matters:** when keepers < 200, the old code falls back to all multi-CDS PASA models, partial ones included. Examples: the pilot Colletotrichum (96 keepers) and Drepanopeziza, the genomes that lost BUSCO.
  - Production sample: of 400 random PASA-trained genomes (`qc_training_sweep_20260925/stageA.tsv`), **95 (24%) had fewer than 200 keepers** and took the fallback path; 305 had at least 200.
- **Consequence for D11:** on clean genomes with good hints (A. nidulans), expect old ≈ oldnew in the predict scores. The selection effect should appear mainly on low-keeper genomes. The next predict arms should include a low-keeper genome (e.g. a production genome from the 95).

### D30: Add two low-keeper production genomes (good reads, collapsed PASA training) to the experiment
- **Owner / status:** SELECT / running.
- **Why:** D29 shows R3 only matters when filterGeneMark keeps fewer than 200 models. In production, several genomes have well-matched reads but almost no usable PASA training models (read map % from `qc_training_sweep_20260925/stageB.tsv`; keepers from their production `funannotate-predict.log`):

  | Genome | Reads mapped | Keepers | Complete PASA | Final genes |
  |---|---|---|---|---|
  | A. niger CBS 101883 | 96.6% | 36 | 74 | 992 |
  | Aureobasidium pullulans EXF-150 | 92.2% | 9 | 16 | 2,817 |
  | Penicillium antarcticum IBT 31339 | 99.4% | 2 | 4 | 11,038 |
  | Cryptococcus neoformans H99 | 98.8% | 74 | 160 | 6,982 |

  So in these genomes PASA itself fails despite matched reads. This is consistent with REVIEW's D24 (R1: 0.8% of spliced minimap2 alignments valid before the fix).
- **Added (RefSeq, chromosome-level or large scaffolds):** Cryptococcus neoformans H99 (GCF_000149245.1, 6,975 coding genes) and Schizophyllum commune H4-8 (GCF_000143185.2, 16,187). They are cases 6 and 7 in `cases.tsv`.
- **Jobs:**
  - PASA baseline 29107599; PASA with R1 (`RESROOT=r1_fix LIBPY_OVERRIDE=r1_fix_library.py`) 29107600.
  - Benchmarks: `refseq_benchmark_lowkeep/` (29107601) and `refseq_benchmark_lowkeep_r1/` (29107602).
  - Predict prep 29107605; 30 arm jobs in `predict_arms/jobs_lowkeep.tsv`.

### D31: gmap work scope (user choice 2026-09-25, asked directly by REVIEW)
- **Owner / status:** REVIEW / accepted by the user; in progress.
- **The user selected all four options:**
  1. **Find the gmap segfault root cause.**
     - Build gmap from source with debug symbols, get a backtrace on the crashing contig `Trinity_GG_5605_c0_g1_i2` (290 bp, tandem repeat, N. crassa OR74A), and find the faulting code.
     - Patch it in a separate build directory (`/bigdata/stajichlab/jstajich/projects/funannotate/gmap_debug`).
     - Verify: the contig no longer crashes, and the output on the other contigs is unchanged.
     - Then decide on an upstream report to the GMAP/GSNAP author.
     - The user also relayed through SELECT: "fix the gmap executable to avoid this crash bug too".
     - The PASA chunked fallback (D22/D27/D28) stays as a safety net.
  2. **Make gmap a robust PASA aligner:** index reuse, threads, and tests for the gmap GFF3 import (chimeras, multiple paths).
  3. **gmap in funannotate's own transcript mapping** (the `--IMPORT_CUSTOM_ALIGNMENTS` input), for example when the identity gate detects divergent reads.
  4. **Wait for the aligner results:** default changes and item 3's routing are decided from the D22 scorecard (`xs/` arms).
- **Order:** item 1 now (a Fable 5.1 subagent does the C debugging); item 2 after the fallback results; item 3 after the D22 scorecard.

### D32: Residual minimap2 failures after R1, and a parser fix for insertions next to introns (R1 v2)
- **Owner / status:** REVIEW / fixed on branch `fix/bam2gff3-cigar` (uncommitted).
  - Frozen as `r1_fix_library_v2.py` (md5 f8d7c5326755…).
  - The running arms keep v1 (`r1_fix_library.py`, dc7ab460). v1 is not overwritten.
- **User question:** do minimap2 splice-site and orientation errors persist?
- **Answer (spliced custom alignments, `alignment.validations.output`):**

  | Failure reason | N. crassa: baseline → v1 | A. nidulans: baseline → v1 |
  |---|---|---|
  | "Splice site validations failed" (the parser's orientation and position errors) | 6,725 → 109 | 525 → 41 |
  | Mismatch within 3 bp of a splice boundary (NUM_BP_PERFECT_SPLICE_BOUNDARY=3) | 183 → 2,665 | 433 → 203 |
  | "Only N % is aligned" (soft clips now counted) | 0 → 1,132 | 0 → 8 |
  | Identity below 95% | 10 → 990 | 1 → 3 |
  | "Incontiguous alignment" | 0 → 274 | 0 → 137 |
  | **Valid** | **0.8% → 63.6%** | **75.2% → 94.6%** |
- **Interpretation:**
  - The orientation and splice-site parsing errors are gone.
  - The N. crassa boundary, identity and aligned-share failures reflect divergent reads (D19); they are not parser errors.
  - "Incontiguous" was a v1 parser artifact. All 137 of 137 A. nidulans cases have an insertion immediately after an intron (CIGAR `…N<k>I…`). v1 left those inserted bases out of both exons, creating a gap in the transcript coordinates.
- **v2 fix:** bases inserted immediately after an intron now start the next exon. A new test covers this: `test_insertion_right_after_intron_keeps_target_contiguous` (24 tests, all passing).
- **Verified on the real BAMs:** spliced records with a transcript-coordinate gap drop from 137 to 0 (A. nidulans) and from 274 to 0 (N. crassa). The record counts do not change (10,448; 27,948).
- **Next:** the final scorecard reruns use v2 plus the gmap fallback. The current v1 arms differ from v2 only in these at most 137/274 alignments.

### D33: R1 merge path (user decision 2026-09-25, option b)
- **Owner / status:** REVIEW / PR open, merge pending.
- **Decision:** push R1 now and open a PR. Merge it after SELECT commits R3/R5 to `target_1.9/rust_EVM_trinity_PASA`, keeping the agreed order.
- **Done:**
  - branch `fix/bam2gff3-cigar` (commit `c134412`) pushed to nextgenusfs/funannotate;
  - PR https://github.com/nextgenusfs/funannotate/pull/1210 against `target_1.9/rust_EVM_trinity_PASA`.
- **PASApipeline:** review files are commit `a3478a1` on `rust_optimize`. R2, F4 and the gmap fallback are branch `fix/unspliced-orient-clustering` (commit `04b63cc`). Both are pushed to hyphaltip/PASApipeline.
- **Note:** the first push/PR attempt was blocked by the Claude Code auto-mode classifier. The user then authorized the push and the PR explicitly.

### D34: Commit of SELECT's gates + R3/R5 is pending the user's direct approval and the Fable review
- **Owner / status:** SELECT / proposed (awaiting the user).
- REVIEW reports the user chose the merge order "SELECT commits R3/R5 first, then R1 (branch fix/bam2gff3-cigar, commit c134412) on top". SELECT will not commit on a relayed decision; it needs the user's instruction in the SELECT session.
- The Fable 5.1 adversarial review of the gates + R3/R5 diff is running (report: `REVIEW_R3R5_fable.md`).
- **Commit scope when approved:**
  - funannotate-live: `library.py`, `train.py`, `predict.py`, `CHANGELOG.md`, `tests/test_training_gates.py`, `tests/test_training_selection.py`;
  - Nextflow (separate repo, `Fungi_BFD/nextflow`): FUNANNOTATE_TRAIN, FUNANNOTATE_PREDICT, `profile_funannotate.config`.
- The D10 locus rule may still change after the B.own predict arms. Committing tx_strand now and changing it later would be a small follow-up commit.

### D35: User decision: commit only after the prediction results settle the locus rule
- **Owner / status:** SELECT / accepted (the user, directly in the SELECT session, 2026-09-26: "wait for prediction results to settle locus rule first").
- **Plan:**
  1. Finish the D11 predict arms (fixed and own evidence) and the low-keeper arms (D30).
  2. Choose tx_strand or cds_blind from the holdout gffcompare and BUSCO scores.
  3. Apply the Fable review findings (`REVIEW_R3R5_fable.md`).
  4. Then ask the user once more before committing.
- R1 (REVIEW) stays unmerged until then, per the agreed order.

### D36: Finding: fixed minimap2 gives accurate splice sites on divergent reads; PASA validation is the bottleneck. Relaxed-validation arms started.
- **Owner / status:** REVIEW / measured; arms running (29107658, 29107659).
- **User direction:** "if minimap2 can be wired to give accurate splice site and alignment info even with some mismatches that would be ideal."
- **Intron accuracy on N. crassa (divergent reads).** Introns come from the spliced alignments in `alignment.validations.output`, compared with RefSeq exon-derived introns:

  | Aligner | Exact RefSeq intron | ≥99% id | 97-99% | 95-97% | <95% |
  |---|---|---|---|---|---|
  | fixed minimap2 (`xs/mm2fix_only`) | **89.5%** (23,609 / 26,368) | 95% | 92% | 88% | 81% |
  | gmap (`xs/gmap_only`) | 61.0% | 83% | 79% | 59% | 51% |
  | blat (`xs/blat_only`) | 60.8% | 86% | 70% | 57% | 35% |
  | old minimap2 parser (`r13_rerun`) | 2.2% | — | — | — | — |
- **Training models against RefSeq (N. crassa):**

  | PASA evidence | Exact chains | RefSeq genes hit | Single-CDS |
  |---|---|---|---|
  | fixed minimap2 + blat | 2,096 | 5,941 | 46.3% |
  | fixed minimap2 only | 2,064 | 5,867 | 46.4% |
  | blat only | 1,363 | 4,550 | 48.1% |
  | gmap only | 1,172 | 4,145 | 58.0% |

  - gmap validated only 28% of its spliced alignments, and its crash fallback skipped 1 transcript in the real run.
  - **gmap is not better for divergent reads in this test.** This contradicts the working hypothesis behind D22.
- **Remaining failures of fixed minimap2** are PASA's rule that the 3 bases at each junction match exactly (`--pasa_num_bp_splice`, default 3; strain SNPs) and `--pasa_min_avg_per_id` (default 95). The minimap2 command in these runs had no `--junc-bed`.
- **New arms:**
  - `xs/mm2fix_blat_relaxed`
  - `xs/mm2fix_only_relaxed`
  - Both use `--pasa_num_bp_splice 0 --pasa_min_avg_per_id 90`, on N. crassa and A. nidulans. A. nidulans is the control: does relaxing admit false alignments in clean data?

### D37: User decisions (grilling, 2026-09-25) and the intron-discordance result
- **Owner / status:** REVIEW / user decisions recorded; measurement done.
- **Identity-gate categories (accepted):**

  | Median read identity | Category | Action |
  |---|---|---|
  | ≥ 99% | same strain | current defaults |
  | 90-99% | divergent | adapt PASA |
  | < 90% | likely another species | skip PASA; train from BUSCO |

  - The user's reasoning: a low-identity transcript mapping may be less informative than BUSCO markers alone.
  - The cut points are judgment calls; no data supports them yet. The user asked whether the ~8,000 production genomes can supply the identity distribution. A scoping agent is running.
- **Divergent action (accepted, pending tests):** do not route to gmap. Keep fixed minimap2 (+ blat) and relax PASA validation (`--pasa_num_bp_splice`, `--pasa_min_avg_per_id`). Values come from the D36 relaxed arms. The user wants less focus on gmap and more on minimap2 splice-site correctness.
- **Exonerate est2genome polishing: not needed** (user asked; REVIEW measured). Non-matching introns of fixed minimap2 on N. crassa, 2,759 of 26,368 (script `code_review_20260925/intron_discordance.py`):

  | Class | Count | Share of all introns |
  |---|---|---|
  | inside a RefSeq gene, different intron | 1,654 | 6.3% |
  | outside RefSeq genes | 801 | 3.0% |
  | shifted 1-10 bp | 210 | 0.8% |
  | non-canonical | 94 | 0.4% |

  - So real minimap2 splice-placement errors are about 1.2% of introns.
  - For comparison: blat non-matching introns are mostly non-canonical (10,856); gmap's are mostly junctions slid within repeats (4,910) and shifts (1,943).
  - 8,397 exact-RefSeq introns sit in minimap2 alignments that PASA rejected. This is the target of the relaxed validation.
  - If shifts ever matter, the first option is funannotate's existing `--junc-bed` input. These runs did not use it.

### D38: R6 single-exon training genes: user chose to test option (b) against (a)
- **Owner / status:** decision by the user (asked by REVIEW). Implementation and the predict arm belong to SELECT (notified).
- **(a)** Current behavior: 0% single-exon genes in training.
- **(b)** Admit single-exon genes that are complete ORFs with protein-homology or BUSCO support, capped at the genome's own single-exon share. The share is estimated from the BUSCO/compleasm gene models funannotate already produces.
- **Evaluation:** a predict arm scored with `predict_scorer.py` on the holdout chromosomes. Compare single-exon sensitivity and precision, and check that multi-exon accuracy does not drop. Make (b) the default only if it helps.
- **Unverified premise:** how Augustus sets its single-exon parameters when it has no single-exon training genes. Nobody has checked this in the Augustus source.

### D39: User decisions: R2/F4 acceptance rule, and gmap work scope reduced
- **Owner / status:** REVIEW / user decisions (grilling, 2026-09-25).
- **R2/F4 acceptance rule (accepted).** Funannotate enables `--UNSPLICED_JOIN_SPLICED` / `--ONE_ALIGNMENT_PER_CDNA` by default only if, compared with R1 alone on both test genomes:
  1. exact RefSeq CDS chains in the training models do not drop, and precision drops by at most 1 percentage point;
  2. at least one of these improves: the single-exon share of PASA models (towards RefSeq), redundant models per RefSeq gene, or the share of single-exon assemblies inside spliced genes;
  3. PASA completes with `--ALT_SPLICE` on all test genomes (N. crassa, A. nidulans, Botrytis, H99, S. commune);
  4. if the two flags disagree, only the one that passes is enabled, and they are tested separately.
  - R2's added value is also measured on top of SELECT's R5, in the selection × PASA cross.
- **gmap root cause: option (a) accepted.** The Fable subagent finishes the root cause, the patch and a draft upstream bug report. No patched gmap goes into the container. The PASA crash fallback (D22/D27) stays.
- **"gmap in funannotate's own mapping" (D31 item 3):** REVIEW proposed dropping it, because the D36/D37 data show gmap is weaker than fixed minimap2 on divergent reads. **Not yet confirmed by the user.**
- **D31 item 2 ("gmap as a robust PASA aligner")** is deprioritized under the user's "less focus on gmap". Only the fallback remains.

### D40: User confirmed: drop "gmap in funannotate's own transcript mapping" (D31 item 3)
- **Owner / status:** REVIEW / decided by the user on 2026-09-25. This closes the open item in D39.
- **Reason:** fixed minimap2 is more accurate than gmap on divergent reads (D36, D37).
- The grilling round is complete: D33 and D37-D39 hold the decisions. The user confirmed the summary.

### D41: Fable 5.1 adversarial review of gates + R3/R5: findings and plan
- **Owner / status:** SELECT / review done; fixes in progress. Report: `REVIEW_R3R5_fable.md` (read-only review, no repo edits; 32/32 tests passed at review time).
- **Verified correct:** protein-based completeness on both strands and multi-mRNA index alignment; cluster_overlapping with nested/containment/unsorted/touching/reversed coordinates; exit 3 propagation; no pipe deadlock.
- **Fix now (test-first):**
  - (2, major, confirmed) selectTrainingModels exits with code 1 when R3 leaves 0 models (empty proteins FASTA → diamond makedb fails) instead of returning 0 so the BUSCO fallback can run.
  - (3, major, plausible) the predict gate runs only when augustus is in "pasa" mode. Gate when any predictor (augustus/snap/glimmerhmm) is in pasa mode.
  - (4, minor, confirmed) is_complete_model does not check CDS length % 3.
  - (8, minor) the gate runs before the augustus.gff3 checkpoint check. On a resume, or with --augustus_gff, it can flip snap to BUSCO. Skip the gate when Augustus output already exists or is supplied.
- **Decide by data, not argument:** (1, major, mechanism confirmed) getBestModel ranks completeness before exon count/length, so a long 5'-partial multi-exon model loses to a short complete single-exon model at the same locus. Whether that is better EVM evidence is empirical. Add ranking variants to the RefSeq benchmark (complete-first = current; structure-first; complete-first only when its CDS is ≥ 80% of the longest) and pick by exact RefSeq matches and redundancy.
- **Deferred minors:** (5) O(k²) in huge stacks; (6) sample_read_map_rate robustness (try/finally, samtools dependency check, FASTQ format assumptions); (7) gate placement after normalization and re-running on resume; (9) the 500 threshold was calibrated on old getBestModel output; (10) defensive checks on codon_start.
- **Effect on running predict arms:** they use the frozen snapshot `predict_arms/code_new/` (pre-fix). Fixes 2, 3, 4 and 8 do not change their outputs (those paths are not reached there). If fix 1 changes the getBestModel ranking, the PASA variant files would change and the affected arms would be rerun.

### D42: R6 option (b) implemented (opt-in) for the D38 test
- **Owner / status:** SELECT / done in the working tree (uncommitted); predict arm to follow.
- **Code:** `library.protein_single_exon_share`, `library.single_exon_supported`, and `selectTrainingModels(..., single_exon_align=None, single_exon_share=None)`. Tests: 8 new cases in `tests/test_training_selection.py` (share estimator; admission, zero cap, unsupported, opposite strand; default unchanged).
- **Deviation from D38 wording:** the single-exon share is estimated from exonerate protein2genome alignments (`predict_misc/protein_alignments.gff3` + `proteins.combined.fa`), not from BUSCO models. Reason: in the PASA-training path, predict produces no BUSCO gene models, but the protein alignments already exist before training.
  - Estimator: among alignments covering ≥90% of the protein, the fraction with no inter-segment gap ≥20 bp.
  - Support: same-strand protein segments cover ≥80% of the single CDS.
  - Cap: share/(1−share) × multi-exon models, rounded down; longest first.
  - Supported singles bypass the filterGeneMark keeper test, which measures intron support.
- **Caveat:** Swiss-Prot evidence contains curated proteins of these model species, so support is strongest for them (the same in all arms).

### D43: R2 + F4 end-to-end results against the D39 acceptance rule: pass on N. crassa and A. nidulans; three more genomes running
- **Owner / status:** REVIEW / criteria 1-2 met on both genomes. Criterion 3 is pending: R1+R2 on Botrytis, H99 and S. commune, job 29107674.
- **Runs:** `r1_fix` (R1) against `r1_r2` (R1 + `--UNSPLICED_JOIN_SPLICED --ONE_ALIGNMENT_PER_CDNA`). Both use the same library (v1, dc7ab460). The R2 files are frozen in `r2_pasa/` (with the D26 fix).
- **Training models against RefSeq** (`training_set_vs_refseq.py`):

  | Genome | Exact chains | Precision | Single-CDS | Redundant models | Exon-count KS |
  |---|---|---|---|---|---|
  | N. crassa: R1 → R1+R2 | 2,096 → 2,213 (+5.6%) | 26.1 → 28.1% | 46.3 → 44.4% (RefSeq 21.4) | 1,831 → 1,697 | 0.250 → 0.230 |
  | A. nidulans: R1 → R1+R2 | 3,220 → 3,221 | 46.3 → 46.4% | 23.9 → 23.9% | 325 → 324 | 0.101 → 0.101 |
- **Assemblies** (`m2_single_exon_overlap.py`, `pasa_assemblies_described.txt`):

  | Genome | Single-exon assemblies overlapping a spliced assembly | Assemblies listing one transcript twice |
  |---|---|---|
  | N. crassa: R1 → R1+R2 | 43.5% → 38.5% | 10,029 → 0 |
  | A. nidulans: R1 → R1+R2 | 26.5% → 25.5% | 9,548 → 0 |

  - F4 removes all duplicate-transcript assemblies. R1 had roughly doubled them, because it made far more custom alignments valid.
- **R2 without R1** (N. crassa): exact chains 1,481 → 1,540 against baseline.
- **--ALT_SPLICE** completed in all 4 runs, which confirms the D26 fix.
- **Criterion 4:** the flags were tested together only. They pass jointly, so separate runs are needed only if a later genome fails.

### D44: R1 on Botrytis (same-strain reads): the old parser also hurt clean data
- **Owner / status:** REVIEW / measured. Jobs 29107589 (baseline) and 29107590 (R1). This adds evidence for PR #1210.
- **Spliced custom alignments valid:** 999 / 8,076 (12%) → 13,406 / 15,192 (88%). Single-exon: 6,589 → 6,584.
- **Training models against RefSeq:**

  | Botrytis | Exact chains | Precision | RefSeq genes hit | Single-CDS (RefSeq 20.4%) | Exon-count KS |
  |---|---|---|---|---|---|
  | baseline | 5,353 | 64.6% | 7,379 | 29.8% | 0.094 |
  | R1 | **5,644 (+5.4%)** | 65.3% | 7,735 | 28.1% | 0.078 |
- **Interpretation:** the Botrytis reads are about 100% identical to the genome (SELECT's check), yet only 12% of spliced minimap2 alignments survived the old parser, against 75% in A. nidulans. So the old bug is not limited to divergent reads.
  - Likely contributors: the SAM-flag motif rule, which drops about half of unstranded contigs, and small indels or mismatches in the Trinity contigs.
  - The breakdown by cause was not measured.
- **R1 exact-chain gain across genomes:** N. crassa +42%, Botrytis +5.4%, A. nidulans +1.5%. Precision rises or stays flat in all three.

### D45: R6 option (b) single-exon share: estimate from GeneMark-ES models, not protein alignments (changes D38 wording; **needs user decision**)
- **Owner / status:** SELECT / proposed. Implemented in the working tree for the D38 test; **the user decides the production source.**
- **Why not BUSCO (D38 wording):** in the PASA-training path, predict produces no BUSCO gene models.
- **Why not protein alignments (first implementation, D42):** REVIEW flagged possible bias. Measured against RefSeq single-CDS genes on the train chromosomes:

  | Genome | RefSeq | GeneMark-ES | Protein, ≥90% coverage | Protein, both ends aligned |
  |---|---|---|---|---|
  | N. crassa | 22.0% | 25.5% | 11.8% | 13.4% |
  | A. nidulans | 14.0% | 20.5% | 11.5% | 12.6% |
  | Botrytis | 21.5% | 22.8% | 7.9% | 10.0% |

  Mean absolute error: GeneMark 3.8 points; protein 8.8 (≥90%) / 7.2 (both ends). The protein estimate is biased LOW (conserved Swiss-Prot hits favor multi-exon genes), the opposite of the inflation REVIEW expected. GeneMark errs slightly high. Only 3 genomes; the direction is consistent.
- **Code:** `library.model_single_exon_share` (GFF3 Parent= or GTF gene_id); predict `--training_single_exon` uses `predict_misc/genemark.evm.gff3` first, then protein alignments as fallback. Tests: `ModelSingleExonShareTests`. Support is still protein homology (`single_exon_supported`).
- **Arms:** frozen snapshot `predict_arms/code_new2/` (diff md5 883fca10…, includes the D41 fixes). Arms `se` (flag on) and `tx2` (flag off, same code), step A + step B fixed, on all five genomes. Job table: `predict_arms/jobs_se.tsv`. Compare se vs tx2: single-exon Sn/Pr on holdout chromosomes, multi-exon accuracy, BUSCO.

### D46: Read identity of test and collapsed genomes (supports the D37 identity-gate categories)
- **Owner / status:** SELECT / done (measurement). Job 29107637; file `qc_training_sweep_20260925/identity_check/identity.out`.
- **Method:** first 50k R1 reads, minimap2 splice:sr, primary MAPQ≥1; identity = 1 − NM/(M+I bases).

  | Genome | Median | 10th pct | ≥99% | D37 category |
  |---|---|---|---|---|
  | A. nidulans FGSC A4 (control) | 100.0% | 100.0% | 0.976 | same strain |
  | N. crassa OR74A (control, strain HJDF reads) | 96.7% | 91.0% | 0.126 | divergent |
  | C. neoformans H99 | 100.0% | 100.0% | 0.986 | same strain |
  | A. niger CBS 101883 | 100.0% | 98.0% | 0.853 | same strain |
  | A. pullulans EXF-150 | 98.0% | 94.1% | 0.322 | divergent |
  | S. commune H4-8 | 94.0% | 89.3% | 0.039 | divergent |
  | P. antarcticum IBT 31339 | 92.7% | 88.8% | 0.011 | divergent (near the 90% skip line) |

- **Conclusions:**
  - The median-identity measure agrees with REVIEW's blat measurement on N. crassa (95-98%).
  - The D30 collapsed genomes split: A. niger and H99 collapsed with same-strain reads, so divergence is not the cause and REVIEW is checking their PASA logs. A. pullulans, P. antarcticum and S. commune have divergent reads, where R1 plus relaxed validation (D36/D37) should help.
  - The benchmark set now covers same-strain clean (A. nidulans, Botrytis), same-strain collapsed (H99) and divergent collapsed (S. commune) genomes.

### D47: Rerun candidate list v2 includes F1-affected production PASA sets (9,374 genomes)
- **Owner / status:** SELECT / done (list). The rerun itself **needs a user decision** (after the D35 commit).
- **Input:** REVIEW's production F1 scan (`identity_production/production_f1_scan.tsv`, job 29107705; 8,010 training sets).
  - The share of models with CDS length not divisible by 3, by `funannotate_train.pasa.gff3` date: **before 2026-06-29: 1,116 sets, max 0.000**; F1 window (06-29 to 09-24): 6,891 sets, median 0.362, 81% > 0.01, 57% > 0.05. No sets are dated after 09-24.
  - Before-window sets are all exactly 0, so any nonzero share in the window is the F1 signal.
- **Criterion:** file dated in the F1 window AND frac_not_mod3 > 0.01 → 5,599 sets. All 5,599 were used by their own predict run (PASA evidence and/or training).
  - Plus 3,052 pretrained genomes inheriting ab-initio parameters from an F1-affected representative (`repr_assignments.tsv`).
- **Merged with the Stage A/B list (D27/D28-era `rerun_candidates.tsv`, 2,572):** 1,849 overlap with F1; 723 are flagged for other reasons.
  - **Total 9,374** in `qc_training_sweep_20260925/rerun_candidates_v2.tsv` (PASA 5,943, pretrained 3,101, BUSCO 315, none 15). The old file is kept unchanged.
- **Implication:** most of what the sweep flagged traces to F1, which rc.1 already fixes (`4376a22`). A rerun with the gate-aware, R3/R5 image plus fixed PASA is the remedy. The remaining 723 (and the 315 BUSCO-trained ones flagged only by the EVM loci ratio) need the EVM review.

### D48: F1 in production: the timeline, and why A. niger and H99 collapsed
- **Owner / status:** REVIEW / measured. SELECT merged the result into `rerun_candidates_v2.tsv` (D47). The rerun itself is a user decision.
- **A. niger CBS 101883 and C. neoformans H99** have same-strain reads (SELECT identity check: median 100%). Their production PASA models are F1-corrupted:
  - A. niger (file dated 2026-08-12): 44% of CDS lengths not a multiple of 3, 91% without a stop codon, 74 complete models out of 6,163.
  - H99 (2026-09-10): 60% not a multiple of 3, 92% without a stop, 160 complete out of 5,448.
  - rc.1 runs show 0% not a multiple of 3.
- **Production-wide scan** (`code_review_20260925/production_f1_scan.py`, 8,007 genomes, `identity_production/production_f1_scan.tsv`):
  - Clean (all genomes below 0.1% of CDSs not a multiple of 3) through 2026-07-04.
  - First affected genomes on 2026-07-06, the same day as PASApipeline commit `ce3bffc` ("fix double name issue"). This is probably when `bce776a` (2026-06-29) first reached production images. The image build was not verified.
  - Median affected share about 0.45 from 2026-07-06 to 09-10. 5,599 of 6,891 genomes in the window have more than 1%.
  - Before the window: 0 of 1,116.
- **Unexplained:** on 2026-09-11, 1,240 of 1,549 genomes show a milder 1-5% signature. No commit in PASApipeline or funannotate on 09-08..09-12 explains it (only funannotate `dde8243`, a docker base-image update). These genomes are still flagged at the 1% threshold. Their cause should be checked on one example before the rerun.

### D49: abu_dhabi nodes (c01-c30) break predict and train with SIGILL: jobs excluded and rerun
- **Owner / status:** SELECT / done.
- **What happened:** jobs submitted with `-p epyc` ran on abu_dhabi nodes c01/c02/c05 (they lack BMI2/AVX2).
  - bam2hints died with "Illegal instruction (core dumped)", exit 132 (reproduced by hand on c01) → 5 predict step B arms failed.
  - kallisto index died in funannotate train → both Schizophyllum PASA runs (29107599_12, 29107600_12) exited 1, while the SLURM job still said COMPLETED.
  - All completed arms ran on r-nodes, so no finished result is affected.
- **Actions:**
  - `#SBATCH --exclude=c[01-30]` added to `predict_arms/arm.sh` and `prep.sh`; `ExcNodeList=c[01-30]` set on all 50 pending jobs under the account (including REVIEW's 29107702).
  - Failed arms resubmitted (`predict_arms/jobs_resub.tsv`). Bad Schizophyllum and H99-R1 outputs moved to `*.bad_abu_dhabi_20260926`.
  - The low-keeper chain was resubmitted: PASA 29107713 and 29107714, benchmarks 29107715 and 29107716, prep 29107717, and new arm ids appended to `jobs_lowkeep.tsv` and `jobs_se.tsv`.
  - REVIEW was notified: its running 29107658/29107659 tasks are on c01/c03.
- **Other workflows currently on c-nodes** (not part of this work, not touched): e.g. `fabench_v1.8.17_conda` (c01), `ed-Mucoromycota`, `polishcap-*`, `serinales-all`. **The user may want to check them.**

### D50: Relaxed-validation arms resubmitted away from the abu_dhabi nodes
- **Owner / status:** REVIEW / resubmitted as 29107756 and 29107757 with `--exclude=c[01-30]`.
- **Why:** SELECT found that `-p epyc` jobs can land on the preemptible abu_dhabi nodes c01-c30. Those nodes lack BMI2/AVX2, so kallisto and bam2hints die with SIGILL, while `run_train_compare.sh` still reports COMPLETED.
  - The first relaxed arms (29107658/59) were running on c01/c03, so they were cancelled and their partial output removed.
- **Check of all earlier REVIEW runs:**
  - All ran on r-nodes (r23-r40), per `sacct`.
  - No capture log contains "Illegal instruction" or "CMD ERROR".
  - All exited 0.
  - So the earlier results (D12-D44) are not affected.

### D51: Aligner comparison complete: gmap adds nothing as a companion aligner; a 90% identity cutoff trades precision for recall
- **Owner / status:** REVIEW / measured (jobs 29107494-96; all on r-nodes, exit 0). This supports D37: no gmap routing, and `minimap2 blat` stays the default.
- **N. crassa (divergent), training models against RefSeq:**

  | PASA evidence | Exact chains | Precision | RefSeq genes hit | Single-CDS | Redundant | Runtime |
  |---|---|---|---|---|---|---|
  | fixed mm2 + blat | 2,096 | 26.1% | 5,941 | 46.3% | 1,831 | 1,520 s |
  | fixed mm2 + gmap | 2,096 | 26.3% | 5,919 | 46.4% | 1,806 | 2,539 s |
  | old mm2 + gmap | 1,220 | 20.6% | 4,574 | 63.5% | 1,125 | 1,996 s |
  | fixed mm2 + gmap, `--pasa_min_avg_per_id 90` | 2,301 | 23.5% | 6,717 | 51.6% | 2,682 | 2,628 s |
- **A. nidulans (clean):** fixed mm2 + gmap 3,238 exact chains against fixed mm2 + blat 3,220 (+0.6%).
- **Conclusions:**
  - gmap as the companion aligner gives no gain and costs 67% more runtime.
  - gmap does not rescue the old parser (1,220 against 1,481 with blat).
  - Lowering the identity cutoff to 90% on divergent reads recovers real genes (+205 exact chains, +776 RefSeq genes hit), but adds single-exon and redundant models (precision -2.8 points before selection).
  - The final value for the divergent category waits for the relaxed arms (29107756/57), and ideally for scoring after SELECT's R3/R5 selection.

### D52: Production identity distribution (8,007 PASA-trained genomes) against the D37 gate categories
- **Owner / status:** REVIEW / measured. Job 29107702. Output `identity_production/production_identity.tsv`. Script `code_review_20260925/production_identity.py`.
- **Source:** column 6 of `genome_annotation_training/<out>/training/transcript.alignments.gff3`, one value per Trinity transcript alignment. Two caveats:
  - These files come from the pre-fix bam2gff3, which dropped alignments below 80% identity and many spliced alignments. A median of 0.67 alignments per transcript survived (p10 0.51, p90 0.91).
  - This is transcript identity, not the read identity the gate is defined on.
- **Median transcript identity per genome:**

  | Median | Genomes | Gate category |
  |---|---|---|
  | ≥ 99.5% | 3,465 | same strain |
  | 99-99.5% | 2,045 | same strain |
  | 98-99% | 837 | divergent |
  | 97-98% | 560 | divergent |
  | 95-97% | 533 | divergent |
  | 90-95% | 521 | divergent |
  | < 90% | 46 | other species |

  - Totals: same strain 5,510 (69%), divergent 2,451 (31%), other species 46 (0.6%).
- **Calibration against read identity** (SELECT's `identity_check/identity.out`):

  | Genome | Read median | Transcript median |
  |---|---|---|
  | A. nidulans | 100% | 100.00% |
  | Botrytis | 100% | 99.58% (p10 98.0%) |
  | N. crassa | 96.7% | 95.85% |

  - Transcript identity is about 0.4-0.9 points lower than read identity, because of Trinity assembly errors.
- **Consequences:**
  - The 98-99% transcript band (837 genomes) is ambiguous; many are likely same-strain by read identity.
  - The best production estimate of truly divergent genomes is about 1,600-2,450 (20-31%).
  - The gate should keep using read identity, as SELECT implemented it, not transcript identity.
  - The below-90% group is small (46). Sending those to BUSCO training affects few genomes.
  - The 99% cut point is supported by the calibration: same-strain reads sit at a median of 100%. No change is proposed.

### D53: Predict arms with fixed PASA inputs (R1, R1+R2, relaxed identity) to settle "BUSCO vs PASA on divergent reads"
- **Owner / status:** SELECT / running. Requested by REVIEW, whose point is valid: the interim result "BUSCO-forced beats PASA-trained on N. crassa (locus 66.0/74.1 vs 64.3/72.6)" used the image's old bam2gff3, which left 55 of 6,973 spliced custom alignments valid.
- **Inputs:** getBestModel (tx_strand, frozen code_new, the same code as the 'tx' arm) on r1_fix/ (R1), r1_r2/ (R1+R2) and xs/mm2fix_gmap_id90/ (R1 + gmap + --pasa_min_avg_per_id 90) for N. crassa, and r1_fix/ for Botrytis. Job 29107761 → `predict_arms/pasa_inputs/<src>/`.
- **Arms:** txR1, txR1R2, txID90 (N. crassa) and txR1 (Botrytis). Each has step A, step B fixed evidence (isolates training) and step B own evidence (production-like). Job table: `predict_arms/jobs_pasafix.tsv`.
- **Decision this informs:** D37's divergent-read action (relaxed PASA vs BUSCO training). Compare with busco-forced on N. crassa (66.0/74.1) and with tx on Botrytis (83.3/82.2). REVIEW's better relaxed arms (29107756/57) will be added when they finish.

### D54: D10 locus rule settled on 3 genomes: keep tx_strand
- **Owner / status:** SELECT / settled for N. crassa, A. nidulans and Botrytis. Confirmation on the low-keeper genomes (H99, S. commune) is pending; only a reversal there would reopen it. Supersedes D10 and D23's "proposed".
- **Evidence (holdout chromosomes, gffcompare CDS level, step B with each rule's own PASA evidence; `predict_arms/scorecard.tsv`):**

  | Genome | tx_strand: locus Sn/Pr | cds_strand | cds_blind |
  |---|---|---|---|
  | N. crassa | 63.8 / 71.8 | 63.5 / 70.9 | 63.5 / 70.8 |
  | A. nidulans | 54.7 / 55.7 | 54.7 / 55.1 | (rerun pending) |
  | Botrytis | 82.9 / 81.4 | 83.2 / 80.3 | 83.2 / 80.2 |

  The CDS-span rules add about 100 predicted genes (mostly single-exon; Botrytis 19.3% → 20.8-20.9%), gain at most 0.3 points of locus sensitivity, and always lose precision (0.6-1.2 points). This matches the evidence-level benchmark (D23: cds_blind adds models with no RefSeq overlap).
- **Consequence for D35:** the locus-rule question the user asked to wait for is answered on three genomes. Commit preparation can start (Fable fixes are done, D41). The user is asked once more before committing.

### D55: getBestModel ranking: guarded complete-first (review finding 1 resolved by data)
- **Owner / status:** SELECT / done in the working tree (uncommitted); predict-arm confirmation running.
- **Rule (`train.pick_locus_model`):** a complete ORF wins a locus only if its CDS is ≥80% of the locus's longest CDS. Otherwise CDS exons, then CDS length, then TPM, then ID decide. Replaces unconditional complete-first (D09).
- **Evidence (`refseq_benchmark/rank_benchmark.tsv`, job 29107681; tx_strand loci):**

  | Source / genome | complete_first: exact CDS / exact intron chain | struct_first | guarded |
  |---|---|---|---|
  | rc1 N. crassa | 1,549 / 2,060 | 1,530 / 2,225 | 1,536 / 2,183 |
  | rc1 A. nidulans | 3,224 / 3,210 | 3,202 / 3,221 | 3,222 / 3,221 |
  | rc1 Botrytis | 5,455 / 5,030 | 5,419 / 5,147 | 5,447 / 5,113 |
  | R1 N. crassa | 2,207 / 3,020 | 2,174 / 3,255 | 2,187 / 3,206 |
  | R1 Botrytis | 5,767 / 5,398 | 5,726 / 5,515 | 5,761 / 5,475 |
  | R1R2 N. crassa | 2,303 / 3,060 | 2,275 / 3,259 | 2,290 / 3,229 |

  - complete_first keeps about 230 more "complete" N. crassa models, but only about 13 more exact RefSeq genes. The rest are mostly short internal ORFs (likely wrong starts) that would enter training via R3.
  - guarded keeps exact matches within 0-20 of complete_first and recovers most of struct_first's intron-chain gain.
- **Tests:** `GuardedRankingTests` (long partial beats short complete; near-full-length complete beats partial; TPM breaks ties last). Full suite passes (31 selection tests).
- **Confirmation arms:** txG (code_new3 snapshot, diff md5 913b342d…) with guarded getBestModel on rc1 PASA; step A + B fixed + B own on N. crassa, A. nidulans, Botrytis. Inputs job 29107776; arm ids in `predict_arms/jobs_pasafix.tsv`. Compare with the tx arms.

### D56: gmap segfault root cause found and patched (D31 item 1, D39 option a)
- **Owner / status:** REVIEW (Fable 5.1 subagent) / done. Materials are in `/bigdata/stajichlab/jstajich/projects/funannotate/gmap_debug/`: README.md, `gmap_segfault.patch`, `reproducer/`, `upstream_bug_report_draft.md`, and the gdb/ASan logs. Nothing is committed or sent. Per D39 (a), no patched gmap goes into the image.
- **Root cause** (the same code is in gmap 2021-12-17, 2023-04-28 and 2025-07-31):
  1. `Pairpool_clean_join` (`src/pairpool.c:688`): the query starts with a 93-bp tandem duplication, so two local alignments map to the same genomic interval. The overlap-fixing loop skips only gap holders and leaves indel pairs (INDEL_COMP/SHORTGAP_COMP) at the junction. Query base 143 then ends up neither clipped nor aligned, so the cigar length check fails (289 != 290).
  2. `Stage3_new` returns NULL, and `Stage3_merge_local` dereferences it without a check (`src/stage3.c:18127`), which gives the SIGSEGV.
  - The crash is not SIMD-dependent.
- **Patch:**
  - `pairpool.c`: the junction-clearing loops also pop indel pairs, using the rule `clean_end_chimera` already applies.
  - `stage3.c`: return NULL from `Stage3_merge_local` when either side is NULL.
- **Verification:**
  - The crashing contig now aligns: NC_026503.1:1333341-1333533 (-), 94% identity, Target 102-290. This is consistent with its sibling isoform.
  - 2,000 other contigs are byte-identical between the unpatched and patched builds (md5 052c1fe2…), with the modified function entered 34 times.
  - Minimal reproducer: a 30-kb genome window plus the first 200 bp of the query.
- **Relation to existing work:** the user's repo `~/projects/gmap-gsnap` (hyphaltip/gmap-gsnap, branch `fix/stage3-cigar-segv`, commit `c02c02d`, PR #4) already has the NULL guards in `Stage3_merge_local` and `Stage3_merge_chimera`, but no `Pairpool_clean_join` change; its notes say the cause was unknown. This work supplies the cause and the `pairpool.c` fix.
- **Not tested:** output formats other than `-f 3`, PMAP mode, and the patch on the 2021/2023 trees. A pre-existing small leak (`pairarray` on the NULL return) is not addressed.
- **Needs a user decision:** add the `pairpool.c` hunk to hyphaltip/gmap-gsnap PR #4, and whether and when to send the upstream report.

### D57: D54 completed: A. nidulans cds_blind (own evidence) = 54.7/55.1, same as cds_strand
- **Owner / status:** SELECT / done. Completes the D54 table; the conclusion (keep tx_strand) is unchanged.
- A. nidulans, own evidence (predicted genes, locus Sn/Pr, intron-chain Pr): tx 4,911, 54.7/55.7, 54.2; cds_strand 4,975, 54.7/55.1, 54.0; cds_blind 4,975, 54.7/55.1, 54.0.

### D58: D45 accepted by the user: the single-exon share for R6 option (b) comes from GeneMark-ES
- **Owner / status:** user decision 2026-09-26 (asked directly by REVIEW). Implementation is SELECT's (already in place per D45).
- **Decision:** GeneMark-ES models (`predict_misc/genemark.evm.gff3`) estimate the single-exon share. Protein alignments are the fallback when GeneMark output is missing. This replaces D38's "BUSCO models" wording.
- **Evidence:** mean absolute error against RefSeq single-CDS share: GeneMark-ES 3.8 points (slightly high); protein alignments 8.8 points (too low). BUSCO models are not available on the PASA-training path.

### D59: Production rerun plan and container sequence (user decisions 2026-09-26)
- **Owner / status:** user decisions (grilling by REVIEW).
- **Rerun timing (Q2 = b):** merge the fixes first, rebuild the images, run a stratified pilot of about 50 genomes, then the full rerun. No rerun with the rc.1 image.
- **Rerun scope (Q3):** all three groups, in two waves:
  - Wave 1: the 5,599 F1-affected genomes plus the 723 non-F1 candidates. The identity gate chooses PASA, relaxed PASA or BUSCO for each.
  - Wave 2: the 3,052 pretrained dependents, re-predicted from the new wave-1 parameters.
  - The pilot includes 09-11 mild-signature genomes, to confirm their cause.
- **User condition:** the fixes must be in the containers before any BFD rerun. The user asked whether they can tag/push PASApipeline, update the pin in `install_scripts/pixi_install_pasa.sh` (now `PASA_COMMIT=v2.6.1-rc.1`), and rebuild `funannotate-base` (`container-base.yml`, workflow_dispatch) and the app image (`container.yml`) through GitHub CI.
- **PASA tag timing (Q4 = b):** wait for the R1+R2 criterion-3 runs on Botrytis and H99 (job 29107674). Then merge `fix/unspliced-orient-clustering` into `rust_optimize` and tag (proposed `v2.6.1-rc.2`).
- **Still open for the user:**
  - who runs the merge, tag and push;
  - whether the `pasa_rust` removal goes into this tag;
  - removing the stale PASA worktrees (`PASA_wt_base`, `PASA_wt_binlookup`, `PASA_wt_sharding`, all merged).
- **Full chain:** PASA merge and tag → user approves the R3/R5 commit (SELECT) → merge R1 PR #1210 → update the PASA pin → rebuild the base, then the app image → pilot → wave 1 → wave 2.

### D60: Correction: no consistent BUSCO-vs-PASA training winner on same-strain genomes
- **Owner / status:** SELECT / done (measurement). Corrects an interim statement sent to REVIEW ("PASA wins on same-strain Botrytis").
- Holdout locus Sn/Pr, fixed evidence: A. nidulans busco-forced 54.8/56.8 vs tx 54.6/55.8; Botrytis busco-forced 82.0/82.0 vs tx 83.3/82.2. The direction differs by genome and the magnitude is about 1 point.
- This matters for D37's rule "<90% identity → BUSCO": that cut point is not challenged, but these data do not show that PASA training beats BUSCO on same-strain reads either. The divergent-read case is re-tested in D53.

### D61: R2/F4 criterion 3 passes; R2/F4 meet the D39 acceptance rule
- **Owner / status:** REVIEW / measured (job 29107674; nodes r27/r35).
- **Criterion 3:** PASA with `--ALT_SPLICE` completed with exit 0, and no splice_variation CHECK errors, for R1+R2 on Botrytis, H99 and S. commune. It had already passed on N. crassa and A. nidulans (D43).
- **Botrytis R1 → R1+R2:** exact chains 5,644 → 5,649; precision 65.3 → 65.6%; redundant models 394 → 384; single-CDS 28.1 → 27.9%. Neutral to slightly positive.
- **H99:** no R1-only reference run (SELECT job 29107600 left no `r1_fix/Cryptococcus_neoformans_H99.rust`), so only criterion 3 is verified there.
- **Conclusion:** R2/F4 pass criteria 1-3 on every genome where they could be compared. Criterion 4 (testing the flags separately) is not needed, since they pass together.
  - funannotate may enable them by default. SELECT's passthrough is `--pasa_unspliced_join_spliced` / `--pasa_one_alignment_per_cdna`.
  - Their value on top of R5 comes from SELECT's rank benchmark.
- **This unblocks the PASA merge and tag** (D59, Q4 = b).
- **Correction to REVIEW's earlier summary to the user:** "PASA beat BUSCO on Botrytis" must not be carried forward. Per SELECT, same-strain genomes show no consistent winner (A. nidulans: BUSCO about +1 precision point; Botrytis: PASA about +1).

### D62: PASApipeline v2.6.1-rc.2 tagged and pushed; stale worktrees removed (user-approved 2026-09-26)
- **Owner / status:** REVIEW / done.
- **Merge:** `fix/unspliced-orient-clustering` merged into `rust_optimize` (merge commit `01c4417`). Perl unit tests pass on the merged tree.
- **Version:** `Launch_PASA_pipeline.pl` $VERSION changed to `2.6.1-rc.2+rust` (commit `1044957`).
- **Tag:** annotated tag **`v2.6.1-rc.2` → `1044957`**, pushed to hyphaltip/PASApipeline. It contains the F1 fix (`4376a22`), R2/F4 (opt-in), the gmap crash fallback (default) and the changelog.
  - `pasa_rust` is **not** removed, per the user. Removal needs a coordinated change, because funannotate's `pixi_install_pasa.sh` checks for `bin/pasa_rust`.
- **After the tag:** commit `8a4e3de` on `rust_optimize` adds review analysis scripts and a `.gitignore` entry for the seqclean build outputs. It contains no pipeline code. The working tree is clean.
- **Worktrees removed:** `PASA_wt_base`, `PASA_wt_binlookup`, `PASA_wt_sharding` (their branches were already merged).
  - Their untracked scripts were backed up first to `/bigdata/stajichlab/jstajich/projects/funannotate/removed_worktree_files_20260926/` (`test_ab_repos.sh`, `test_tmpdir_effect.sh`). Only `bin/` build output was discarded.
  - Local branches were not deleted.
- **Next in the chain (D59):** the user approves the R3/R5 commit (SELECT) → merge R1 PR #1210 → set funannotate `PASA_COMMIT` to `v2.6.1-rc.2` → rebuild the base image, then the app image via CI → pilot → rerun waves.

### D63: funannotate side of the R2/F4 PASA flags (feature-detected, opt-in); version pin stays with REVIEW
- **Owner / status:** SELECT / done in the working tree (uncommitted); part of the planned R3/R5 commit.
- **Code:**
  - `train.pasa_feature_flags` searches the installed `Launch_PASA_pipeline.pl` for each flag name and passes `--UNSPLICED_JOIN_SPLICED` / `--ONE_ALIGNMENT_PER_CDNA` only if present; otherwise it warns. PASApipeline v2.6.1-rc.2 (D62) has them.
  - New train options `--pasa_unspliced_join_spliced` and `--pasa_one_alignment_per_cdna` (off by default); `runPASAtrain(extra_flags=...)`.
  - Nextflow: params `train_pasa_unspliced_join_spliced` and `train_pasa_one_alignment_per_cdna` (default false; passed only when true). `nextflow lint`: no errors (the 5 warnings predate these changes).
- **Tests:** `PasaFeatureFlagTests` (4). The full funannotate suite passes: gates 25, selection 31.
- **Pin:** REVIEW does `install_scripts/pixi_install_pasa.sh` PASA_COMMIT v2.6.1-rc.1 → v2.6.1-rc.2 as a separate commit after R1 merges (user's D59 order). SELECT does not touch it.

### D64: R6 option (b) result: single-exon sensitivity +5 to +8 points; multi-exon Sn −0.3 to −0.9, Pr +1.2 to +2.4 (**needs user decision** under D38's rule)
- **Owner / status:** SELECT / measured on 3 genomes; low-keeper se arms pending. D38 said to make (b) the default only if it helps without a drop in multi-exon accuracy. Multi-exon Sn drops slightly, so the user decides.
- **Method:** `predict_arms/single_exon_score.py` → `predict_arms/single_exon_scores.tsv`. Holdout chromosomes, exact CDS-chain match to RefSeq protein-coding mRNAs; Sn per RefSeq gene, Pr per predicted model. se = `--training_single_exon`; tx2 = same code_new2 snapshot without it; step B uses fixed evidence.

  | Genome | Single Sn | Single Pr | Multi Sn | Multi Pr | Net exact genes |
  |---|---|---|---|---|---|
  | A. nidulans | 71.3 → 78.5 | 58.6 → 55.6 | 47.7 → 47.4 | 50.6 → 51.8 | +35 |
  | Botrytis | 64.2 → 71.8 | 74.0 → 73.2 | 79.5 → 79.2 | 77.2 → 79.0 | +81 |
  | N. crassa | 49.4 → 54.7 | 65.9 → 60.6 | 59.9 → 59.0 | 65.4 → 67.8 | +19 |

- **Training effect:** A. nidulans admitted only 38 supported single-exon models (the cap, 289, did not bind). A small number of single-exon training genes changes Augustus' single-exon calls a lot.
- **The busco-forced arm shows the same single-exon pattern** (A. nidulans single Sn 76.9, Botrytis 67.0, N. crassa 57.0). BUSCO training sets contain single-exon genes, which supports the mechanism behind R6: with zero single-exon training genes, Augustus under-calls them.
- **Gene-level holdout scores (gffcompare):** A. nidulans locus Sn/Pr 54.7/55.9 (tx2) → 55.7/56.7 (se). BUSCO C 98.5 → 98.3.

### D65: Relaxed PASA validation (--pasa_num_bp_splice 0 --pasa_min_avg_per_id 90): +26% exact chains on divergent reads, no harm on clean reads
- **Owner / status:** REVIEW / measured. Jobs 29107756/57 on r-nodes, exit 0. This sets the D37 divergent-category values, pending a user decision and the predict arms.
- **Spliced custom alignments valid:** N. crassa 64% → 88% (blat 34% → 40%); A. nidulans 95% → 97%.
- **Training models against RefSeq:**

  | Arm | Exact chains | Precision | RefSeq genes hit | Redundant | Single-CDS |
  |---|---|---|---|---|---|
  | N. crassa fixed mm2 + blat | 2,096 | 26.1% | 5,941 | 1,831 | 46.3% |
  | N. crassa fixed mm2 + blat, relaxed | **2,650** | 25.4% | **7,154** | 2,921 | 47.1% |
  | N. crassa fixed mm2 only, relaxed | 2,610 | 25.0% | 7,137 | 2,928 | 47.4% |
  | N. crassa fixed mm2 + gmap, id 90 only (D51) | 2,301 | 23.5% | 6,717 | 2,682 | 51.6% |
  | A. nidulans fixed mm2 + blat | 3,220 | 46.3% | 6,282 | 325 | 23.9% |
  | A. nidulans fixed mm2 + blat, relaxed | **3,260** | **46.6%** | 6,327 | 329 | 23.8% |
- **Interpretation:**
  - On divergent reads, dropping the 3-bp junction rule gives most of the gain: +349 exact chains beyond the identity change alone.
  - The added redundancy is what R5's one-model-per-locus removes.
  - On clean reads, relaxed validation is slightly positive, not harmful. So relaxed validation could become the default for all genomes, not only divergent ones. Botrytis, H99 and S. commune have not been tested.
- **Next:**
  - SELECT predict arms on `xs/mm2fix_blat_relaxed/Neurospora_crassa_OR74A.rust` (relaxed PASA against busco-forced at 66.0/74.1).
  - Optionally, relaxed runs on Botrytis, H99 and S. commune before any "relaxed for all" default.

### D66: Predict arms txRel (relaxed PASA validation, REVIEW D65) added to the fixed-PASA series (D53)
- **Owner / status:** SELECT / running.
- **Inputs:** `xs/mm2fix_blat_relaxed/<genome>.rust` (R1 v1 library, fixed mm2 + blat, `--pasa_num_bp_splice 0 --pasa_min_avg_per_id 90`) → getBestModel tx_strand with code_new (the same code as tx/txR1, for a like-for-like comparison) → `predict_arms/pasa_inputs/REL/` (job 29107874).
- **Arms:** txRel on N. crassa (main candidate for D37's divergent category) and A. nidulans (clean control). Each has step A, B fixed and B own. Job ids in `predict_arms/jobs_pasafix.tsv`.
- **Comparators:** N. crassa busco-forced 66.0/74.1 (locus Sn/Pr, fixed evidence), txR1 and tx; A. nidulans tx 54.6/55.8.

### D67: Divergent reads (N. crassa): fixed PASA improves the EVM evidence more than the training; BUSCO's training edge comes from single-exon genes
- **Owner / status:** SELECT / measured; two follow-up arms running.
- **Holdout locus Sn/Pr (gffcompare CDS level) and exact-match single/multi-exon Sn/Pr (`single_exon_score.py`):**

  | Arm | Evidence | Training models | Locus Sn/Pr | Single Sn/Pr | Multi Sn/Pr |
  |---|---|---|---|---|---|
  | tx (old parser) | fixed | 465 | 64.3/72.6 | 50.3/66.7 | 59.9/65.6 |
  | txR1R2 | fixed | 704 | 64.6/73.0 | 49.5/66.0 | 60.6/66.2 |
  | txID90 | fixed | 739 | 64.5/72.6 | 48.9/65.7 | 60.6/65.9 |
  | busco-forced | fixed | – | 66.0/74.1 | 57.0/61.3 | 59.6/68.5 |
  | tx | own | 465 | 63.8/71.8 | | |
  | txR1 | own | 676 | 65.4/73.1 | | |
  | txR1R2 | own | 704 | 65.8/73.6 | | |

- **Reading:**
  1. With fixed evidence (training only), the PASA fixes add just +0.3/+0.4. With their own evidence, txR1R2 gains +2.0/+1.8 over tx. So the PASA fixes mostly improve the EVM evidence.
  2. busco-forced keeps a training-only lead of about 1.4/1.1. Its edge is single-exon sensitivity (57.0 vs about 49-50) and multi-exon precision. PASA training has 0% single-exon genes (D23).
- **Follow-up arms** (`predict_arms/jobs_pasafix.tsv`):
  - **seR1R2:** R1R2 PASA with `--training_single_exon` (code_new2). Does closing the single-exon gap make PASA training match BUSCO?
  - **buscoR1R2:** BUSCO-forced training with R1R2 PASA as EVM evidence (step B own). This is what the predict gate's BUSCO fallback does in production, and it may combine both strengths.
- These decide D37's divergent-category action. Pending also: txRel (relaxed validation) and txG (guarded ranking).

### D68: Guarded ranking (D55) is not confirmed by predictions on 2 of 3 genomes: mixed, under 1 point
- **Owner / status:** SELECT / interim; the A. nidulans txG arm is pending. Decision on D55 deferred until then.
- **Holdout locus Sn/Pr, tx (complete-first) vs txG (guarded), rc1 PASA:**
  - Botrytis: fixed 83.3/82.2 vs 83.3/82.1; own 82.9/81.4 vs **83.3/82.0**.
  - N. crassa: fixed **64.3/72.6** vs 63.9/71.9; own 63.8/71.8 vs 63.6/71.6.
  - Training models: 1,284 → 1,251 and 465 → 441 (guarded selects fewer complete models, as designed).
- **Reading:** the evidence-level gain (D55: +117 to +235 exact intron chains) does not carry through to the final predictions consistently. The direction differs by genome and the magnitude is under 1 point.
- **Options once A. nidulans is in:**
  - (1) keep guarded (evidence-level better, prediction-neutral);
  - (2) revert to complete-first (what every other predict arm used);
  - (3) make the 0.8 ratio a parameter.
  No data yet favors (1) or (2).

### D69: getBestModel default ranking reverted to complete-first; guarded kept as a parameter (supersedes D55's default)
- **Owner / status:** SELECT / done in the working tree; **proposed for user review** (the user may prefer guarded).
- **Evidence (holdout locus Sn/Pr, tx complete-first vs txG guarded):**

  | Genome | Fixed evidence (training only) | Own evidence |
  |---|---|---|
  | A. nidulans | 54.6/55.8 vs 54.6/55.8 | pending |
  | Botrytis | 83.3/82.2 vs 83.3/82.1 | 82.9/81.4 vs 83.3/82.0 |
  | N. crassa | 64.3/72.6 vs 63.9/71.9 | 63.8/71.8 vs 63.6/71.6 |

  With training only, guarded is never better (tie, tie, −0.4/−0.7) and gives fewer training models (1,122→1,111; 1,284→1,251; 465→441). With own evidence it is mixed and under 1 point.
- **Why revert:** no measured prediction benefit. Also, every other predict arm (the locus rule D54/D57, R6 D64, the PASA-fix series D67) used complete-first, so the committed default matches what was measured.
- **Code:** `train.pick_locus_model(..., complete_min_frac=0.0)`. 0 = complete-first (default); 0.8 = guarded. Tests: `GuardedRankingTests` (default complete-first, plus two guarded cases). Full selection/gate suites pass (32/25).

### D70: D69 completed: A. nidulans own evidence, tx 54.7/55.7 vs txG 54.7/55.8 (neutral)
- **Owner / status:** SELECT / done. Completes the D69 table; the conclusion (complete-first default, guarded opt-in) is unchanged.

### D71: Interim: fixed PASA (R1+R2) plus single-exon training (seR1R2) is the best N. crassa arm; relaxed validation does not beat R1+R2
- **Owner / status:** SELECT / interim. Pending: buscoR1R2 (the like-for-like comparator: BUSCO training + the same R1R2 evidence) and seR1R2 with fixed evidence.
- **N. crassa holdout locus Sn/Pr (intron-chain Sn/Pr):**
  - seR1R2 own **67.1/75.0** (64.3/76.0);
  - txR1R2 own 65.8/73.6 (65.0/73.1);
  - txRel own 65.8/72.9 (64.8/73.6);
  - busco-forced fixed 66.0/74.1 (62.6/75.6).
- **Relaxed validation (txRel, REVIEW D65):** training models 861 (vs 704 for R1R2), but predictions are not better: fixed 64.5/72.6 vs 64.6/73.0; own 65.8/72.9 vs 65.8/73.6. A. nidulans control: txRel fixed 54.5/55.7 vs tx 54.6/55.8 (no harm).
- **Tentative reading for D37's divergent category:** fixed PASA (R1 + R2) with single-exon training, not relaxed validation and not BUSCO training. To be confirmed by buscoR1R2 and seR1R2 fixed.

### D72: REVIEW withdraws the relaxed-validation proposal (D65) after the predict results
- **Owner / status:** REVIEW / withdrawn. Supersedes the D65 "next" items.
- **Evidence (SELECT predict arms, N. crassa holdout, locus Sn/Pr, own evidence):**

  | Arm | Locus Sn / Pr |
  |---|---|
  | txRel (relaxed validation) | 65.8 / 72.9 |
  | txR1R2 | 65.8 / 73.6 |
  | seR1R2 (R1 + R2 + single-exon training) | 67.1 / 75.0 |
  | busco-forced (fixed evidence) | 66.0 / 74.1 |

  - Fixed evidence: txRel 64.5/72.6 against R1R2 64.6/73.0.
  - A. nidulans control: txRel 54.5/55.7 against 54.6/55.8.
- **Conclusion:** relaxed validation adds training models (861 vs 704) and exact chains at the evidence level (+26%, D65), but does not improve predictions, and it costs about 0.7 precision points.
  - Keep the standard PASA validation (`--pasa_num_bp_splice 3`, `--pasa_min_avg_per_id 95`) for all categories.
  - Do not run the "relaxed for all" test.
  - The D37 "adapt PASA" action for divergent reads becomes, tentatively: R1 + R2 + R6 (b) single-exon training. The final call waits for buscoR1R2 (the like-for-like comparator) and seR1R2 fixed.
- **Lesson:** evidence-level gains (exact chains in training models) do not reliably carry through to prediction accuracy. Keep scoring at the prediction level before choosing defaults.

### D73: Divergent-read result (N. crassa): fixed PASA + single-exon training ≥ BUSCO training > fixed PASA alone
- **Owner / status:** SELECT / done (1 genome; margin under 1 point). Supersedes the interim D71. **Informs a user decision** on D37's divergent category and on the R6 (b) default (D64).
- **Holdout locus Sn/Pr:**

  | Arm | Fixed evidence (training only) | Own evidence |
  |---|---|---|
  | seR1R2 (R1+R2 PASA, `--training_single_exon`) | **66.4/74.9** | **67.1/75.0** |
  | buscoR1R2 (BUSCO-forced training, R1R2 PASA evidence) | 66.1/74.2 | 66.8/74.4 |
  | busco (BUSCO-forced, rc1 PASA) | 66.0/74.1 | – |
  | txR1R2 (R1+R2 PASA, no single-exon training) | 64.6/73.0 | 65.8/73.6 |
  | se (rc1 old-parser PASA + single-exon) vs tx2 | 65.1/73.6 vs 64.0/72.3 | – |

  Own evidence, exon Pr / intron-chain Pr: seR1R2 83.5/76.0; buscoR1R2 82.7/75.2; txR1R2 82.0/73.1.
- **Conclusions:**
  1. Single-exon training is the decisive factor for PASA training (+1.1 to +1.9 locus Sn), on both old-parser and fixed PASA.
  2. With it, fixed-PASA training matches or slightly beats BUSCO training (+0.3/+0.6-0.7).
  3. Without it, BUSCO training is better (by 1.5/1.1).
  4. BUSCO training plus fixed-PASA evidence (the predict gate's fallback) is a close second: safe as a fallback.
- **Suggested action for D37's divergent category (90-99% identity):** fixed PASA (R1, R2 flags) + R6 (b) single-exon training; do not relax validation (D71: txRel not better); BUSCO training as the fallback when the gate fails. **The user decides.**

### D74: Low-keeper arms resubmitted (harness path bug) and extended with R1+R2 arms
- **Owner / status:** SELECT / running.
- **Bug:** the first low-keeper step A jobs failed immediately. `arm.sh` reads the PASA variants from `refseq_benchmark/`, but the low-keeper benchmark wrote them to `refseq_benchmark_lowkeep/`, so every step B stayed `DependencyNeverSatisfied`. The PASA runs, benchmarks and prep were fine.
- **Fix:** symlinked `refseq_benchmark_lowkeep/*.best.*.gff3` into `refseq_benchmark/`. Moved the failed arm dirs to `*.failed_missing_bm_20260926`.
- **Resubmitted (50 jobs, `predict_arms/jobs_lowkeep_v2.tsv`):** H99 (same-strain reads, F1-collapsed in production) and S. commune (divergent reads, 94% identity).
  - Arms: old, oldnew, tx, cdsS, cdsB, busco, tx2, se.
  - New: txR1R2 and seR1R2 on REVIEW's r1_r2 PASA (getBestModel job 29108020), testing the D73 best configuration on collapsed genomes.

### D75: H99: with rc.1 PASA (F1 fixed), the production collapse disappears; selection is neutral
- **Owner / status:** SELECT / done (the H99 arms that have finished).
- **Keepers:** production H99 (F1-affected PASA, D47/REVIEW) had 74 filterGeneMark keepers and 160 complete models. The rc.1 rerun has 892 keepers (826 multi-CDS) and 2,002 complete of 2,704 PASA models. So H99 is not a low-keeper genome once F1 is fixed. This directly confirms REVIEW's diagnosis that F1 caused the collapse, and supports the D47 rerun list.
- **Holdout locus Sn/Pr:**
  - fixed evidence: old 80.1/80.9, oldnew 80.1/81.0, tx 79.9/80.8, se 80.1/80.9, cdsS 79.9/80.9;
  - own evidence: tx 80.5/81.4, cdsS 79.8/78.8.
  - Selection and R6 (b) are neutral here. The locus rule holds (cds_strand own −2.6 Pr).
- **Consequence:** a genuinely low-keeper genome (after F1 is fixed) still needs to be tested for R3's effect. S. commune (divergent reads) is pending. A candidate from the production low-keeper list should be re-checked against the F1 scan first.

### D76: Correction to D29: most low-keeper genomes are F1-affected; with F1 fixed, expect about 10% low-keeper genomes, not 24%
- **Owner / status:** SELECT / done (measurement).
- **Method:** a new random sample of 400 production PASA-trained genomes (`qc_training_sweep_20260925/stageA.tsv`). filterGeneMark keeper count from each `funannotate-predict.log`; F1 status from REVIEW's `production_f1_scan.tsv` (in the window 06-29 to 09-24 with frac_not_mod3 > 0.01).
- **Result:** F1-affected 275 genomes, 82 with fewer than 200 keepers (30%); clean 125 genomes, 13 with fewer than 200 keepers (10%).
- **Consequence:** the 24% in D29 was mostly the F1 bug (confirmed on H99, D75). After a rerun with the fixed image, R3 changes the training set in about 10% of genomes. Most of the benefit of the planned reruns comes from the F1 fix (rc.1 PASA) and R1; R3/R5 act as safeguards for the remaining low-keeper genomes.

### D77: H99 adds a 4th genome to D64 (R6 b) and to the BUSCO-vs-PASA comparison (D60)
- **Owner / status:** SELECT / done. Adds evidence for the pending user decisions D64 and D73.
- **R6 option (b), se vs tx2, fixed evidence:**
  - single Sn/Pr 50.6/34.5 → 55.8/32.8 (only 77 RefSeq single-CDS genes on the holdout, so noisy);
  - **multi Sn 74.2 → 74.2 (no drop)**, multi Pr 76.2 → 76.6;
  - locus 79.9/80.8 → 80.1/80.9.
  - Across 4 genomes: single Sn +5 to +8 in all; multi Pr up in all; multi Sn −0.9 to 0.0.
- **BUSCO vs PASA training, same strain:** H99 busco-forced 78.1/80.6 vs tx 79.9/80.8 (multi Sn 71.3 vs 74.2). PASA training is better here. The same-strain tally with F1-fixed PASA is now PASA ahead on H99 (+1.8) and Botrytis (+1.3), and BUSCO ahead on A. nidulans (+0.2 Sn / +1.0 Pr).

### D78: Interim: on a genuine low-keeper genome (S. commune), R3 without the gate makes predictions worse; the gate is essential
- **Owner / status:** SELECT / interim; busco, se, txR1R2 and seR1R2 arms pending.
- **S. commune with rc.1 PASA** (divergent reads, 94% identity): 822 PASA models, 93 complete; 30 filterGeneMark keepers. So it is low-keeper even with F1 fixed.
- **Holdout locus Sn/Pr, fixed evidence (the gate is disabled in these arms, `--min_pasa_complete_models 0`):**
  - old (partial models kept, 552 training models) **29.2/44.2**, intron chain 34.0/43.8;
  - oldnew 24.5/35.5 (87 models); tx 24.9/36.0 (92); cdsS/cdsB 25.1-25.2/35.7-35.8 (109).
- **Reading:** with about 90 complete genes, Augustus/SNAP train worse than with 552 mostly partial multi-exon fragments. R3 alone is harmful at very low counts. In production the D07 gate (93 complete < 500) would switch this genome to BUSCO training. So R3 must ship together with the gate, and the gate threshold matters. The busco arm will show whether that fallback rescues it.
- **Implication for the commit:** R3 and the predict gate must not be separated; the D07 default 500 stays on.

### D79: S. commune validates the predict gate: BUSCO fallback beats every PASA option by about 8 points
- **Owner / status:** SELECT / done. Completes D78.
- **S. commune (93 complete PASA models, 30 keepers, divergent reads), holdout, fixed evidence:**

  | Arm | Locus Sn/Pr | Intron chain Sn/Pr | Multi Sn/Pr | Proteome BUSCO C |
  |---|---|---|---|---|
  | busco-forced (= D07 gate fallback) | **37.4/49.2** | **42.3/48.6** | **37.2/42.7** | **97.5** |
  | old (production path) | 29.2/44.2 | 34.0/43.8 | 29.0/37.3 | 86.9 |
  | txR1R2 | 28.8/39.2 | 32.5/39.0 | 27.7/33.2 | 85.0 |
  | tx (gate off) | 24.9/36.0 | 28.1/36.2 | 23.7/30.5 | 79.4 |

- **Conclusions:**
  - With very few complete PASA models, BUSCO training is far better than any PASA training, including the old production path. The D07 gate (93 < 500 → BUSCO) makes the right call here.
  - R3 without the gate is harmful (D78), so R3 and the gate ship together.
  - The 500 threshold separates this genome (93) from the ones where PASA training is fine (N. crassa 2,146; A. nidulans 5,181; Botrytis about 5,000; H99 2,002 complete). No genome tested falls between 157 and 2,000, so the exact threshold is not pinned down.

### D80: User decision on D64: R6 option (b) (`--training_single_exon`) is ON by default
- **Owner / status:** user decision 2026-09-26 (asked directly by REVIEW). Implementation: SELECT sets the default to on in the R3/R5 commit.
- **Basis:**
  - N. crassa: best arm (seR1R2 own 67.1/75.0).
  - A. nidulans and Botrytis: +35 and +81 exactly correct genes; multi-exon Pr +1.2 to +1.8; multi-exon Sn −0.3.
  - H99: neutral. S. commune: no single-exon genes admitted.
  - The small multi-exon Sn drop (0.3-0.9, no confidence interval) is accepted in exchange for the precision and net correct-gene gains. No bootstrap was requested.

### D81: User decisions: divergent category implementation, and calibration of the complete-model gate (experiments A then B, conservative rule)
- **Owner / status:** user decisions 2026-09-26 (asked by REVIEW). Experiments are a joint REVIEW + SELECT effort.
- **Implementation accepted (D37/D73):**
  - No special PASA settings for divergent reads (standard validation; relaxed validation withdrawn in D72).
  - Every genome runs fixed PASA (R1 + R2) + R6 (b) (D80).
  - Two gates route to BUSCO training: read identity below 90% (D37), and complete PASA models below the threshold (D07; now 500).
  - The identity category is recorded per genome; it acts only below 90%.
- **Threshold calibration: both experiments, A first; conservative rule.** The user's reason: 500 was set on F1-corrupted output, and counts after the fix are much higher (H99: 74 → 2,002).
- **Experiment A: within-genome titration** (run first).
  - Genomes: A. nidulans, Botrytis, H99 and N. crassa. Use fixed PASA output (`r1_r2/` where available).
  - Subsample complete PASA training models to N = 50, 100, 200, 300, 500, 750, 1000 and 2000, with 5 random draws for N ≤ 500.
  - Train Augustus/SNAP from each subset with fixed evidence. Score the holdout chromosomes with `predict_scorer.py` (locus/exon/intron-chain F1) and proteome BUSCO.
  - Compare with busco-forced training on the same genome.
  - Also record the keepers after filterGeneMark per subset (complete models or keepers as the gate variable?).
- **Experiment B: cross-genome validation** (with the fixed-image rerun pilot, D59).
  - About 40 of the 254 RefSeq and PASA-trained BFD genomes, stratified by post-fix complete-model count and by identity category.
  - Each is run PASA-trained and busco-forced, with the same evidence.
  - Fit the difference in holdout locus F1 against log(complete models), with bootstrap 95% CIs over genomes. Test the effects of identity and genome size.
- **Decision rule (accepted: conservative):** the new threshold is the smallest N where the lower 95% bound of (PASA-trained − BUSCO-trained) holdout locus F1 is at least 0.
- **Roles:**
  - SELECT: the train/predict arm machinery, subsampling hooks, submission and post-run collection.
  - REVIEW: the experimental design, the scoring and CI/fit analysis (bootstrap and crossover estimate), and a check of the collected data.
  - All jobs exclude c[01-30] and use frozen code snapshots.

### D82: D80 implemented (R6 b on by default); experiment A (titration) set up and running
- **Owner / status:** SELECT / done (code, uncommitted) and running (experiment).
- **D80:** `predict --training_single_exon` is now on by default (`--no_training_single_exon` turns it off). Tests pass (selection 32, gates 25).
- **Experiment A (D81)**, directory `predict_arms/titration/`:
  - Code snapshot `predict_arms/code_new4` (diff md5 be240f23…): complete-first getBestModel, R3/R5, D41 fixes, R6 (b) on, R2/F4 passthrough.
  - Pools (`make_inputs.py`, job 29108788): getBestModel on `r1_r2/<g>.rust`, complete models on train chromosomes: N. crassa 1,999; A. nidulans 2,872; Botrytis 3,574; H99 2,210.
  - Tasks (`tasks.tsv`): N = 50/100/200/300/500 with 5 draws, 750/1000/2000 with 3 draws (REVIEW's request), skipping N above the pool (N. crassa: no 2000). 133 tasks. Seeds from `titration_lib.seed_for(genome, N, draw)`.
  - Each task (`run_task.sh`, one job, --exclude=c[01-30]):
    - step A trains on the train chromosomes from the subset (gate off: `--min_pasa_complete_models 0`, R6 (b) on);
    - step B predicts the full genome with the same fixed evidence as the fixed arms (PASA = rc1 old getBestModel), image-native code;
    - predict_scorer.py on the holdout chromosomes + proteome BUSCO.
    - Row → `rows/<g>_N<N>_d<draw>.tsv` with genome, N, draw, n_subset, n_complete_used (models passing selection), n_keepers (filterGeneMark), n_single_admitted, locus/exon/intron_chain Sn/Pr, pred_genes, busco_c, exit codes, runtime_s, node.
  - Jobs: smoke test (task 0) then the main array 1-132 (%40) with afterok (`titration/jobs.tsv`).
  - Comparator rows (N="busco") come from the existing busco.B.fixed arms, with the same evidence.
- **Aggregation:** `rows/*.tsv` → `titration_scores.tsv` (REVIEW's column names) when done. REVIEW does the analysis (`PASApipeline/code_review_20260925/titration_analysis.py`).

### D83: Training-decision audit log in funannotate; methods write-up for the paper
- **Owner / status:** SELECT / done in the working tree (uncommitted; part of the pending commit). Requested by the user: the log files must show what training data is used, where it comes from, which thresholds are crossed, and why. A paragraph with empirical tables must document the testing for the paper.
- **Code:**
  - `library.set_training_decision_log`, `record_training_decision` and `training_decision_summary` write `logfiles/training_decisions.tsv` (command, stage, decision, value, threshold, outcome, reason, timestamp), echo `TRAINING-DECISION ...` lines, and log a summary table before Augustus training.
  - Recorded stages: rnaseq_gate, pasa_options (train); training_mode_initial, pasa_gate / pasa_gate_override, training_mode_switch, busco_training, single_exon_share, select_complete_orf, select_keeper_filter, select_multi_cds, select_single_exon_support, select_single_exon_admit, select_redundancy, select_overlap, select_final, min_training_models, final_training_source (predict).
  - Tests: `TrainingDecisionLogTests` (5) and `SelectionDecisionLogTests` (1). Full suite passes (gates 30, selection 33).
- **Real-run check** (`predict_arms/audit_check/`, job 29108801, code frozen from the working tree):
  - A. nidulans records 2,846 complete ≥ 500 → PASA; 847 partial removed; 1,408 keepers → keeper filter on; 38 single-exon admitted (cap 289); 1,160 final.
  - S. commune records 93 < 500 → switched to BUSCO; 647 BUSCO models.
  - One misleading outcome was fixed afterwards: select_multi_cds now says "excluded unless protein-supported (R6)" when R6 is active.
- **CHANGELOG.md** updated (Added / Changed / Fixed).
- **Methods write-up:** `training_data_selection_methods.md` (copied to funannotate-live `docs/`). Narrative from `methods_narrative.md`; Tables M1-M8 generated by `methods_tables.py` from the result files (no hand-typed numbers). Experiment A results will be added as a new table when done.
- **Note for Nextflow:** train and predict both write `training_decisions.tsv`. FUNANNOTATE_TRAIN copies train's logfiles into genome_annotation/<out>/logfiles and predict then appends to it. A train rerun after predict would overwrite the predict rows (edge case; not changed).

### D84: D83 Nextflow note resolved: train appends to training_decisions.tsv
- **Owner / status:** SELECT / done (uncommitted). FUNANNOTATE_TRAIN now appends train's `training_decisions.tsv` rows (without a header) when the target file exists, instead of overwriting predict's rows. Other log files are still copied with cp -f. `nextflow lint`: no errors.

### D85: SELECT changes committed to funannotate-live (66c75ea, local, not pushed)
- **Owner / status:** SELECT / done. The user approved in the SELECT session ("commit all those changes in funannotate, add to the changes file what was added and info as well in docs folder").
- **Commit:** `66c75ea` on `target_1.9/rust_EVM_trinity_PASA`; parent `03c5f46` (PASA pin to v2.6.1-rc.2, already on the branch). 11 files, +2,127/−72:
  - funannotate/library.py, predict.py, train.py;
  - CHANGELOG.md, CHANGES.md (new branch section);
  - docs/training.rst (new user guide, linked from index.rst and predict.rst);
  - docs/training_data_selection_methods.md;
  - tests/test_training_gates.py, tests/test_training_selection.py.
- **Not included:** unrelated untracked files in the repo, and the Fungi_BFD Nextflow changes (still uncommitted; they need the user's instruction).
- **Next per D59:** REVIEW rebases R1 (fix/bam2gff3-cigar) onto 66c75ea → PR #1210 → image rebuild. Push only on the user's instruction.

### D86: Nextflow changes committed to Fungi_BFD (c5ce230, main, local, not pushed)
- **Owner / status:** SELECT / done. The user approved in the SELECT session. Files: FUNANNOTATE_TRAIN/main.nf, FUNANNOTATE_PREDICT/main.nf, conf/profile_funannotate.config (+109/−11). Excluded by the user: k8/*, rescue_earlgrey_partials.sh, conf/test_swissprot_production.config.

### D87: funannotate branch published and base image build launched (user-approved 2026-09-26)
- **Owner / status:** REVIEW / done, with the base build in progress.
- **R1 rebased** onto SELECT's R3/R5 commit `66c75ea` (parent `03c5f46`, the PASA pin to `v2.6.1-rc.2`, committed by the user).
  - The only conflict was `CHANGELOG.md`; both entries were kept.
  - New R1 commit: `8422006`. All 10 committed test files pass on the combined tree.
- **Pushed to nextgenusfs/funannotate:**
  - `fix/bam2gff3-cigar` force-pushed (`c134412` → `8422006`, with lease);
  - `target_1.9/rust_EVM_trinity_PASA` fast-forwarded and pushed (`41a2fd7` → `8422006`).
  - This publishes the pin, R3/R5 and R1 in the agreed order. **PR #1210 was merged** (GitHub, 2026-09-26 16:18 UTC).
- **Base image:** `container-base.yml` dispatched manually on `target_1.9/rust_EVM_trinity_PASA` with `tag=v1.9.0-rc.2`. GitHub Actions run 36254997032, commit `8422006`.
- **Next (user-approved):** after the base run succeeds, tag `target_1.9/rust_EVM_trinity_PASA` (`8422006`) as `v1.9.0-rc.2` and push it. `container.yml` then builds the app image FROM `funannotate-base:latest` (and `latest-norust`).
  - The tag is pushed only after the base succeeds, because a `v*` tag push triggers the app build immediately.
- **Not in this release:**
  - the Nextflow repo changes (SELECT, uncommitted);
  - the gate threshold (still 500, pending experiment A under the conservative rule);
  - the read-identity category logging (D81), unless SELECT's commit already has it.

### D88: funannotate v1.9.0-rc.2 tagged; app image build triggered
- **Owner / status:** REVIEW / done, with the app build in progress.
- **Base image:** Actions run 36254997032 succeeded (build rust/norust, sign rust/norust). Published as `funannotate-base:latest` and `:v1.9.0-rc.2`, from commit `8422006`, with PASA `v2.6.1-rc.2`.
- **Extra commit:** `d4b2655` adds `tests/test_check_rust_engines.py` (the user's request; 9 tests pass), pushed to `target_1.9/rust_EVM_trinity_PASA`.
- **Tag:** annotated `v1.9.0-rc.2` → `d4b2655`, pushed to nextgenusfs/funannotate. This triggered "Build and Sign Container Image" (run 36256284819), which builds FROM `funannotate-base:latest` and `latest-norust`.
- **Next:**
  1. The app build succeeds.
  2. Pull the `v1.9.0-rc.2` .sif into `/bigdata/stajichlab/shared/singularity_cache`.
  3. Verify in the image: PASA version `2.6.1-rc.2+rust`, the R1 `bam2gff3` present, R6 (b) on by default, and the R2/F4 flags detected.
  4. Rerun pilot (about 50 genomes, stratified, including 09-11 mild-F1 genomes), plus experiment B.
  5. Wave 1, then wave 2 (D59).

### D89: User approves removing predict_arms run data once experiment A is complete and analyzed
- **Owner / status:** SELECT / accepted by the user (SELECT session): "once you have completed the run and you know the stats from the run, we should remove the predict_arms run data".
- **Size:** predict_arms/ = 238 GB (per-arm predict outputs and intermediates, per-arm Augustus config copies, genome/BAM prep copies).
- **Plan:**
  1. Wait for all 133 titration tasks.
  2. Aggregate titration_scores.tsv (plus N="busco" rows) and hand it to REVIEW for the analysis.
  3. Update the methods doc/page and the Fungi_BFD/docs copy.
  4. Then clean up with a manifest: KEEP per-arm score.tsv, busco_short_summary.txt, run_status.tsv, arm_info.txt, logfiles/ (incl. training_decisions.tsv), parameters.json + ab_initio_parameters, and the final predict_results GFF3/proteins (compressed); REMOVE genome/BAM copies, augustus config copies, predict_misc intermediates and scratch. Record the sizes before and after.

### D90: Cleanup keep-list widened: keep every log file and the evidence inputs used (user request)
- **Owner / status:** SELECT / accepted by the user: "keeping logfiles may be important also when doing the cleanup so we know how evidence was used". Amends the D89 plan.
- **Per arm/task directory, KEEP (compressing text files over 1 MB with gzip):**
  - all of `out/logfiles/` (funannotate-predict.log, funannotate-EVM.log, augustus/busco/p2g logs, training_decisions.tsv, predict_training_gate.tsv);
  - predict.capture.log / A.capture.log / B.capture.log, busco.log, arm_info.txt, run_status.tsv, score.tsv, busco_short_summary.txt;
  - the PASA GFF3 actually given to the run (`pasa.gff3` / `pasa.subset.gff3`);
  - `predict_misc/weights.evm.txt`, `final_training_models.gff3`, `ab_initio_parameters/`, and `predict_results/*.parameters.json`;
  - final `predict_results` GFF3 and proteins.
- **Per genome prep dir, KEEP:** split.*_chroms.txt, genemark.*.gtf and the gmes log (compressed).
- **REMOVE:** genome FASTA/BAM copies (regenerable from input_clean_genomes and rnaseq_data), per-arm augustus config copies, the remaining predict_misc intermediates, and the frozen code snapshots (their diffs are kept in Fungi_BFD/docs).
- A manifest of kept and removed paths, with sizes before and after, goes into Fungi_BFD/docs/training_selection_2026-09/.

### D91: v1.9.0-rc.2 images pulled and verified; the version label is wrong
- **Owner / status:** REVIEW / done (job 29108912, node r41). The version label needs a user decision.
- **App build:** run 36256284819 succeeded (build, sign, provenance and verify, rust and norust).
- **Images** (new files; nothing overwritten):
  - `/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.2.sif`
  - `/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.2-norust.sif`
  - The two files differ: the norust image has the `*_rust` binaries renamed to `.disabled`. Their identical size is a coincidence of the image format.
- **Verified inside the rc.2 image:**
  - PASA `$VERSION = "2.6.1-rc.2+rust"`.
  - R2/F4 flags present in Launch_PASA_pipeline.pl; the gmap fallback is present.
  - F1 fix: one `print_GFF_row` call.
  - R1 `parse_minimap2_splice_record` and the `bam2gff3(strand=)` parameter are present.
  - `--no_training_single_exon` exists, so R6 (b) is on by default.
  - `train.pasa_feature_flags` detects both PASA flags.
  - The gates are present.
- **Problem:** funannotate reports `1.9.0-beta.13.dev0+gd4b2655`. `funannotate/__version__.py` hard-codes `PRERELEASE = "beta.13"`, and the label comes from that constant, not the tag. The code is correct (`d4b2655` = tag `v1.9.0-rc.2`); only the label is misleading in logs and provenance.
- **Options for the user:**
  - (a) Accept it, and record `gd4b2655` as the identifier.
  - (b) Set `PRERELEASE = "rc.3"`, tag `v1.9.0-rc.3`, rebuild the app image only (the base is unchanged).
  - (c) Move the `v1.9.0-rc.2` tag. Not recommended: the tag is already pushed.

### D92: User chose D91 option (b): v1.9.0-rc.3 fixes the version label
- **Owner / status:** REVIEW / tag pushed; app build in progress (run 36262116189).
- **Change:** commit `cd1b5ee` sets `PRERELEASE = "rc.3"` in `funannotate/__version__.py`. The code is otherwise identical to `v1.9.0-rc.2` (`d4b2655`).
  - Locally the tag resolves to `1.9.0-rc.3`, and the gate tests pass.
- **Tag:** annotated `v1.9.0-rc.3` → `cd1b5ee`, pushed to nextgenusfs/funannotate. The app image builds FROM the unchanged `funannotate-base:latest` (rc.2 base, PASA v2.6.1-rc.2).
- **Next:** pull `funannotate-1.9.0-rc.3.sif` and `-norust.sif` into the singularity cache, confirm the version reports `1.9.0-rc.3`, and use rc.3 for the rerun pilot. The rc.2 images stay in the cache as reference.

### D93: v1.9.0-rc.3 images pulled and verified: ready for the rerun pilot
- **Owner / status:** REVIEW / done. CI run 36262116189 succeeded; pull job 29111159.
- **Images:**
  - `/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.3.sif`
  - `/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.3-norust.sif`
- **Verified in the image:**
  - funannotate version `1.9.0-rc.3` (rust and norust); the label is fixed.
  - PASA `2.6.1-rc.2+rust` (the base is unchanged from rc.2; verified on rc.2 in D91). The `#…BLEEDING_EDGE` line printed by the check script is a commented-out line.
  - R2/F4 flags present; gmap fallback present; F1 fix present.
  - R1 parser and the strand parameter present.
  - R6 (b) on by default; R2/F4 feature detection works; the gates are present.
- **Next (D59/D81):** the stratified rerun pilot (about 50 genomes, including 09-11 mild-F1 genomes) and experiment B, both with rc.3. Then wave 1 and wave 2. Experiment A is still running.

### D94 — Experiment A: H99 step B failed on a dangling symlink; bind fixed and 34 tasks resubmitted (SELECT, 2026-09-26)
- Result of job 29108790: 98 of 133 array tasks completed (99 with the test task 29108789). All 34 failures (tasks 99-132) are H99. Step A passed (exit 0). Step B stopped with "…H99.best.old.gff3 is not a valid file".
- Cause: refseq_benchmark/Cryptococcus_neoformans_H99.best.old.gff3 is a symlink to refseq_benchmark_lowkeep/ (valid, 9,027,228 bytes). run_task.sh did not bind refseq_benchmark_lowkeep into the container, so the link was dangling inside apptainer. REVIEW reported that the file was never generated. That is not correct: the file exists.
- Fix: run_task.sh line 29 BINDS now includes refseq_benchmark_lowkeep. No inputs changed. The 34 tasks were resubmitted as job 29112383 (exclude c[01-30]).
- Aggregation: aggregate.py writes titration_scores.tsv (per-draw rows, status column, one N=busco row per genome from busco.B.fixed) and titration_summary.tsv (mean/SD per genome x N). Sent to REVIEW for the three complete genomes.
- Caveats: each N=busco comparator is one run (no draw variance). N. crassa has no N=2000 (pool 1,999). Step A code is the code_new4 snapshot (version label 1.9.0-beta.13.dev0+g6f3eaad), not the rc.3 image.

### D95: Experiment A interim (3 genomes): BUSCO training wins below about 300 models; the difference is under 1 point above 500
- **Owner / status:** REVIEW / interim analysis. H99 is being rerun (29112383); experiment B is pending. **No threshold change yet.**
- **Analysis:** `code_review_20260925/titration_analysis.py` on `predict_arms/titration/titration_scores.tsv` (status==ok). Output `predict_arms/titration/titration_analysis_{locus,exon,intron_chain}.tsv`.
- **Holdout locus F1, (PASA-trained − BUSCO-trained)**, mean [95% bootstrap over draws]:

  | N | A. nidulans | Botrytis | N. crassa |
  |---|---|---|---|
  | 50 | −2.7 | −3.0 | −4.6 |
  | 100 | −1.5 | −1.5 | −2.6 |
  | 200 | −1.0 | +0.1 | −1.2 |
  | 300 | −0.6 | −0.4 | −1.5 |
  | 500 | −0.4 [−0.7, −0.2] | +0.7 [+0.3, +1.2] | −0.5 [−1.0, +0.05] |
  | 750 | −0.6 | +0.2 | −1.0 |
  | 1000 | −0.3 [−0.4, −0.1] | +0.9 [+0.7, +1.2] | +0.0 [−0.3, +0.4] |
  | 2000 | +0.1 [−0.04, +0.3] | +1.7 | — (pool 1,999) |
- **Conservative N\* (D81 rule, locus):** A. nidulans none; Botrytis 1000; N. crassa none. Exon and intron-chain F1 give the same pattern (A. nidulans N* = 2000).
- **Gate variable:** complete models track F1 slightly better than keepers (Spearman 0.84-0.95 against 0.80-0.88). Keep counting complete models.
- **Interpretation:**
  - BUSCO training is clearly better below about 300 complete models (−1 to −5 points). So 500 is not too high.
  - Above 500, |difference| < 1 point except Botrytis at 2000 (+1.7).
  - Literally applied, the conservative rule sends A. nidulans and N. crassa to BUSCO at any tested N. That would be a large policy change on sub-1-point differences.
- **Caveats:**
  - The BUSCO comparator is a single run; its variance is not in the intervals, so they are too narrow.
  - 3 draws at N ≥ 750.
  - 3 species only.
- **Next:** add H99 when its tasks finish; run experiment B (about 40 RefSeq genomes, rc.3). Then decide whether the threshold moves to about 1000 or BUSCO training becomes the default, per the conservative rule with the comparator's variance included (for example, 3 repeat BUSCO runs per genome).

### D96 — Experiment A: three BUSCO-forced comparator repeats per genome with code_new4 (SELECT, 2026-09-26)
- Request (REVIEW): the N=busco comparator was one run per genome, so its run-to-run noise is missing from the PASA − BUSCO intervals.
- Seeds: funannotate has no random-seed option, and predict/library do not use Python random. Repeats therefore measure only run-to-run nondeterminism (for example thread order in BUSCO, Augustus or EVM). If the repeats are identical, the comparator variance is zero, and that result is reported as such.
- Code: the earlier busco.B.fixed arms used the code_new snapshot, and the titration step A used code_new4. To compare like with like, the repeats use code_new4 (new arm busco4 in arm.sh; REP=1..3 writes busco4_r<k>.A / .B.fixed). Backup: arm.sh.bak_20260926.
- Jobs: 12 step A + 12 step B jobs (4 genomes x 3 repeats, H99 included), listed in predict_arms/jobs_busco4_reps.tsv.
- aggregate.py: N=busco rows now come from busco4_r1-3 (draws 1-3). The earlier run is kept as N=busco_code_new, draw 1. There is a new source column.
- The same three-repeat BUSCO-forced arm is planned for experiment B on rc.3.

### D97: Caveat on the D95 interim: its comparator used a different code snapshot
- **Owner / status:** REVIEW / note.
- **Problem:** the D95 BUSCO comparator was `busco.B.fixed` built with the `code_new` snapshot, but the titration arms ran `code_new4` (SELECT, D96). So the D95 differences compare across code versions, as well as across training sources.
- **Fix:** 3 fresh `busco4_r1-3` repeats with `code_new4` for all 4 genomes. They are labelled `N="busco"`; the old run is kept as `N="busco_code_new"`.
- **Script:** `titration_analysis.py` now ignores non-numeric reference rows (commit on `rust_optimize`). It bootstraps the repeats together with the PASA draws.
- **Consequence:** D95's numbers are provisional. The analysis will be rerun when `aggregate.py` includes `busco4` and H99.
  - If funannotate is deterministic, the repeats may be identical (SELECT: there is no seed option). In that case the comparator variance is zero, and the intervals are as reported.
