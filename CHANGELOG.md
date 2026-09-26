# Changelog

## Unreleased

### Added
- **Evidence-quality gates for RNA-seq training.** RNA-seq that does not come
  from the genome being annotated (dual RNA-seq of infected host tissue, a
  mislabeled species, a stale read file) still yields a Trinity assembly and a
  PASA training set, but the PASA models are mostly partial ORFs. Augustus/SNAP
  trained on them call too few genes and outvote GeneMark in EVM. On an
  18-genome pilot, genomes whose reads mapped at 0.1-2% to their own assembly
  lost up to 60 BUSCO points (e.g. *Drepanopeziza brunnea* 97.1% -> 46.8%).
  - `funannotate train --min_rnaseq_map_rate` (default 10%): maps the first
    `--rnaseq_gate_reads` (default 200,000) left/single reads to the genome
    with minimap2 `-x splice:sr`; below the rate, train stops before
    Trinity/PASA with exit code 3 (`library.RNASEQ_GATE_EXIT`). Report:
    `logfiles/train_rnaseq_gate.tsv`.
  - `funannotate predict --min_pasa_complete_models` (default 500): counts
    complete-ORF models (ATG start, stop codon, CDS length divisible by 3) in
    `--pasa_gff`; below the minimum, Augustus/SNAP train from BUSCO instead.
    PASA models are still EVM evidence. Report:
    `logfiles/predict_training_gate.tsv`.
  - Both thresholds come from one pilot: failing genomes mapped 0.1-2% and had
    4-281 complete PASA models; the clearly good PASA case mapped 97% with
    1,258. Set either option to 0 to disable. The 500 threshold is being
    re-calibrated on PASA output without the F1 defect (see
    `docs/training_data_selection_methods.md`).
- **Training-decision audit log.** train and predict write every decision about
  which data trains Augustus/SNAP to `logfiles/training_decisions.tsv`
  (command, stage, decision, value, threshold, outcome, reason, timestamp), echo
  each one to the run log as `TRAINING-DECISION ...`, and predict prints a
  summary table before Augustus training. Recorded stages: RNA-seq gate, PASA
  options, initial training source per predictor, PASA gate, BUSCO switch,
  single-exon share, each `selectTrainingModels` step (complete ORFs, keeper
  filter, multi-CDS requirement, single-exon support/admission, redundancy,
  overlap, final count), minimum-model check and final training source.
- **Single-exon training genes (on by default; `--no_training_single_exon`).**
  Complete single-exon PASA models with ≥80% same-strand protein2genome
  coverage are admitted, capped at share/(1−share) × multi-exon models; the
  share is the single-CDS fraction of GeneMark-ES models (fallback: protein
  alignments). Before, the multi-exon requirement removed every single-exon
  gene. On 4 RefSeq genomes: single-exon Sn +5.2 to +7.6, multi-exon Pr +0.4
  to +2.4, multi-exon Sn 0.0 to −0.9 (holdout chromosomes).
- **PASA opt-in flags** `train --pasa_unspliced_join_spliced` and
  `--pasa_one_alignment_per_cdna`, passed only when the installed
  Launch_PASA_pipeline.pl supports them (PASApipeline ≥ v2.6.1-rc.2).

### Changed
- **Training-set selection (`selectTrainingModels`) uses only complete ORFs**
  (protein-based test, correct on both strands) and removes overlaps
  transitively (union-find clusters). Must be used with the PASA gate: without
  it, a genome with 93 complete models lost 4.3 points of holdout locus Sn.
- **`getBestModel` keeps one model per transcript-span locus** (same strand,
  ≥`--pasa_alignment_overlap` of the shorter model, transitive), ranked by
  complete ORF, CDS exons, CDS length, then TPM. The old one-sided,
  strand-blind, TPM-first test discarded distinct genes (about 920 complete
  loci in *Meyerozyma guilliermondii*). A guarded ranking
  (`pick_locus_model(..., complete_min_frac=0.8)`) is available; it was
  prediction-neutral on 3 RefSeq genomes, so it is not the default.

### Fixed
- `selectTrainingModels` no longer exits when no model survives filtering; it
  returns 0 so the BUSCO fallback runs (an empty FASTA made `diamond makedb`
  fail).
- The PASA gate now applies when any predictor trains from PASA, and is
  skipped on resume or with `--augustus_gff`.
- **PASA dependency (`install_scripts/pixi_install_pasa.sh`): pinned to
  `v2.6.1-rc.1` (hyphaltip/PASApipeline, `rust_optimize` branch)**, which fixes
  a duplicate-output bug traced to `PASA_transcripts_and_assemblies_to_GFF3.dbi`.
  Commit `bce776a` ("fix N+1 query in GFF3 output", 2026-06-29) added a new
  batch-query print loop to replace the old per-alignment N+1 loop, but never
  deleted the old loop — both ran unconditionally for GFF3/GTF output, so every
  validated alignment segment was printed twice. This silently doubled
  `valid_${aligner}_alignments.gff3`, `failed_${aligner}_alignments.gff3`, and
  `pasa_assemblies.gff3` (confirmed empirically: 55,205 unique lines -> 110,410
  total in a real run), degrading PASA's transcript-assembly structure and
  cascading into fewer multi-exon training models downstream (compounding the
  `getBestModel`/`selectTrainingModels` issues below). BED output was
  unaffected (both loops wrote into the same hash key, so it self-deduped),
  which is why the bug went undetected until GFF3 structure was compared
  directly against PASA 2.5.3 under matched inputs. Not an upstream PASA bug —
  traced via `git blame` to a commit on this project's own `rust_optimize`
  fork; fixed there (`PASApipeline@4376a22`), not reported upstream. Also fixed
  a checkpoint-filename typo in `Launch_PASA_pipeline.pl` that collided the
  failed-alignments `.bed` writer's checkpoint with the `.gff3` writer's,
  silently skipping the `.bed` writer on every run.
- `aux_scripts/trinity.py`: after aggregating genome-guided Trinity's per-partition
  output, `trinity.py` counted the final transcripts and logged the count but
  never checked whether it was zero before returning — so if every partition
  failed (e.g. from the `salmon_runner.pl` FASTA/FASTQ incompatibility fixed
  upstream in trinityrnaseq, or any other total-failure cause), this script
  still exited 0. `funannotate train`'s own `checkannotations()` check already
  stops the pipeline on an empty `trinity.fasta`, so this wasn't letting a
  zero-transcript run through downstream — it was this script's own exit code
  disagreeing with what actually happened, producing the self-contradictory
  log line `Trinity de novo assembly failed (trinity.py exit code: 0)`. Now
  exits 1 with a pointer to the underlying `Trinity-gg.log` when the transcript
  count is 0.
- `aux_scripts/funannotate-p2g.py`: fixed `spawn()`'s error-reporting path, which
  crashed with `AttributeError` on any exonerate failure (`p.communicate()`
  returns a `(stdout, stderr)` tuple; with stdout redirected to a file, the whole
  tuple was being assigned to `stderr` and then `.decode()`'d), masking the real
  exonerate error behind an unrelated crash.
- `aux_scripts/iprscan-local.py`: replaced fragile manual XML-header line-splicing
  (which raised `NameError` and lost all completed InterProScan chunk results
  whenever a chunk's header didn't match either hardcoded format) with the
  existing, already-correct `combine_xml()` helper that was defined but never
  called. Also fixed `os.errno.ENOENT` (`os.errno` doesn't exist in Python 3) to
  `errno.ENOENT`, which previously turned the friendly "Docker is not installed"
  message into an unrelated `AttributeError`/`NameError`.
- `aux_scripts/fasta2agp.py`: fixed a crash on any SPAdes-style contig header
  (`NODE_1_length_...`, a common assembler output format) — the fallback regex
  branch matched the same pattern twice instead of falling back to `numnamepat`,
  and `m.match(1)` (not a valid method on a `Match` object) should have been
  `m.group(1)`.
- `aux_scripts/augustus_parallel.py`: fixed float-division chunk math when
  splitting contigs >500kb for parallel Augustus prediction. Python 3's `/`
  produces floats, which both got embedded directly into Augustus's
  `--predictionStart`/`--predictionEnd` CLI args and caused the last chunk of a
  contig to fall short of the actual contig end (silently dropping up to ~10-20%
  of the contig, and any genes in it, from prediction). Now uses integer/ceiling
  division throughout.
- Added return-code/output checks (previously silently discarded) around
  subprocess calls in `aux_scripts/funannotate-runEVM.py` (EVM partition
  workers), `aux_scripts/trinity.py` (Trinity genome-guided partition workers),
  `aux_scripts/iprscan-local.py` (Docker/local InterProScan workers),
  `aux_scripts/phobius-multiproc.py` (remote/local Phobius workers, plus a
  result-count sanity check against input count), `aux_scripts/hmmer_parallel.py`
  (Pfam/dbCAN hmmsearch/hmmscan, plus a guard against `combineHmmerOutputs`
  crashing with `IndexError` if every chunk failed), `aux_scripts/augustus_parallel.py`
  (the `join_aug_pred.pl` merge step), and `aux_scripts/enrichment_parallel.py`
  (`find_enrichment.py`). All previously failed silently, leaving partial/missing
  results with no error surfaced to the user.
- `aux_scripts/funannotate-BUSCO2.py`: fixed a `NameError` (`file.close()` should
  be `nucl_file.close()`) that masked the intended "please provide a nucleotide
  file" error message; fixed a missing `s` conversion character in two
  `%(config)` format strings (`_augustus_rerun`) that raised `ValueError` on
  every Phase-2 Augustus retraining pass; fixed a wrong-index bug in
  `_get_coordinates` (`coords[busco_name][contig][2][1] = contig_end` should be
  `coords[busco_name][contig][1] = contig_end`) that corrupted the alignment
  position deque and could crash `_check_overlap` on multi-HSP tblastn hits; and
  fixed a lexical (string) version comparison for tblastn (`'2.10.0' > '2.2.31'`
  is `False`) that misclassified nearly every currently-deployed BLAST+ version
  as "not newer than 2.2.31", causing the wrong tblastn thread-count decision.
- `update.py`: fixed `getBestModels` locus grouping — the `else` branch looked up
  the previously-seen-locus map by `geneID` instead of `Loc`, so any transcript
  after the first one seen at a coordinate got grouped under a `None` locus,
  merging TPMs/transcripts from unrelated loci that happened to share a physical
  location and corrupting "keep best/alt transcript per locus" selection.
- `update.py`: fixed the exact-duplicate-coordinate gene tie-breaker in
  `GFF2tblCombinedNEW`, which checked a global `ExpressionValues` dict that was
  never populated, so duplicate-coordinate gene models were never resolved by
  expression and both survived into the final annotation. `getBestModels` now
  populates `ExpressionValues` with each locus's best TPM.
- `update.py`: fixed the same malformed minimap2 flag bug as `train.py`/`predict.py`
  (`'-ax' 'map-ont'` missing a comma, collapsing into the invalid flag
  `-axmap-ont`), which broke the long-read-only expression-mapping fallback.
- `update.py`: the minimap2 | samtools pipe used to map long reads to PASA
  transcripts now checks both processes' return codes and validates the output
  BAM before passing it to `mapCount`, instead of silently proceeding on failure.
- `annotate.py`: fixed `Gene2ProdFinal` naming — when a gene had multiple
  transcripts sharing one locus, the name written to
  `annotations.genes-products.txt` had no `_N` suffix, but `Gene2ProdFinal` was
  always populated with a suffixed name, so gene names in the curation report
  files (`must-fix.txt`, `need-curating.txt`, etc.) didn't match anything in the
  actual annotation and a `--fix` file built from them silently failed to apply.
- `annotate.py`: fixed `getEggNogHeaders`' falsy-zero bug — `if not IDi:` treated
  a valid header column index of `0` (the common case, since `query_name` is
  always the first column) the same as "not found", so the parsed v1 eggNOG
  header was always discarded in favor of hardcoded guessed column indices.
- `annotate.py`: fixed an `if/elif` that dropped InterPro term collection for any
  gene that also had a `note` annotation (EggNog/MEROPS/CAZy/BUSCO), since the two
  conditions are independently present per gene, not mutually exclusive.
- `annotate.py`: added return-code/output checks after the Pfam, dbCAN, and
  Phobius subprocess calls, which previously failed silently and just logged
  "0 annotations added" — indistinguishable from a legitimately empty result.
- `train.py`: fixed malformed minimap2 command (`'-ax' 'map-ont'` collapsed into
  the invalid flag `-axmap-ont` due to a missing comma), which broke long-read-only
  (ONT/PacBio, no short reads) training.
- `train.py`: `getBestModel` now tracks per-gene expression correctly. The dedup
  guard previously checked `transcriptID` (always unique) against a dict keyed by
  `geneID`, so it never prevented overwrites — for genes with multiple PASA
  transcripts, the "best" transcript was effectively chosen by file order rather
  than expression. Now keeps the transcript with the highest TPM per gene.
- `train.py`: fixed a header-skip typo (`targed_id` -> `target_id`) when parsing
  kallisto `abundance.tsv`.
- `train.py`: added an explicit exit when neither short nor long reads are
  available for training, instead of logging an error and continuing into
  Trinity/PASA with no input. Downgraded the short-reads-only message to `info`
  since long-reads-only training is a valid path.
- `predict.py`: fixed "HiQ" Augustus promotion, which never fired because the
  gene/mRNA `ID=` attribute prefix was never stripped before the code tried to
  reduce it to a bare gene id. Hint-supported Augustus models were silently
  missing their intended 2x EVM consensus weight.
- `predict.py`: fixed the EVM gene-separator blank-line writer, which could drop
  the blank line between two adjacent `gene` records when `previous_line` was
  unset, breaking EVM partitioning on that input.
- `predict.py`: the protein-to-genome (`funannotate-p2g.py`/Exonerate) subprocess
  call now checks its return code and output; on failure, protein evidence is
  dropped instead of silently being passed to EVM as if it succeeded.
- `predict.py`: added a missing check after the protein-mode BUSCO validation
  run so a failed run is caught immediately instead of failing further downstream
  with a less clear error.

### Performance
- `train.py`: `pOverlap` now computes interval overlap with O(1) min/max
  arithmetic instead of materializing full `set(range(...))` objects for both
  intervals, avoiding O(feature length) time/memory per call in nested loops
  over genome-wide PASA/InterLap hits.

### Changed
- `predict.py`: routed the GeneMark contig-recovery subprocess chain
  (`hmm_to_gtf.pl` / `reformat_gff.pl` / GTF-to-GFF3 conversion) through
  `lib.runSubprocess` for consistent error logging instead of unchecked
  `subprocess.call`.
- `predict.py`: collapsed the ~100-line duplicated `elif` chain that copies
  pretrained Augustus parameter files by suffix into a loop over a suffix list.
- `annotate.py`: `AllProts`/`SMgenes` are now built as sets directly instead of
  lists with O(n) membership checks that were immediately converted to sets
  anyway; also dropped a redundant unused file handle when building
  `mibig_fasta` (`SeqIO.parse` already reopens the path directly).

### Known issues (flagged, not fixed)
- `update.py:908-957`: when reusing an existing, validated MySQL PASA database,
  the alignment/import step is skipped entirely, so newly generated
  Trinity/long-read transcript evidence for this run never gets loaded into the
  database. Needs confirmation of original intent before changing — only affects
  the MySQL PASA backend, not the default SQLite path.
- `annotate.py:1712-1739` (`--renumber_antismash`): locus-tag renumbering scans
  the whole antiSMASH cluster file for any 3-9 digit run to align with cluster
  count; desyncs once a genome has ≥100 antiSMASH clusters (an extra 3-digit
  match appears from `Cluster_100`+). Pre-existing, not introduced by this pass.
- `aux_scripts/funannotate-BUSCO2.py`: ~20 call sites use `subprocess`/`p_open`
  with `shell=True` and unquoted string interpolation of `self.mainout`
  (derived from `--species`/`-o` values) and `--augustus_parameters`, a
  shell-injection-shaped pattern. Fixing this properly means auditing/rewriting
  ~20 call sites to `shell=False` with argument lists (or `shlex.quote()`), a
  larger scope better suited to its own pass. Also has an O(n²) directory-listing
  progress-check pattern (`_process_augustus_tasks`/`_process_hmmer_tasks`
  re-`os.listdir()` the shared output directory after every completed task) worth
  replacing with an atomic counter if this file gets revisited.
- `aux_scripts/funannotate-runEVM.py`: `solve_partitions()` is dead code (defined,
  never called) with raw `print()` debug statements; `RangeFinder` silently drops
  gene/evidence blocks that straddle a partition boundary with no warning logged;
  `next_start` isn't updated on one branch of the chunk-splitting loop, which
  could produce an invalid/overlapping partition in an edge case.
- `aux_scripts/funannotate-p2g.py`: `--maxintron`/`--contig_expand`/
  `--exonerate_pident` argparse args are missing `type=int` (currently latent,
  since no caller passes them explicitly, but any int arithmetic on
  `args.contig_expand` would `TypeError` if one did); an unconditional
  `tblastn -version` check runs even in the default `--filter diamond` mode,
  requiring `tblastn` on `PATH` when it's never actually used; offset parsing via
  `filename.split('__')[1]` is fragile against protein/contig IDs that themselves
  contain `__`.
- `aux_scripts/iprscan-local.py`: `basename.rstrip(".fa")` strips any trailing
  characters in the set `{.,f,a}`, not the literal suffix — a classic
  `str.rstrip()` footgun, currently harmless since it's only used for temp
  filenames.
- `aux_scripts/trinity.py`: per-partition command parsing does a naive
  `input.split(' ')` instead of `shlex.split()`, which would corrupt any
  path/argument that legitimately contained a space; all parallel workers also
  append to one shared logfile without locking, risking interleaved/garbled log
  output under concurrency.
- `aux_scripts/runIPRscan.py`: appears to be an unused legacy EBI REST client
  (the pipeline only dispatches to `iprscan-local.py`) and is broken under
  Python 3 in multiple places (bytes/str comparison mismatches in response
  polling). Candidate for deletion.
- `aux_scripts/xmlcombine.py`: `ElementTree.tostring()` without
  `encoding='unicode'` prints raw `bytes` reprs; currently dead/orphaned code
  (referenced as a path in `remote.py` but never invoked), so latent rather than
  active.
- `aux_scripts/iprscan2annotations.py`: uses `etree.iterparse` but never calls
  `elem.clear()`, defeating its memory-saving purpose on large InterProScan XML;
  neither this nor `xmlcombine.py` use `defusedxml`, a latent XXE hardening gap
  if that XML is ever externally influenced.
- `aux_scripts/funannotate-BUSCO2-py2.py` (3340 lines): unreachable dead code —
  both call sites (`predict.py:1980`, `library.py:6855`) select it only when
  `sys.version_info <= (3, 0)`, which cannot occur since the package requires
  Python ≥3.6. Candidate for deletion.

### Tests
- Added `tests/test_train_poverlap.py` covering `train.pOverlap`: full/no/partial
  overlap, asymmetric intervals, touching (non-overlapping) boundaries, and a
  large-interval case exercising the new O(1) arithmetic path.
