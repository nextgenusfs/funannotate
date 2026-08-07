# Changelog

## Unreleased

### Fixed
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

### Tests
- Added `tests/test_train_poverlap.py` covering `train.pOverlap`: full/no/partial
  overlap, asymmetric intervals, touching (non-overlapping) boundaries, and a
  large-interval case exercising the new O(1) arithmetic path.
