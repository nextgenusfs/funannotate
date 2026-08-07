# Changelog

## Unreleased

### Fixed
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

### Tests
- Added `tests/test_train_poverlap.py` covering `train.pOverlap`: full/no/partial
  overlap, asymmetric intervals, touching (non-overlapping) boundaries, and a
  large-interval case exercising the new O(1) arithmetic path.
