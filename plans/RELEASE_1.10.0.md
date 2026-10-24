# funannotate 1.10.0: planned features

Status 2026-10-24. A list to add to; nothing here is promised until it has a test or an experiment behind it.
Items marked (staged) are written but not yet released; (open) are decisions or work not started.

## 1. Python 3.10 or newer (deferred from 1.9.x; work kept on branch `release/1.10-python310`)

- Biopython 1.87+ (outside the Bio.Entrez XXE advisory range, `<=1.86`) needs Python 3.10, so clearing the advisory
  forces the move; 3.8 and 3.9 are past end of life.
- Deferred on 2026-10-24: master and `release/1.9` stay on Python 3.8/3.9 and Biopython `<1.80`; the Dependabot alert
  (GHSA-x3vf-39hj-gxr4) is dismissed as not used, since funannotate does not import `Bio.Entrez`. The changes are in
  commit 72de224 on local branch `release/1.10-python310`:  `pixi.toml`, `environment*.yml`, `conda-recipe/meta.yaml`, `setup.py` (`>=3.10,<3.12`), CI
  matrix (3.10, 3.11), `pixi.lock` solved on both platforms (Python 3.11; 3.10 also solves, 3.9 does not).
  The `Bio.pairwise2` to `PairwiseAligner` change is separate and is on master (it also runs on Biopython 1.79).
- Open before release: build the base image and run `funannotate test` (predict, update, annotate, compare) in it;
  confirm Trinity, PASA, EVM and bowtie2 source builds on Python 3.11 (not in the lock); `sklearn` `MDS` call in
  `library.py` warns on scikit-learn 1.9 and breaks on 1.10 (pin `<1.10` or port); `distutils` import in `annotate.py`
  (removed in 3.12, so the `<3.12` cap stays).
- Python 3.12 support is a separate item (`distutils`, goatools, Trinity/PASA scripts).

## 2. Split `compare` out of funannotate (open)

- Design and size estimates: `plans/COMPARE_EXTRACTION_SPIKE.md` (compare is 1,662 lines plus a 79-name, 2,266-line closure
  in `library.py`).
- Decisions needed: package name and repository; whether `funannotate compare` stays as a thin wrapper that calls the new
  package when installed; what happens to the phylogeny/orthology dependencies (proteinortho, PhyML, RAxML, trimAl,
  MAFFT) in the main image; how outgroups (`funannotate outgroups`) and the compare database files move.
- Acceptance: `funannotate test -t compare` passes identically before and after in the image.

## 3. Training-source and Augustus gates (done on master; paper work continues)

- PASA vs BUSCO decision functions, held-out Augustus accuracy report, `augustus_accuracy` row in
  `training_decisions.tsv`; `--optimize_augustus` experiment D pending; PASA+BUSCO union untested.
  See `docs/training_data_selection_methods.md`.

## 4. Translation tables and GeneMark (done on master)

- `--table` validation and warnings, GeneMark sidecar GTF handling, `funannotate check --strict --allow-missing`.
- Open: alt-start-codon `partialStart` handling; end-to-end `predict --table 12` run in the image.

## 5. Packaging and images (in progress)

- NCBI `table2asn` direct install for macOS arm64 and Linux x86_64 (`install_scripts/install_table2asn.sh`), used in
  `Dockerfile.base`; needs an image build to confirm.
- Open: published pruned images, Apptainer SIF, PyPI trusted publisher, linux-aarch64 solve.
