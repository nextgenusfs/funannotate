# funannotate 1.10.0: planned features

Status 2026-10-24. A list to add to; nothing here is promised until it has a test or an experiment behind it.
Items marked (staged) are written but not yet released; (open) are decisions or work not started.

## 1. Python 3.10 or newer (deferred from 1.9.x; work is on branch `release/1.10`, renamed from `release/1.10-python310`)

- Biopython 1.87+ (outside the Bio.Entrez XXE advisory range, `<=1.86`) needs Python 3.10, so clearing the advisory
  forces the move; 3.8 and 3.9 are past end of life.
- Deferred on 2026-10-24: master and `release/1.9` stay on Python 3.8/3.9 and Biopython `<1.80`; the Dependabot alert
  (GHSA-x3vf-39hj-gxr4) is dismissed as not used, since funannotate does not import `Bio.Entrez`. The changes are in
  commit 72de224 on that branch:  `pixi.toml`, `environment*.yml`, `conda-recipe/meta.yaml`, `setup.py` (`>=3.10,<3.12`), CI
  matrix (3.10, 3.11), `pixi.lock` solved on both platforms (Python 3.11; 3.10 also solves, 3.9 does not).
  The `Bio.pairwise2` to `PairwiseAligner` change is separate and is on master (it also runs on Biopython 1.79).
- Open before release: build the base image and run `funannotate test` (predict, update, annotate, compare) in it;
  confirm Trinity, PASA, EVM and bowtie2 source builds on Python 3.11 (not in the lock); `sklearn` `MDS` call in
  `library.py` warns on scikit-learn 1.9 and breaks on 1.10 (pin `<1.10` or port); `distutils` import in `annotate.py`
  (removed in 3.12, so the `<3.12` cap stays).
- Python 3.12 support is a separate item (`distutils`, goatools, Trinity/PASA scripts).

## 2. Split `compare` out of funannotate (decided; the new package is built and tested)

- New package: `funannotate-compare` 0.9.0, https://github.com/hyphaltip/funannotate-compare (stays under `hyphaltip`). Design, sizes and the
  decision log: `plans/COMPARE_EXTRACTION_SPIKE.md` here and `docs/MIGRATION.md` there.
- Done in the new repo: the compare code and the functions it needs, own reference-file `setup`, `outgroups` listing, PHYling staging, pixi
  environment with a `dnds` feature, end-to-end CI (3 genomes; dN/dS `estimate`), BSD replacement for the GPL plot module, `pal2nal.pl`
  unbundled, `CITATION.cff`, `CHANGES.md`. Checked by hand: `--run_dnds full` and `--ml_method raxml`.
- Gate for the funannotate side: funannotate 1.10 removes `compare` only after funannotate-compare 0.9.0 is on **PyPI and Bioconda**, so the
  message can say `conda install funannotate-compare`. Release workflow (tag, build, test the wheel, PyPI trusted publishing, GitHub release) is
  in the new repo's `.github/workflows/release.yml`; the PyPI trusted publisher has to be added on pypi.org by the owner, and the Bioconda recipe
  can be written once the sdist is on PyPI.
- To do in funannotate (1.10 branch), when the gate is met: `compare` and `outgroups` become stubs that print the pointer and exit non-zero;
  same for `database --show-outgroups` and `setup -i busco_outgroups`; delete `compare.py`, `outgroups.py`, `stackedBarGraph.py`,
  `html_template/` and the compare-only names in `library.py` (re-run the call analysis first); drop `compare` from `test -t`, the test data
  table and `image-test.yml`; remove `proteinortho`, `raxml`, `iqtree`, `trimal`, `mafft`, `ete3` from `pixi.toml`, `environment*.yml` and
  `conda-recipe/meta.yaml` after checking `seaborn`, `scikit-learn`, `matplotlib` for other users, and re-lock; remove `docs/compare.rst` and the compare
  and outgroups sections of `commands.rst`, `index.rst`, `databases.rst`; tests for the stubs; `CHANGES.md`.

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
