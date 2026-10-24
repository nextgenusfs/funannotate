# Changes

## Unreleased

- Removed `funannotate.library.gb2output`: nothing called it, and it failed with a `TypeError` for any GenBank file with an mRNA feature.
- Python 3.10 or 3.11 is required (3.12 is still not supported), and Biopython is `>=1.87`. Biopython 1.87 and later are
  outside the range of the Bio.Entrez XXE advisory (GHSA-x3vf-39hj-gxr4); funannotate does not use `Bio.Entrez`. Biopython
  1.87 needs Python 3.10. `pixi.lock` is solved with Python 3.11 on linux-64 and osx-arm64.
- Biopython 1.87 removed `nofuzzy_start` and `nofuzzy_end`; every use (`library.py`, `annotate.py`, `update.py`,
  `utilities/contrast.py`) now takes `int(location.start)` and `int(location.end)`. Without this, any GenBank parsing failed with an
  `AttributeError` on Biopython 1.87 or newer.
- `funannotate util gff2prot` and `funannotate util gff-rename` no longer crash: they passed `table=` to `gff2dict`, whose argument is
  `transl_table`, so every run raised a `TypeError`. `--table` now takes effect in both.
- `funannotate annotate`: version checks of emapper and diamond use `packaging.version` through a new `parse_version()` instead of
  `distutils.version.LooseVersion` (distutils is removed in Python 3.12). A version string such as `2.1.4-5-gabc1234` or `v2.0.11` still
  compares by its leading numbers, and a missing or unparseable string counts as 0. The check that refuses emapper versions above
  2.0.0 and below 2.0.5 compared a string with a version object; it now compares two versions.
- `funannotate predict`: the choice between PASA models and BUSCO models for training, the minimum-training-models check and the
  Augustus held-out test set size are now functions in `library.py` with unit tests; behaviour is unchanged. The held-out accuracy of the
  trained Augustus parameters is written to `training_decisions.tsv` (row `augustus_accuracy`); it is a report and does not switch the
  training source. The accuracy table printed after training had the Sensitivity and Specificity headings swapped; they now match the values.
- `funannotate predict --genemark_gtf` with contig names that are not in the genome no longer crashes in the recovery code: the models on
  those contigs are dropped and the contig names are logged as an error.
- `funannotate predict --table` and `--mtable`: an unassigned NCBI table id is reported before any work starts, and a table that reassigns
  stop codons (4, 6, ...) logs a warning that SNAP, GlimmerHMM and CodingQuarry cannot use it. `genetic_codes.is_valid_table` rejects
  fractional ids such as `1.5`.
- `funannotate check --strict --allow-missing a,b` accepts a list of programs that may be missing (for example `gmes_petap.pl,GENEMARK_PATH`
  for images without GeneMark); other missing programs still fail.
- NCBI's native `table2asn` is installed by `install_scripts/install_table2asn.sh` (macOS arm64 and Linux x86_64; it replaces
  `install_table2asn_macos.sh`) and by `Dockerfile.base`. The download's SHA-256 is logged, and `TABLE2ASN_SHA256` pins it. `funannotate
  annotate` already prefers `table2asn` and falls back to `tbl2asn`.
- The development environment (`pixi` default environment) includes `pytest`, so the pytest-based tests run with
  `pixi run -e default python -m pytest tests`. The `base` environment used for the container image is unchanged.
- `pairwiseAlign` in `funannotate update` and `funannotate util contrast` uses `Bio.Align.PairwiseAligner` instead of the
  deprecated `Bio.pairwise2`, so it keeps working when Biopython removes `pairwise2`. The scores are identical, and
  both run on the pinned Biopython 1.79.

## 1.9.1-beta.3 (2026-10-22)

Release 1.9.1-beta.3

## 1.9.1-beta.2 (2026-10-14)

- `funannotate predict --evm_tmpdir DIR` runs EVM in a new directory under `DIR` (for example the node-local
  `$SCRATCH`) instead of `<out>/predict_misc/EVM`. EVM makes one directory per contig, which is slow on a shared
  filesystem when a genome has hundreds of thousands of contigs. Only the final EVM GFF3 is written to the output
  directory. The EVM directory is removed after a successful run and kept, with its path in the log, after a failed
  one. Without the option nothing changes: EVM still runs in `predict_misc/EVM`.
- `funannotate-runEVM.py` passes `--threads` (from `--cpu`), `--skip_empty` and `--concat_out` to the Rust EVM tools
  `recombine_evm_outputs` and `convert_EVM_outputs_to_GFF3`, each only if the installed build lists the option in
  `--help`. Older Rust EVM builds and the Perl engine run as before. With `--concat_out` the converter writes the
  combined GFF3 itself and the script no longer walks the tree to copy every `evm.out.gff3`. Gene order in the
  combined file then follows the partitions listing and can differ from the old filesystem-order walk; the content is
  the same. A zero-byte `evm.out` no longer gets an empty `evm.out.gff3`. An empty `evm.out` can also mean the EVM
  run failed; `--skip_empty` does not change how that is reported. The command builders are in the new
  `funannotate/evm_cmds.py` and are covered by `tests/test_evm_cmds.py`.
- The matching EVidenceModeler changes (parallel convert and recombine, `--concat_out`, `--threads`, `--skip_empty`,
  buffered output, a fixed order that removes a dependence on hash order) are on its `rust_optimize` branch and are not
  yet in a release. The funannotate changes have no effect until an EVM build with them is installed, except
  `--evm_tmpdir`, which works with any build.
- `funannotate predict` runs tRNAscan-SE in parallel chunks when the genome has 2,000 or more contigs with a mean
  length of 100 kb or less, and `--cpus` is 2 or more. tRNAscan-SE spends its time per contig (a temp file and one
  `cmsearch` process per covariance model for each contig), and its `--thread` option only sets `cmsearch --cpu`, so it
  does not help on short contigs. The new code deals whole contigs round robin to one tRNAscan-SE process per CPU
  (`--thread 1`, its own `TMPDIR`), then merges the tables under one header. Contigs are never split, so the tRNA calls
  are the same. On 4,000 random contigs of Massospora_platypediae (4 CPUs) the whole-genome run took 194 s with
  `--thread 4` and 144 s with `--thread 1`; four chunks took 65 s. Chunk temp files go to `$TMPDIR` when it is set.
  Other genomes use the old single call. The code is in `funannotate/trna_chunks.py`, covered by
  `tests/test_trna_chunks.py`.
- `tests/benchmark_evm_convert.sh` times the converter, old build against new, on shared storage and on `$SCRATCH`.
- The pixi route for macOS is merged from `pixi_multiplatform`: the `osx-arm64` platform and the macOS fixes as pixi tasks
  (`install-externals`, `fix-blat-signature`, `install-table2asn`, `macos-smoke-checks`), a smoke-check script, and a weekly
  `macos-smoke` workflow. Tested on Apple Silicon: `check`, `annotate`, `predict`, the smoke checks, and a full RNA-seq run
  (1,709 gene models, 1,708 on the conda route).

## 1.9.0 (2026-10-11)

Release 1.9.0. This is a summary of everything since 1.8.17 (493 commits, 254 files changed, 8
contributors). The per-branch and per-release-candidate details (rc.1 to rc.6) are in the sections
below and in `CHANGELOG.md`.

### Highlights

- **Rust-accelerated externals.** EvidenceModeler, PASA and Trinity are built from Rust-optimized
  forks, selected automatically when present (`FUNANNOTATE_EVM_ENGINE=rust|perl` for EVM) and reported
  by `funannotate check`. Measured effect: see Performance. Each component has a Perl/C++ fallback.
- **Better training data for the ab initio predictors.** RNA-seq evidence-quality gates, complete-ORF
  training sets, one model per locus, single-exon training genes, and an audit log of every training
  decision (`logfiles/training_decisions.tsv`).
- **Genetic code (translation table) support** in all pipeline commands and utilities.
- **Containers.** Signed images on GHCR (a base image with the dependencies and a thin application image),
  Rust and non-Rust variants, built from pinned dependency tags.
- **Installation.** pixi environments, rewritten conda instructions, and an experimental macOS (Apple
  Silicon) route with scripts for the problems found there.
- **A large bug-fix and cleanup pass** over `library.py`, `train.py`, `predict.py`, `update.py`,
  `annotate.py` and the helper scripts (see Fixed and Code cleanup).

### Performance (Rust-optimized tools against the Perl/C++ originals)

Measured on one dataset (*Cordyceps militaris*: a 33.6 Mb genome of 7 sequences, 4.05 million normalized
paired-end RNA-seq read pairs, 25,189 Trinity transcripts), 32 threads, Intel Xeon E5-2683 v4 nodes, one whole
node per run, three replicates per engine, all on the same CPU model. Run-to-run variation within an engine is
under a minute. Trinity transcript counts and PASA gene counts are the same in every run; the final `predict`
gene count varies by about 0.5% between runs (10,596 to 10,652) with either engine.

| step | with Rust | without Rust | speedup | what changes |
|---|---|---|---|---|
| Trinity (read preparation, hisat2 and genome-guided assembly) | 1:05:10 | 1:09:24 | 1.06x (4 min 14 s) | Rust utilities in Trinity's alignment-preparation and coverage scripts (`TRINITY_NO_RUST=1` turns them off) |
| `funannotate predict`, EVM stage | 2 min | 6 min | 4 minutes shorter (stage times have minute resolution) | Rust EVM engine; the whole `predict` step takes 1:02:34 against 1:06:03 (1.06x) |
| PASA stage of `train`, as shipped (Rust `slclust`, `cdbyank` and `faidx`; C++ assembler) | 35:26 | 35:25 | none (1.00x, ranges overlap) | the small Rust tools do not change the time at this scale |
| PASA stage with the Rust assembler (`PASA_ASSEMBLER=pasa_rust`, not the default) | 34:03 | 35:25 | 1.04x (1 min 22 s) | Rust assembler |

The PASA rows were measured on the inputs of a real `train` (the StringTie GTF and splice junctions are
used). A first PASA benchmark that used `train --trinity`, which skips hisat2 and StringTie and gives PASA a
smaller job, gave the same result (1.01x as shipped, 1.05x with the Rust assembler).

Whole `train`, measured end to end (3 runs per engine, all valid): 1:41:32 with the Rust tools (range
1:40:53 to 1:42:03) and 1:46:59 without (1:46:46 to 1:47:03), a speedup of 1.05x (about 5 minutes; the ranges do
not overlap). Put together, `train` takes about 102 minutes with the Rust tools and 107 without, and `predict` 63 minutes against 66, so the Rust tools save about 9 minutes out of about 170 for `train` plus
`predict` (about 5%), with the Rust PASA assembler adding another minute and a half when it is enabled. The
cost of a run is dominated by tools that are not replaced: genome-guided Trinity assembly (61% of `train`),
PASA alignment and assembly, CodingQuarry and protein mapping (21 minutes each in `predict`). The planning
documents' projected speedups for individual Trinity utilities (5 to 50 times) were not measured here.

Stage profile of `predict` (minutes, same for both engines unless noted): protein evidence and Augustus
training 21, Augustus 6, CodingQuarry 21, SNAP 2, GlimmerHMM 9, EVM and filtering 2 (6 with the Perl EVM), final
files 1. Stage profile of `train` (one full run with the Rust tools, 1:43): hisat2 2, genome-guided Trinity 63,
StringTie and alignments 1, PASA alignment and assembly 26, PASA TransDecoder training set 9, kallisto and
model choice 2. In the replicated whole-`train` runs (3 per engine, with and without Rust) the stages were: hisat2 2
minutes in both, genome-guided Trinity 62 (Rust) against 68, StringTie and alignments 1 and 1, PASA alignment and
assembly 26 and 26, PASA TransDecoder training set 9 and 9, kallisto and model choice 1 and 1. The whole-step
saving comes from the Trinity stage.

Equivalence of results: `train` is not deterministic. Two runs of the same engine on the same input differ in
about 4% of PASA feature lines (3,300 to 4,050 of about 81,000 lines, with 9,156 or 9,157 genes in every run), and
two `predict` runs differ in about 1,100 to 1,700 of about 74,500 lines. Differences between engines are 1.0 to
1.3 times the differences between replicates of one engine, with overlapping ranges (PASA: 4,305 to 5,357 against
3,325 to 4,053 differing lines; `predict`: 1,390 against 1,106 and 1,704). With three replicates per engine
this does not show a difference in output between the Rust and the original tools, and it does not rule out a
small one. In the whole-`train` runs (about 81,300 lines), replicates of one engine differ in 564 to 2,008 lines (median 704, 6
pairs) and runs of different engines in 953 to 2,513 (median 2,182, 9 pairs): the median is 3.1 times larger between
engines, the ranges overlap, and every run has 25,188 or 25,189 Trinity transcripts. The replicate counts are too small to
say whether the larger median is a real engine effect; it was not investigated further. The benchmark kit and results are in `tests/benchmark_*` and `benchmark_*.sbatch`.

### Added

- **Genetic codes.** `funannotate/genetic_codes.py` and `--table`/`--mtable` through `predict`, `update`,
  `annotate`, `fix` and the utilities (`gff2tbl`, `gff2prot`, `gff_reformat`, `gbk2parts`), the auxiliary
  scripts (`--gcode`/`--mgcode`), and validation of the table recorded in a parameters file. `annotate` defaults
  the table for `--genbank`/`--gff` input.
- **RNA-seq training.** `--stop_after_trinity` (reuse one Trinity assembly for several strains), fastp trimming,
  a gzip-compressed genome FASTA and gzip-compressed reads, normalized reads compressed in place, PASA options
  (`--pasa_remove_contained`, `--pasa_max_isoforms`, `--pasa_fl_accs`, `--pasa_alt_splice`, opt-in PASA flags),
  the RNA-seq mapping-rate gate, a read-identity measure, a fallback to BUSCO training when PASA gives too few
  transcripts, and `PASACONF` from the environment or the command line.
- **`predict`.** `--min_pasa_complete_models`, `--no_training_single_exon`, `--auto-skip-genemark` for
  fragmented assemblies, `--augustus_BUSCO_timeout`, a standalone tRNA prediction utility, skipping training
  genes with a failed CDS translation, and Braker GFF3 parsing.
- **`annotate`.** DeepTMHMM parsing (an alternative to Phobius and the retired TMHMM for transmembrane regions), tRNAscan-SE 2.x
  pseudo-gene handling, `--allow_ec_without_genename`, level-by-level EC reduction for EggNOG, SignalP 6
  results without a cleavage site, gzip-compressed `--iprscan`, `--antismash` and cached BLAST results.
- **`check`.** A Rust-engine report (EVM, PASA, Trinity utilities), a warning for an incompatible EVM 2.x, and
  names for programs that exist under more than one name (`fasta`/`fasta36`/`fasta3`, `table2asn`/`tbl2asn`).
- **Databases.** The BUSCO, repeat and outgroup downloads are mirrored on a GitHub release (`db-mirror`);
  `setup -u` compares checksums and follows HTTP 308 redirects; the BUSCO setup hint says `-i busco -b <lineage>`.
- **Versions.** `funannotate --version` reports the release, with the git commit when running from a checkout.
- **Release tooling.** `scripts/release_version.py`, a tag guard (`.githooks/pre-push`), and a release workflow.

### Changed

- Training-set selection keeps complete ORFs only, one model per locus (same strand, transitive overlaps), and
  admits single-exon genes with protein support. `funannotate train` no longer runs PASA's alternative-splicing
  analysis by default. `--aligners minimap2` no longer adds blat.
- Product names follow the NCBI protein nomenclature guidelines more closely, with word-boundary rules
  (this fixed substring damage such as `biogenesis`).
- Python 3.6 to 3.11 are supported (3.12 is not). Version parsing uses `packaging`.
- Base image Debian trixie (glibc 2.39 or newer); `mamba` and the Anaconda `defaults` channel are no longer used.
- `funannotate test` reports a failed command and stops the RNA-seq test when `train` fails.

### Fixed

- **minimap2 alignments in PASA (present since 2018).** `bam2gff3` and `bam2ExonsHints` now take exon
  coordinates from the CIGAR string. The old walk never advanced the genome position on mismatches, counted
  soft-clipped bases as aligned, and checked the intron motif only on the last intron. In *N. crassa* only 55
  of 6,973 spliced custom alignments passed PASA's validation before the fix; after it, 9,045 of 14,215 did.
- **PASA fork.** A duplicate-output bug (every GFF3 line written twice) was fixed in PASA v2.6.1-rc.1, and a
  checkpoint name that skipped the failed-alignment BED output was fixed.
- **Resumed `train`.** After a failed PASA step, a resumed run gave PASA the Trinity-only transcripts although the
  imported alignments included the long reads, so it did not reproduce a fresh run.
- **`table2asn`.** The `-Z` flag is a switch, not a file name; the sequencing-center `REFERENCE` feature is no
  longer written for `table2asn` (it made `table2asn` exit with code 2); the preflight check accepts either tool.
- **Deadlocks and crashes.** Two deadlocks in the BUSCO step, an out-of-memory deadlock in soft-mask parsing, a
  `Pool.join()` hang, a crash when no hints are passed to Augustus, a crash when `AUGUSTUS_BASE` was a renamed
  config directory, and the EVM-failure cleanup in `predict`.
- **Silent failures.** Return codes and outputs are now checked after many helper-script subprocess calls
  (EVM partitions, Trinity partitions, InterProScan, Phobius, hmmer, Augustus merge, enrichment), which used to
  leave partial results without an error.
- **Numerical and logic bugs found in review** in `library.py`, `train.py`, `predict.py`, `update.py`,
  `annotate.py`, BUSCO and the helper scripts (for example: Augustus chunk arithmetic that dropped up to 20% of a
  contig, a lexical version comparison of tblastn, `getBestModel` ranking by file order, a locus-grouping key in
  `update`, EggNOG header index 0, and `avg_exon_length` always reported as 0).
- **Downloads.** HTTP 308 redirects (osf.io), the dbCAN and MEROPS locations, and invalid UniProt gzip cache files.

### Code cleanup

- Python 3 compliance in the helper scripts (`queue` import, `cElementTree`, `getchildren`, the reserved word
  `async`, `os.errno`), the removal of dead code (`GAGprotClean`, `download_buscos`, `DBURL`, `busco_links`,
  duplicated `elif` chains), `pOverlap` in constant time, set-based membership tests, an N50 rewrite, quoting
  of file names in shell pipes, and a shared path for the dependency check.
- 267 tests in 26 test files, including real minimap2 and DIAMOND runs, a minus-strand regression test,
  stub-program tests of the benchmark harness, and tests that fail when an install script pins a branch.
- Known issues that were found and not fixed are listed in `CHANGELOG.md`.

### Contributors

Jason Stajich, Jon Palmer, Peter Cock, Adalberto Garcia Garces, Tayab Soomro, and automated contributions
from GitHub Copilot.

### Pinned dependency versions

The source-built tools are pinned to release tags, so a rebuilt image or environment
does not change when a branch is pushed to:

| tool | tag | repository |
|---|---|---|
| Trinity (Rust utilities) | `v2.16.2_rust` | hyphaltip/trinityrnaseq |
| PASA | `v2.6.1-rc.4` | hyphaltip/PASApipeline |
| EVidenceModeler (Rust) | `v2.5.0_rust` | hyphaltip/EVidenceModeler |

- EVM was pinned to the moving branch `rust_optimize`; it is now tag `v2.5.0_rust` (the
  program prints `EVidenceModeler-v2.5.0`). `tests/test_install_pins.py` fails if any
  install script defaults to a branch name.
- `pixi_install_evm.sh` and `pixi_install_pasa.sh` print the commit of a local
  `../EVidenceModeler` or `../PASApipeline` checkout and warn when it is not the pinned
  ref (`PIXI_INSTALL_STRICT=1` makes that an error).

### Docs: conda installation rewritten

- `docs/conda.rst` replaces the old Python 2.7 instructions. A conda or mamba
  environment provides the external tools, `pip install funannotate==1.9.0`
  provides the package, and EVM, PASA and Trinity are built with
  `install_scripts/pixi_install_*.sh`.
- `environment.yml` (Linux) is generated from the `base` environment in
  `pixi.toml`. Regenerate it with `pixi workspace export conda-environment`
  after `[dependencies]` changes, and remove the exported `variables:` block.
- `README.md` and `docs/install.rst` now say that the bioconda package is built separately and may be older than
  this release, and point to `docs/conda.rst` and the container image for 1.9.0.

### Experimental: macOS (Apple Silicon)

A macOS route is documented in `docs/conda.rst`, with `environment-macos.yml` (Python 3.9) and scripts for the problems
found by running it on one Apple Silicon Mac. Tested there: the environment, `funannotate check`, the unit tests, the
EVM, PASA and Trinity builds, `funannotate test -t annotate` and `-t predict` (1,690 genes with the Rust and with the Perl
EVM, 1,500 to 1,800 expected), and the RNA-seq pipeline (`train`, `predict` and `update`: 1,696 and 1,707 genes,
1,630 to 1,830 expected; `train` was resumed after a `pblat` fix, so it was not one uninterrupted run). Not tested:
Intel Macs, `compare`, and a single uninterrupted RNA-seq run on the final code. Not available on macOS arm64: CodingQuarry,
RepeatModeler, legacy `tbl2asn` (NCBI's native arm64 `table2asn` works), and the licensed GeneMark and SignalP.

- `install_scripts/install_augustus_gcc_macos.sh`: `augustus --proteinprofile` crashes in every clang build (bioconda
  included), and `predict` runs that test before BUSCO training. The script builds Augustus 3.5.0 with Homebrew GCC. The cause
  is a `std::map` iterator decremented past `begin()` in `PP::SubstateModel::clearLowScoring`; reported upstream
  (Gaius-Augustus/Augustus#451).
- `install_scripts/install_pblat_shim_macos.sh`: bioconda's arm64 `pblat` crashes with two or more threads, and PASA runs
  it with `-threads=$CPU`. The script runs it with one thread (byte-identical output, seconds); reported upstream
  (icebert/pblat#25).
- PASA's C++ assembler had no `rpath` for `libomp` on macOS and stopped at the first assembly.
  `install_scripts/pixi_install_pasa.sh` adds it (macOS only; the Linux binary is unchanged).
- `install_scripts/pixi_install_trinity.sh`: prepends `SHELL := /bin/bash` without GNU-only `sed -i`, once; on macOS only it passes the
  conda compilers to `make` and defines `_LIBCPP_ENABLE_CXX17_REMOVED_RANDOM_SHUFFLE` for ParaFly. The Linux build is unchanged.
  From `v2.16.2_rust` the Trinity repository does this itself (`trinity-plugins/Makefile` keeps `CXX`/`CC` from the environment
  on macOS, and ParaFly uses `std::shuffle`), so the script no longer passes either.
- PASApipeline `v2.6.1-rc.3` builds `cdbfasta`, `cdbyank` and `psx` with clang (removed the C++17 `register` keyword and added
  missing `int` return types). `v2.6.1-rc.4` links `pasa` with an `rpath` on macOS (`pasa_cpp/Makefile`) and passes
  `-msse4.2` only on x86-64; the install script's rpath repair remains as a fallback for builds made with older tags.
- Python code: pool workers use the `fork` start method (macOS defaults to `spawn`, which re-ran the helper scripts in every
  worker); `fasta`/`tbl2asn` requirements accept `fasta36` and `table2asn`; the REFERENCE feature is not written for `table2asn`;
  the `interpro.tsv` check runs before `annotate` starts.
- A pixi route for macOS (`osx-arm64` platform, tasks for these fixes, a smoke-check script and a weekly workflow) is on the
  branch `pixi_multiplatform` and is planned for 1.9.1.

## 1.9.0-rc.6 (2026-10-03)

Validated against RefSeq on held-out chromosomes of 5 fungal genomes. The design
and data tables are in `docs/training_data_selection_methods.md`; the user guide
is `docs/training.rst`.

### Release tooling: version tags must match `__version__.py`

- `scripts/release_version.py` reads and bumps the release label in
  `funannotate/__version__.py` (`VERSION`, `PRERELEASE`). Release candidates are
  tagged `v1.9.0-rc.N` with `PRERELEASE = "rc.N"`; the final release is tagged
  `v1.9.0` with `PRERELEASE = ""` and reports `1.9.0`.
- `.githooks/pre-push` refuses to push a `v*` tag that does not match the label
  at the tagged commit. Enable it once per clone with
  `git config core.hooksPath .githooks`.
- `container.yml` runs the same check before building a tagged image.
- `publish_release.yml` (CD/Release) bumps `__version__.py`, `CITATION.cff` and
  `CHANGES.md`, tags, tests, creates the GitHub release and starts the image and
  PyPI workflows. Bumps: `prerelease` (rc.4 -> rc.5; `preid rc` moves beta -> rc.1),
  `final`, `patch`/`minor`/`major`, `prepatch`/`preminor`/`premajor`.

### Feature: evidence-quality gates

**`funannotate/train.py`**
- `--min_rnaseq_map_rate` (default 10) and `--rnaseq_gate_reads` (default 200000):
  maps sampled reads to the genome with minimap2 (`-x splice:sr`); below the rate,
  train exits with code 3 before Trinity/PASA. Report `logfiles/train_rnaseq_gate.tsv`.

**`funannotate/predict.py`**
- `--min_pasa_complete_models` (default 500): counts complete-ORF models in
  `--pasa_gff`; below it, predictors that would train from PASA train from BUSCO
  (PASA stays EVM evidence). Applies when any predictor trains from PASA; skipped
  on resume or with `--augustus_gff`. Report `logfiles/predict_training_gate.tsv`.

### Feature: training-decision audit log

**`funannotate/library.py`, `train.py`, `predict.py`**
- `set_training_decision_log`, `record_training_decision`, `training_decision_summary`.
  Every choice about training data is written to `logfiles/training_decisions.tsv`
  (command, stage, decision, value, threshold, outcome, reason, timestamp), echoed as
  `TRAINING-DECISION ...` log lines, and summarized in a table before Augustus training.

### Feature: single-exon training genes (default on)

**`funannotate/library.py`, `predict.py`**
- `selectTrainingModels(..., single_exon_align, single_exon_share)` admits complete
  single-exon PASA models with ≥80% same-strand protein2genome coverage, capped at
  share/(1−share) × multi-exon models. The share comes from GeneMark-ES models
  (`model_single_exon_share`), with protein alignments as fallback
  (`protein_single_exon_share`). `--no_training_single_exon` turns it off.

### Change: training-set selection and one model per locus

**`funannotate/library.py`**
- `selectTrainingModels` keeps only complete ORFs (`is_complete_model`, protein-based,
  correct on both strands) and removes overlaps transitively (`cluster_overlapping`).
- Returns 0 instead of exiting when no model survives filtering.
- New helpers: `count_complete_orf_models`, `is_complete_cds`, `pasa_gate_applies`.

**`funannotate/train.py`**
- `getBestModel` keeps one model per same-strand transcript-span locus (transitive
  clusters), ranked by complete ORF, CDS exons, CDS length, TPM (`pick_locus_model`;
  optional guarded ranking via `complete_min_frac`).
- `pasa_feature_flags`, `--pasa_unspliced_join_spliced`, `--pasa_one_alignment_per_cdna`:
  opt-in PASA flags, passed only when the installed launcher supports them.

### Tests
- `tests/test_training_gates.py` (30 tests) and `tests/test_training_selection.py`
  (33 tests), including a minus-strand regression test and real minimap2/DIAMOND runs.

## Branch: enable_augustustimeout_busco

- Pull request #1173 also fixed two deadlocks in the BUSCO step that hang funannotate predict (@atiweb)

### Fix: Augustus jobs in BUSCO training can now be killed if they hang

Individual Augustus jobs called during BUSCO-based ab-initio training had no timeout,
so a single stuck job would block its worker thread indefinitely and hang the entire run.

**`funannotate/aux_scripts/funannotate-BUSCO2.py`**
- Added `--augustus_timeout` argument (default 300 s). Each per-gene Augustus job is now
  killed with SIGKILL (via `os.killpg` on its own process group) if it exceeds this limit,
  and the run continues with that gene model skipped.
- Per-job elapsed times are collected across all worker threads and logged as a summary
  after the Augustus pass (`mean`, `median`, `max`, timeout count) so the threshold can
  be calibrated empirically.
- `_extract` calls are now wrapped in try/except so a truncated output file left by a
  killed job does not abort the protein-extraction step.

**`funannotate/predict.py`**
- Added `--augustus_BUSCO_timeout SECONDS` argument (default 300) that is forwarded to
  `funannotate-BUSCO2.py` as `--augustus_timeout`.

## Branch: auto_skip_genemark

### Feature: `--auto-skip-genemark` flag for `funannotate predict`

Adds a new `--auto-skip-genemark` argument to `funannotate predict`. When set, if the
assembly appears highly fragmented (longest scaffold < 50 kb), GeneMark weight is
automatically reset to 0 and GeneMark execution is skipped with a warning, rather than
logging an error and continuing with a likely-failing GeneMark run.

Without the flag the existing error message is shown (now augmented with a hint to use
`--auto-skip-genemark`), and behaviour is unchanged.

**`funannotate/predict.py`**
- Added `--auto-skip-genemark` argument (store_true, default False).
- In the fragmented-assembly check: when `--auto-skip-genemark` is set, logs a warning
  and sets `StartWeights["genemark"] = 0`; otherwise logs the existing error.
- GeneMark execution block now checks `StartWeights["genemark"] == 0` first and skips
  all GeneMark steps when the weight was zeroed (covers both `--auto-skip-genemark` and
  any explicit `--weights genemark:0` passed by the user).

### Fix: `table2asn -Z` flag treated as boolean, not file-path argument

`table2asn` (the NCBI replacement for `tbl2asn`) takes `-Z` as a boolean switch that
writes the discrepancy report alongside the output as `<input>.dr`; it does **not**
accept a filename after `-Z`. The legacy `tbl2asn` still takes `-Z <file>`. The previous
code always passed `-Z discrepency.report.txt`, which caused `table2asn` to interpret the
filename as a positional argument and silently misbehave.

**`funannotate/library.py`** (`build_tbl2asn_cmd`)
- When `dialect == "table2asn"`, appends only `-Z` (no filename).
- When `dialect == "tbl2asn"`, appends `-Z <discrepancy_file>` as before.

**`funannotate/aux_scripts/tbl2asn_parallel.py`** (`tbl2asn_runner`)
- Same dialect split: `table2asn` gets `-Z` alone; `tbl2asn` gets `-Z <file>`.
- Added `DEBUG` print of the final command to aid future troubleshooting.
- Discrepancy-report collection now also matches files ending in `.dr` (the filename
  `table2asn` produces) in addition to the legacy `discrepency.report.txt`.

---

## Branch: stoptrain_after_trinity

### Feature: `--stop_after_trinity` flag for `funannotate train`

Adds a new `--stop_after_trinity` argument to `funannotate train` that exits cleanly
after the Trinity genome-guided assembly step, before PASA runs.  This supports
pipelines that process multiple strains per species and need to reuse a single Trinity
assembly across annotation runs.

**`funannotate/train.py`**
- Added `--stop_after_trinity` argument (store_true).
- After Trinity transcripts are confirmed/built, logs the output path and calls
  `sys.exit(0)` when the flag is set.

**`funannotate/funannotate.py`**
- Added `--stop_after_trinity` to the `train` subcommand help text.

### Feature: Normalized FASTQ files compressed in place after Trinity normalization

After `in silico` read normalization, the output `.norm.fq` files are now gzip-compressed
in place (via `lib.Fzip_inplace`) to reduce disk usage.  Cache detection on re-runs is
updated to look for `.norm.fq.gz` instead of symlinks to plain `.norm.fq`.

**`funannotate/train.py`** (`runNormalization`)
- After normalization, `left.norm.fq`, `right.norm.fq`, and `single.norm.fq` are each
  compressed to `.gz` and the function returns the `.gz` paths.
- Re-run checks replaced `os.path.islink()` with `lib.checkannotations()` against the
  `.gz` paths.

## Branch: eggnog_geneprod_issue

### Fix: `funannotate --version` always reported base version in egg installs

`_git_version()` in `funannotate/__version__.py` returned the static base version
(`1.8.17`) whenever `.git` was absent — which is always the case inside a built egg or
wheel.  The `_version.txt` fallback described in earlier CHANGES was never implemented.

**`funannotate/__version__.py`**
- Added `_version.txt` read between the `.git` check and the `_base` fallback: if
  `funannotate/_version.txt` exists, its contents are returned as the version string.

**`setup.py`**
- Writes the resolved version to `funannotate/_version.txt` immediately after
  `exec()`-ing `__version__.py`, so every `pip install` / `python setup.py bdist_egg`
  bakes the correct version into the package.
- Added `package_data={"funannotate": ["_version.txt"]}` so the file is included in
  eggs and wheels.

### Feature: IPRscan XML and UniProt BLAST XML automatically compressed after use

Large intermediate XML files are now gzip-compressed in place after they are fully
consumed, reducing disk usage in `annotate_misc/`.  Subsequent runs detect the `.gz`
cache and skip re-running the expensive tools.

**`funannotate/annotate.py`** (`SwissProtBlast`)
- After parsing `uniprot.xml`, the plain file is compressed to `uniprot.xml.gz` and
  removed.  The existing `.gz` cache check (from the previous feature) means diamond
  is not re-run if `uniprot.xml.gz` already exists.

**`funannotate/annotate.py`** (IPRscan block)
- After `iprscan2annotations.py` runs, `iprscan.xml` is compressed to `iprscan.xml.gz`
  and removed (skipped when the input was already a user-supplied `.gz` file).
- On subsequent runs with no `--iprscan` argument, the block now checks for
  `iprscan.xml.gz` if the plain `iprscan.xml` is absent, so parsing is not repeated.

### Feature: IPRscan XML and UniProt BLAST XML support gzip-compressed files

**IPRscan XML (`--iprscan`)**

`funannotate/aux_scripts/iprscan2annotations.py` opened the XML with plain `open()`,
rejecting `.gz` inputs.  `annotate.py` copied the file to `iprscan.xml` before parsing,
so pre-computed `.gz` results also couldn't be passed directly.

- `iprscan2annotations.py`: added `import gzip`; the XML file is now opened with
  `gzip.open(..., "rt")` when the path ends in `.gz`, plain `open()` otherwise.
- `annotate.py` (`--iprscan`): when the supplied file ends in `.gz`, `IPRCombined` is
  pointed directly at the source (no copy step), so `iprscan2annotations.py` receives
  the `.gz` path and decompresses on the fly.  Also fixed the `samefile` guard to
  pre-check `os.path.isfile` so it doesn't raise `FileNotFoundError` on a first run.

**UniProt BLAST XML (`uniprot.xml`)**

`SwissProtBlast` in `annotate.py` cached blast results to `uniprot.xml` and re-read
them with plain `open()`.  If the cached file was manually compressed to `uniprot.xml.gz`
to save space, the cache was missed and blast was re-run.

- The cache check now also looks for `uniprot.xml.gz`; if found, blast is skipped.
- The result file is opened with `gzip.open(..., "rt")` when `.gz` is present, plain
  `open()` otherwise, so `SearchIO.parse` receives a text stream in both cases.

### Feature: `--antismash` accepts gzip-compressed input (`.gbk.gz`)

`ParseAntiSmash` and `antismash_version` in `funannotate/library.py` previously opened
the antiSMASH GenBank file with `open(..., "r")`, rejecting compressed inputs.

**`funannotate/library.py`**
- Added `_open_maybe_gzip(path, mode)` helper: returns `gzip.open(path, "rt")` for
  `.gz` paths and plain `open(path, mode)` otherwise.
- `antismash_version` and `ParseAntiSmash` now call `_open_maybe_gzip` instead of
  `open`, so both `.gbk` and `.gbk.gz` inputs are handled transparently.

**`funannotate/annotate.py`**
- When `--antismash` is a `.gz` file, `antismash_input` is set to the supplied path
  directly (no copy/decompression step needed since parsing is now gzip-aware).
- Also fixed the `os.path.samefile` guard to pre-check `os.path.isfile(antismash_input)`
  so it does not raise `FileNotFoundError` on a first run when the destination doesn't
  exist yet.

### Fix: `SameFileError` when `--signalp/--phobius/--eggnog/--iprscan/--antismash` points to the existing output file

`shutil.copyfile(src, dst)` raises `SameFileError` when src and dst resolve to the same
inode.  This happened when users re-ran `funannotate annotate` in the same output
directory and passed the already-generated results file as the `--signalp` (or other)
argument.  The `phobius`, `eggnog`, and `antismash` variants were also at risk of
silently deleting the file (they called `os.remove(dst)` before copying).

**`funannotate/annotate.py`**
- All five pre-computed result copy sites (`--signalp`, `--phobius`, `--eggnog`,
  `--iprscan`, `--antismash`) now guard with `os.path.samefile(src, dst)` and skip
  the copy when src and dst are the same file.  The `--iprscan` string-equality guard
  is replaced with the reliable `samefile` check.

---



### Fix: Python 3 compliance in `funannotate/aux_scripts/`

Four scripts contained Python 2-only constructs that would fail on Python 3.9+.

**`funannotate-BUSCO2.py`**
- Fixed reversed `queue` import: `import queue.Queue as queue` (invalid syntax) was
  the primary branch, with the actual Python 3 module as the fallback. Corrected to
  `import queue` (Python 3) / `import Queue as queue` (Python 2) order.

**`xmlcombine.py`**
- `from xml.etree import cElementTree` — `cElementTree` was deprecated in Python 3.3
  and removed in Python 3.9. Changed to `from xml.etree import ElementTree as cElementTree`.

**`iprscan-local.py`**
- `Element.getchildren()` was deprecated in Python 3.8 and removed in Python 3.9.
  Replaced `data.getchildren()` with `list(data)`.

**`runIPRscan.py`**
- `urllib2.__version__` referenced the Python 2 `urllib2` module, which is never
  imported in the Python 3 version of this file. Replaced with `platform.python_version()`
  (`platform` was already imported).
- `--async` / `options.async`: `async` is a reserved keyword in Python 3.7+, causing
  a `SyntaxError` at compile time. Renamed the flag to `--asyncjob` (dest `asyncjob`).

---



### Fix: EggNog product descriptions leaving dangling prepositions (e.g. "to Saccharomyces cerevisiae gds1")

EggNog mapper descriptions such as "Homolog to Saccharomyces cerevisiae gds1" were being
cleaned to "to Saccharomyces cerevisiae gds1" because the `\bhomolog\b` rule removed only
the word "homolog" without consuming the trailing preposition "to". The same latent problem
existed for "homologue of …" and "Ortholog to …" forms.

**Fix:** the homolog/homologue/ortholog patterns now first match the trailing `to`/`of`
phrase (`\bhomou?logs?\s+(?:to|of)\s+`), then a bare-form fallback removes any remaining
standalone occurrence. The result for the example above is "Saccharomyces cerevisiae gds1".

### Model-organism descriptions converted to "Protein <gene>" form

EggNog and UniProt descriptions that refer to a characterised gene in a model organism
(e.g. "Homolog to Saccharomyces cerevisiae gpr1") are now converted to a clean
"Protein \<gene\>" form (e.g. "Protein gpr1") rather than leaving a garbled remnant.

**Recognised model organisms** (defined in `_MODEL_ORGS`):
- *Saccharomyces cerevisiae*
- *Candida albicans*
- *Aspergillus nidulans*
- *Aspergillus fumigatus*
- *Cryptococcus neoformans*

**Capture patterns** (applied before bare-removal fallbacks):

| Input pattern | Output |
|---|---|
| `Homolog(ous) to <organism> <gene>` | `Protein <gene>` |
| `Homologue of <organism> <gene>` | `Protein <gene>` |
| `Ortholog(ous) to/of <organism> <gene>` | `Protein <gene>` |
| `similar to <organism> <gene>` | `Protein <gene>` |

When the organism is not in the recognised list (e.g. *Neurospora crassa*), the homolog
phrase is stripped and the organism + gene token are left as-is for downstream curation.

**"homologous recombination" preserved**: bare `\bhomologous\b` is intentionally *not*
removed because "homologous recombination" is a valid biological term. Only the
`homologous to/of <X>` phrase form is removed.

### Expanded product-name cleaning rules per NCBI International Protein Nomenclature Guidelines

The `rep` substitution table in `funannotate annotate` (`annotate.py`) has been substantially
revised to align with the
[NCBI International Protein Nomenclature Guidelines](https://www.ncbi.nlm.nih.gov/genbank/internatprot_nomenguide/).

**Rules removed**

| Old rule | Reason for removal |
|---|---|
| `\buncharacterized\b` → `putative` | "uncharacterized protein" is a valid NCBI product name; only the British spelling variant (`uncharacterised`) is now corrected |
| `\bfamily\b` → `""` | NCBI recommends "X family protein" as the correct format for family-similarity descriptions |

**Rules fixed**

| Rule | Change | Reason |
|---|---|---|
| `\bhomolog\b` / `\bhomologue\b` → `""` | Replaced with compound `\bhomou?logs?\s+(?:to\|of)\s+` first, then bare form | Prevents dangling preposition |
| `\binactivated\b` → `" "` | Now → `"inactive"` | NCBI: "inactive" is reserved for proteins with altered catalytic residues; "inactivated" implies a process |
| `\bgene\b` → `"protein"` | Now anchored `\bgene\b\s*$` (terminal only) | Prevents corruption of "gene expression regulator", "gene ontology", etc. |
| Organism patterns | Consolidated and made case-insensitive; added `human`, `murine`, `mouse` | Consistency |

**New rules added**

*Multi-word phrases (applied before single-word rules):*
- `similar to` → `""` (NCBI: eliminate)
- `ortholog(ue) (to/of)` → `""` then bare ortholog form (same fix as homolog)
- `protein of unknown function` → `"hypothetical protein"`
- `conserved hypothetical` / `hypothetical conserved` → `""`
- `cell surface`, `surface antigen` → `""`
- `involved in`, `implicated in`, `also known as` → `""`

*Single terms (NCBI "eliminate entirely" list):*
`novel`, `fragment`, `partial`, `truncated`, `dubious`, `doubtful`, `expressed`,
`secreted`, `WGS`, `ORF`, `GO`

*`conserved` as standalone qualifier* — removed unless immediately followed by `domain`
(negative lookahead `(?!\s+domain)` preserves "conserved domain-containing protein").

*Post-substitution cleanup:*
- Strip leading prepositions (`to`, `of`, `by`, `with`, `from`, `for`) left by removals
- Collapse `protein protein` → `protein`
- `-ase protein` → `-ase` (NCBI: enzyme names ending in -ase do not need trailing "protein")

*Greek letter normalisation (NCBI: spell out in lowercase):*
Alpha → alpha, Beta → beta, Gamma → gamma, Kappa → kappa, Sigma → sigma.
Delta is intentionally excluded (capitalised in steroid/fatty acid formal names).

*Arabic numerals preferred over informal Roman numerals:*
type I–IV → type 1–4. Established formal names such as "RNA polymerase II" are not affected.

### Dynamic git version tag in `funannotate --version` and "Running …" log lines

`funannotate --version` and all "Running funannotate v…" log lines previously reported the
static version baked into installed package metadata (`importlib.metadata.version()`),
ignoring local git state.

**`funannotate/__version__.py`**
- Added `_git_version()` which calls `git describe --tags --dirty --always --long` to
  produce PEP 440 strings such as `1.8.17.dev107+gfca1081.dirty`.
- Priority: git describe (when `.git` present) → `_version.txt` (baked at install) → static base.
- Short-circuits immediately if no `.git` directory is present (zero subprocess cost for
  installed/production environments).

**`funannotate/_version.txt`** (generated, not committed)
- Written by `setup.py` at `pip install` time so wheel/sdist installs report the exact
  version string without needing git at runtime.

**`funannotate/funannotate.py`**
- Replaced `importlib.metadata.version()` with a direct import from `funannotate.__version__`.
- Restored `import importlib` (still needed for `importlib.import_module` dispatch at line 712).

**`funannotate/library.py`** (`get_version`)
- Replaced `importlib.metadata.version()` with import from `funannotate.__version__`
  so all "Running …" log lines show the full git-describe version string.

**`setup.py`**
- Pre-populates `__file__` in the `exec()` context so `_git_version()` can locate `.git`
  during `pip install`; writes the resolved version to `funannotate/_version.txt`.

### Debug logging for product-name cleaning

Added `lib.log.debug()` calls throughout the product-name cleaning block in `annotate.py`
so that running `funannotate annotate --debug` now logs every substitution that fires,
making it straightforward to trace where any description is being transformed.

---

## Branch: bug1160_annotationregex_fix

### Fix: product-name cleaning over-aggressively replaced substrings (bug #1160)

The `rep` substitution table in `funannotate annotate` (annotate.py) was built without
word-boundary anchors, so every key matched as a bare substring anywhere in a product
description. The most visible symptom was `"gene"` → `"protein"` corrupting words that
merely *contain* `gene` as a substring (e.g. `biogenesis` → `bioproteinsis`,
`xenogeneic` → `xenoproteinlic`). The same latent bug affected `"EC"` (→ ECM, PECAM),
`"frame"` (→ frameshift, framework), `"homolog"` (→ homologous), and `"COG"` (→ COGnate).

**Fix:** replaced the single compiled-regex dict approach with a list of
`(compiled_pattern, replacement)` tuples, each using explicit `\b` word-boundary
anchors so substitutions only fire on whole words.

**`funannotate/annotate.py`** (`annotate.py:1202–1229`)
- Rewrote the product-name cleaning block to use `\b`-bounded patterns; the old
  `re.escape(k)` dict lookup trick is removed in favour of sequential `_pat.sub()` calls.

### New flag: `--allow_ec_without_genename` in `funannotate annotate`

`tbl2asn` rejects an `EC_number` qualifier on a CDS that has no `gene` qualifier (i.e.
no gene name was assigned). By default funannotate now suppresses EC_number output for
any gene model that did not receive a gene-name annotation from UniProt or EggNog.

Pass `--allow_ec_without_genename` to restore the previous behaviour and write EC numbers
regardless of gene-name assignment.

**`funannotate/annotate.py`**
- Added `--allow_ec_without_genename` flag (default off).
- `lib.updateTBL()` is called with `require_gene_for_ec=not args.allow_ec_without_genename`.

**`funannotate/library.py`** (`updateTBL`)
- Added `require_gene_for_ec=True` keyword argument.
- When `True`, `EC_number` qualifiers are skipped for any locus whose annotation dict
  has no `"name"` entry (i.e. no gene name assigned).

---

## Branch: add_translationtable_support

### Translation table support — consistent `transl_table` naming

Standardized the NCBI genetic code / translation table parameter across all modules.
The canonical keyword argument is now `transl_table` throughout, matching the GenBank
`transl_table` qualifier name. All functions and their call sites have been updated.

**`funannotate/library.py`**
- `translate()`: renamed parameter `table` → `transl_table`
- `extend2stop()`: updated internal `translate()` calls to use `transl_table=`
- `convertgff2tbl()`: fixed `NameError` — `table_table` was undefined; corrected to `transl_table`
- `tbl2dict()`: fixed `translate()` call keyword from `table=` to `transl_table=`
- `zff2gff3()`: updated two `translate()` calls to use `transl_table=`
- `translatemRNA()`: updated `translate()` call to use `transl_table=`
- `gff2dict()`: calls to `translate()` already used `transl_table=`; now correct after rename above

**`funannotate/predict.py`**
- `lib.RunGeneMarkES()` / `lib.RunGeneMarkET()` calls: `table=` → `transl_table=`
- `lib.gff2dict()` call: `table=` → `transl_table=`
- `lib.GFF2tbl()` call: `table=` → `transl_table=`
- `lib.tbl2allout()` call: fixed syntax error and typo (`transl_tabel =` with no value) → `transl_table=args.table`

**`funannotate/annotate.py`**
- `lib.convertgff2tbl()` call: `table=` → `transl_table=`
- Two `lib.gb2parts()` calls: added missing `transl_table=args.table`
- `lib.tbl2allout()` call: added missing `transl_table=args.table`

**`funannotate/update.py`**
- `gff2pasa()`: renamed parameter `table` → `transl_table`; fixed internal `lib.gff2dict()` call
- `GFF2tblCombinedNEW()`: renamed parameter `table` → `transl_table`; fixed `lib.gff2interlapDict()` calls, added `transl_table=` to `lib.translate()` call, fixed `lib.dicts2tbl()` call
- `gff2interlap()`: renamed parameter `table` → `transl_table`; fixed internal `lib.gff2dict()` call
- All call sites of the above local functions updated from `table=` to `transl_table=`
- `lib.tbl2allout()` call: added missing `transl_table=args.table`

**`funannotate/fix.py`**
- `lib.tbl2allout()` call: added missing `transl_table=args.table`

**`funannotate/utilities/gbk2parts.py`**
- `lib.gb2parts()` call: `table=` → `transl_table=`

### Fix: Product annotations lost during `funannotate update` PASA reannotation

When updating an annotation from a GenBank file, the `product` qualifier (e.g.
"cytochrome P450", "hypothetical protein") was parsed correctly from the GBK input
but silently dropped during the PASA reannotation step. `getBestModels()` rewritten
mRNA lines without preserving the `product=` GFF3 attribute, so `GFF2tblCombinedNEW()`
defaulted every gene to "hypothetical protein" in the TBL passed to tbl2asn/table2asn.
The final GFF3, GBK, proteins, and transcripts all inherited the lost products.

**Fix:** product annotations are now persisted to a TSV file
(`update_misc/products.tsv`) during input parsing, then merged back into the gene
dictionary immediately before the TBL is written.

**`funannotate/update.py`**

- **`saveProductTable(genes, productfile)`** (new, line 110): writes
  `locus_tag<TAB>product` for each gene whose product is not "hypothetical protein".
  Called from both `gbk2pasaNEW()` and `gff2pasa()` when `productfile` is supplied.

- **`gbk2pasaNEW()`** (line 126): added `productfile=None` parameter; calls
  `saveProductTable()` before returning so that original products from the GBK file
  are captured.

- **`gff2pasa()`** (line 360): added `productfile=None` parameter; calls
  `saveProductTable()` before returning so that original products from a GFF3+FASTA
  input are captured.

- **`GFF2tblCombinedNEW()`** (line 1747): added `productfile=None` parameter. After
  the gene-renaming loop and before `lib.dicts2tbl()`, reads `products.tsv` and
  overwrites `product` entries in `renamedGenes` for any locus tag present in the
  table. Genes not in the table (e.g. novel PASA models) retain their existing
  "hypothetical protein" default.

- **`main()`**: defines `productfile = os.path.join(tmpdir, "products.tsv")` and
  passes it to `gbk2pasaNEW()`, `gff2pasa()`, and `GFF2tblCombinedNEW()`.

**Data flow:**

```
GBK/GFF3 parsed → products saved to update_misc/products.tsv
                                                         │
PASA runs → getBestModels drops products from GFF3 ◄──────┤
                                                         │
GFF2tblCombinedNEW generates TBL → products merged back ─┘
                                                         │
tbl2asn converts TBL → GBK ◄─────────────────────────────┘
tbl2allout converts TBL → final GFF3/proteins/transcripts
```
