# Image and installer distribution plan (v1.9.1 line)

Status as of 2026-10-23. Owner: Jason Stajich. Branch: `release/1.9` (kept equal to `master`).

## Goal

Give HPC users (Linux, Apptainer) and Mac users the right artifact, built reproducibly in CI, smaller than
today, and tested by running the tools, not only `--version`.

## Artifacts

| # | Artifact | Platform | Status | Size (compressed) | Audience |
|---|---|---|---|---|---|
| 1 | `funannotate:<tag>` slim (no databases) | linux/amd64 | published; prune not yet in a published image | 4.16 GB measured; about 2.6 GB expected | HPC via Apptainer, cloud, Intel Macs |
| 2 | `funannotate:<tag>-norust` | linux/amd64 | built | similar to #1 | comparing Perl/C++ against Rust code paths |
| 3 | `funannotate:<tag>-db` (databases bundled) | linux/amd64 | published; `chown` double-layer fixed, not rebuilt | 10.46 GB measured; about 7 GB expected | quick start, teaching |
| 4 | `.sif` built from #1 | linux/amd64 | **not built**; `scripts/build_overlay_sif.sh` exists | about the size of #1 | HPC |
| 5 | `funannotate:<tag>` as linux/arm64 | linux/arm64 | **blocked** (see Phase 5) | unknown | Apple Silicon Macs, native speed |
| 6 | One multi-arch tag (#1 + #5) | amd64 + arm64 | depends on #5 | | everyone |
| 7 | `install-funannotate-macos.sh` (native pixi build, not a container) | macOS arm64 | **done**; verified in CI on `macos-14` for `v1.9.1-beta.3` | a few KB script; about 30 min to install | Apple Silicon Mac users |

Not provided: Intel Mac native install (no `osx-64` lock, untested; Intel Macs run the amd64 container natively).
Licensed tools (GeneMark, Phobius, SignalP) are never bundled.

## Decisions

- Canonical image is the slim linux/amd64 one. Databases stay out of it; `-db` is an optional convenience.
- Mac users get the native installer, not a copied environment. A copied pixi environment is not relocatable
  (absolute shebangs, rpaths and links into the build path, editable funannotate install, no activation); the
  installer rebuilds from the lock on the user's Mac and CI runs that same script before attaching it.
- MySQL/MariaDB stay in the Linux image (PASA `--pasa_db mysql`, local runs). The conda `mysqld` is 15-19 MB
  compressed; the apt `mariadb-server` is installed separately in `Dockerfile.base`. Both are tested.
- BLAST+ (`blast` 2.16.0, `lib/ncbi-blast+`) and rmblast (2.14.1, `lib64/ncbi-blast+`) both stay: different versions
  of the same libraries, not duplicates.
- Release tags: `:vX.Y.Z`, `:X.Y`, `:latest` (final releases only). A prerelease never moves `:latest`.
- A different `diamond` pin is acceptable if needed to solve linux-aarch64 (Phase 5). Check first that it does not
  change annotation results.

## Changes from the review (Fable, 2026-10-23)

Verified by us against the repo: items 1, 3 (database path, setup.py dependencies, `Rscript` use) below. The rest is
the reviewer's reading and still to be confirmed.

1. **Release plumbing is a deterministic failure, not a race (Phase 4).** `publish_release.yml:260` dispatches
   `container.yml`, `release.yml`, `release-triggered.yml`, `python-publish.yml` and `macos-release.yml` in one loop, and
   `release-triggered.yml` fails at once if `ghcr.io/...:<tag>` does not exist yet. Fix: take `release-triggered.yml` out of
   the loop and start it from `container.yml` after its verify job (needs `actions: write`), or use `workflow_run`.
   Also: `release.yml` builds a second `:<tag>-slim` image from the same Dockerfile (unsigned), and
   `triggered-build.yml` / `Dockerfile2:5` refer to a `:slim` tag that `container.yml` never pushes.
2. **Pin the base image.** `container.yml` builds `FROM funannotate-base:latest` and `container-base.yml` only pushes
   `:latest`, so a code-only release silently takes whichever base was pushed last and nothing records the digest.
   Push an immutable base tag (for example the `pixi.lock` hash), pin it in `Dockerfile` and `container.yml`, and record
   `org.opencontainers.image.base.digest`. The SLSA predicate has empty `resolvedDependencies`.
3. **The prune test does not yet test what the release does.**
   - `image-prune-test.yml` installs with `pip install --no-deps`; `Dockerfile` resolves `setup.py`'s requirements.
     `setup.py` lists seaborn, psutil, requests and packaging, none of which `pixi.toml` names. After the prune there are
     no compilers or headers, so any requirement pip must build from source will fail. Make the Dockerfile install
     `--no-deps` too and add the missing packages to `pixi.toml`.
   - `funannotate test -t predict` stops at `predict.py:511` unless `$FUNANNOTATE_DB/<busco_db>` exists, and the slim
     image's database directory is empty. The test step needs `funannotate setup -i busco -b dikarya` first.
   - The slim image sets `FUNANNOTATE_DB=/opt/databases/funannotate_db` (`Dockerfile.base:249`) but the `-db` image uses
     `/opt/databases` (`Dockerfile2:28`). Make them consistent.
   - The prune deletes `__pycache__` in a root-owned environment, so every Python start recompiles numpy, pandas and
     sklearn. Run `python -m compileall` at the end of the prune instead.
   - R/Bioconductor (about 570 MB in `lib/R`, more with packages) is kept "for TransDecoder", but no Python in
     `funannotate/` calls `Rscript`. It could go if TransDecoder works without it; test before removing.
   - Add to `image_smoke_tests.sh`: gmap, blat/pblat, fasta36, minimap2, `TransDecoder.LongOrfs`, `emapper.py --version`,
     tRNAscan-SE, CodingQuarry, snap, glimmerhmm, hmmsearch, tbl2asn, exonerate, miniprot, RepeatModeler, and one
     `Launch_PASA_pipeline.pl` run with SQLite.
4. **SIF:** `apptainer pull docker://...@sha256:<digest>` is enough; CI only needs to verify the OCI-to-SIF conversion.
   Check `ENV USER=funannotate` and `LD_LIBRARY_PATH` (`Dockerfile.base:248-251`) under Apptainer.
5. **pixi versions differ:** `Dockerfile.base` uses v0.71.3, the installer and macOS workflow use v0.81.0, and the lock
   is format 7. Confirm 0.71.3 still reads it.

## Phases

### Phase 1 - Prune and measure the Linux image (in progress)

Done:
- `install_scripts/prune_runtime_env.sh`: removes the Rust toolchain, compilers, cmake, headers, JDK-only parts, pandoc,
  docs, static libs, sanitizer libs, qdoc/qmake, and LLVM/clang/LTO libraries nothing else needs. It checks with `ldd`
  that no library goes missing and runs java, mysqld, blastn, python and R afterwards. `PRUNE_TOOLCHAIN=0` keeps the toolchain.
- `pixi_install_trinity.sh` prunes cargo output, `.git`, objects, `trinity_ext_sample_data`, the wiki and bamsifter's htslib tests.
- Measured on the `base` env: 7985 MB -> 5448 MB (first prune), then 7890 MB -> 4810 MB (second prune).
- `Dockerfile2`: the database `chown` is in the same layer as the setup (it stored the 3.15 GB databases twice).
- `pixi.toml`: removed the hardcoded `CONDA_PREFIX`.

In progress:
- `image-prune-test.yml` (run 39984355569): builds the full `Dockerfile.base` without pushing, runs
  `install_scripts/image_smoke_tests.sh` (Trinity on its sample data, conda `mysqld` and apt `mariadbd` with DBD::mysql,
  libraries), installs funannotate from the checkout, runs `funannotate check --show-versions`, then
  `funannotate test -t predict rna-seq`.

Exit criteria: image-prune-test passes, including `predict` and `rna-seq` inside the tolerance windows (1,500-1,800 and
1,630-1,830 genes).

### Phase 2 - Publish the pruned images

- Run `container-base.yml` (pushes `:latest`), then `container.yml` and `release-triggered.yml` for the beta tag.
- Remember `release-triggered.yml` needs the base image to exist first (race seen on `v1.9.1-beta.3`).
- Measure the real compressed sizes and update the table.

### Phase 3 - Apptainer `.sif` (artifact #4)

- Build the `.sif` from #1 in CI and attach it to the release, or document `apptainer pull docker://...` with a tested
  one-liner. Decide which after checking size against the 2 GiB release-asset limit (it will exceed it; a GHCR/ORAS push
  is likely the route, not a release asset).
- Test on the HPC (hpcc): run `funannotate test` from the `.sif` with bind-mounted databases.

### Phase 4 - Release plumbing

- PyPI trusted publisher mismatch (`invalid-publisher`, fails for `v1.9.0` and `v1.9.1-beta.3`): fix the publisher entry
  on PyPI to match `python-publish.yml` with environment `pypi`.
- `publish_release.yml` already starts the container, PyPI and macOS workflows by `workflow_dispatch`; confirm
  `release-triggered.yml` waits for the base/app image (retry or order it after `container.yml`).
- Restructure the release dispatch as in 'Changes from the review' item 1.
- Optional hardening for `funannotate check`: `--strict` exists (exit 1 on any `ERROR:`); do not make it the default
  (macOS legitimately reports missing CodingQuarry/GeneMark).

### Phase 5 - Native arm64 Linux image (artifacts #5, #6)

- The first test solve with `linux-aarch64` failed at `eggnog-mapper` (`diamond >=2.0.11,<2.1`), with `augustus`,
  `transdecoder` and `ete3` conflicts. The review (not yet re-checked by us) says the `diamond` pin is a red herring: the
  lock's eggnog-mapper 2.1.13 build needs only `diamond >=2.0.11`, and aarch64 builds exist for diamond, augustus,
  codingquarry, gmap, salmon, kallisto, openjdk, tbl2asn and repeatmasker (`blat` is linux-64 only).
- The suspected real cause: the `python` pin lives in `[target.linux-64.dependencies]` (`pixi.toml:160-194`), so on
  aarch64 python is unconstrained, bioconda's aarch64 builds want python 3.10+, and eggnog-mapper needs `python <3.12`.
  The `ete3` and `transdecoder` conflicts would follow from that.
- Cheapest attempt: a scratch copy of the repo with `platforms = ["linux-aarch64"]` and a
  `[target.linux-aarch64.dependencies]` block mirroring the `osx-arm64` one (python ~3.10, loosened samtools, fastp,
  stringtie, fasta3), then `pixi lock` (no arm hardware needed). Only if that still fails, try another `diamond` pin.
- If `diamond` or `eggnog-mapper` versions change: results are stable but not bit-identical, and different diamond
  versions between the amd64 and arm64 images could make a shared `FUNANNOTATE_DB` reject existing `.dmnd` files
  (`CheckDiamondDB`, `library.py` around line 610). Keep one diamond version across architectures.
- Build the arm64 base on `ubuntu-24.04-arm`, run the same `image-prune-test` checks, publish one multi-arch manifest.
- Exit criteria: `funannotate test -t predict rna-seq annotate` passes on arm64.

### Phase 6 - Gating test suite (separate, parallel)

Fast (pytest, seconds, every push) and slow (image-level, on demand) gates for decisions where a regression is silent:
- PASA vs BUSCO training choice when the model count is low (boundary cases: just below and above the threshold).
- Augustus species parameters created, found and used; the `--proteinprofile` check.
- GeneMark/BRAKER install works from the image (or fails with a clear message).
- Translation-table support (branch `add_translationtable_support`).
- Output stability: gene-count windows, plus recorded timings to catch slowdowns.
- A local `nf_funannotate1` run against the pushed `prune-test-<variant>` tag (`push_test_tag` input of image-prune-test).

## Open risks

- The prune is validated by `ldd` and a short tool list, not by every code path; the image-level tests are the safety net.
  No compilers in the image means users cannot `install.packages()` or build extensions in it.
- Trinity layer cache is invalidated whenever `pixi_install_trinity.sh` changes (full rebuild, over an hour).
- On an M-series Mac the amd64 image runs under emulation (slow); the native installer (#7) is the recommended Mac path
  until #5 exists.
- GeneMark/Phobius/SignalP availability on arm64 is unknown.
- `docs/conda.rst` still describes the bundle in places; keep it in step with the installer.
