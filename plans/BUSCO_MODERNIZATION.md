# BUSCO modernization for funannotate 1.10 (vendored BUSCO2 vs upstream BUSCO v5/v6)

Worktree reviewed: `/Users/jstajich/projects/funannotate-1.10`, branch `release/1.10-python310`
(HEAD 58f5f21, `compare`/`outgroups` already removed). Read-only review; no tracked file changed.
Date: 2026-10-28. Supersedes the "(verify)" items of `docs/BUSCO_REVIEW.md`.

Legend. **Code facts** carry `file:line` from this worktree. **Upstream facts** (section 1) were
checked against the BUSCO GitLab repo, user guide, bioconda recipe, the busco-data index and the
OrthoDB pages, or by a local experiment (section 7); each row says what was confirmed and what was not.
Nothing from memory is presented as verified.

## 0. Answer in brief

- The vendored script (`funannotate/aux_scripts/funannotate-BUSCO2.py`, 3,462 lines, `VERSION = '2.0'`
  at :2990; py2 copy 3,340 lines) does two different jobs: (a) **training-set generation** for Augustus,
  SNAP and GlimmerHMM in `predict` (genome mode, tblastn + `augustus --proteinprofile` per marker, two
  passes with retraining), and (b) **functional annotation** in `annotate` (proteins mode, hmmsearch
  only, writes `note BUSCO:<id>` qualifiers into the GenBank/annotation tables).
- Job (b) can be replaced by upstream BUSCO or, more cheaply, by keeping the vendored proteins mode and
  swapping the lineage data. Job (a) cannot be replaced by upstream BUSCO in 1.10 without losing
  behaviour the project has invested in during 2026 (per-job Augustus timeouts, two deadlock fixes,
  `-db_gencode` for alternative genetic codes, the run-local Augustus config, the protein re-validation
  pass, the ploidy>1 duplicate selection) and without a benchmark that the methods draft
  (`docs/training_data_selection_methods.md`) would need to be re-run against.
- The lineage data is the easier and more valuable modernization: `fungi_odb10` from busco-data.ezlab.org
  has the same layout the vendored script reads (`hmms/`, `prfl/`, `scores_cutoff`, `lengths_cutoff`,
  `ancestral*`, `dataset.cfg`) and the vendored script runs on it unchanged (local test, section 7:
  94.3 % Complete on the dataset's own ancestral sequences, 6 s). **odb12 datasets do not work with the
  vendored script** (no `lengths_cutoff`; HMM `NAME` lacks the `at<taxid>` suffix; section 7).
- The BUSCO ids change namespace (odb9 `EOG0926...` -> odb10/odb12 `<n>at4751`), which breaks
  comparison with every existing `funannotate annotate` output and the outgroup proteomes of
  `funannotate-compare` (keyed by `EOG` ids), and `funannotate-compare outgroups` imports
  `aux_scripts/funannotate-BUSCO2.py` from the installed funannotate package
  (`funannotate-compare/src/funannotate_compare/outgroups.py:9-20, 37-60, 154-166`). Deleting the script
  in 1.10 breaks that package.
- Upstream state (section 1.1): BUSCO 6.1.0 (2026-06) still has the Augustus genome mode with
  retraining output (U7, U8), so the earlier review's "v6 drops Augustus" is wrong (U33); miniprot has
  been the default since 5.7.0 (U10); there is no genetic-code flag for eukaryotic genome mode and the
  miniprot call cannot take one (U11); Augustus mode needs a writable `$AUGUSTUS_CONFIG_PATH/species`
  (U9); the bioconda package is noarch with `python >=3.3` and a `sepp==4.5.5` pin (U4) and is not in
  the pixi lock (L9). The default dataset generation is now odb12.2 (U19), whose fungal tarballs are
  0.4-1.6 GB (U28).
- **Recommendation**: 1.10 = hygiene (section 6) + odb10 lineage support behind the existing
  `--busco_db` name space (new names `fungi_odb10` etc., odb9 names kept) + an id-namespace tag in the
  annotation note; keep the vendored script as the engine for both jobs. 1.11 = optional external BUSCO
  for `annotate` (proteins mode, any version) and a `--busco_results` import for `predict`, after the
  benchmark in section 4.4. Removal of the vendored script is not a 1.10 or 1.11 item.

## 1. Verified fact table (upstream BUSCO)

See section 1.1 for the web-verified rows and section 1.2 for the rows verified by local experiment.
Status: **confirmed** (primary source fetched and quoted), **partly** (source fetched but the exact
claim is not stated there), **unconfirmed** (could not be established from a primary source).

### 1.1 Upstream rows (web sources)

Sources (abbreviated): GL = https://gitlab.com/ezlab/busco (`-/tags`, `-/raw/master/CHANGELOG`,
`-/raw/master/setup.py`, `-/raw/master/src/busco/run_BUSCO.py`, `.../busco_tools/{augustus,miniprot,metaeuk,hmmer}.py`,
`.../BuscoDownloadManager.py`); UG = https://busco.ezlab.org/busco_userguide.html;
FV = https://busco-data.ezlab.org/v5/data/file_versions.tsv; LIN = `https://busco-data.ezlab.org/v5/data/lineages/<name>.tar.gz`;
BC = https://raw.githubusercontent.com/bioconda/bioconda-recipes/master/recipes/busco/meta.yaml (6.1.0) and commit
34e08fb of the same file (5.8.3); AC = https://anaconda.org/bioconda/<pkg>/files; IUC = real v6.1.0 test output at
https://raw.githubusercontent.com/galaxyproject/tools-iuc/main/tools/busco/test-data/genome_results_miniprot/full_table.

| # | Claim (from docs/BUSCO_REVIEW.md or the maintainer's questions) | Status | Evidence | Source |
|---|---|---|---|---|
| U1 | Latest releases: v6.1.0 (2026-06-05), v6.0.0 (2025-07-09), v5.8.3 (2025-03-07); 5.8.0 2024-10-11, 5.7.0 2024-03-15 | confirmed | tags page; `_version.py` = 6.1.0 | GL/-/tags |
| U2 | CHANGELOG 5.5-6.1: 5.5.0 miniprot pipeline (beta); 5.6.1 `--update-data` removed (default; `--offline` to disable); **5.7.0 miniprot default for eukaryote genomes, `--metaeuk` optional, GFF added to miniprot output**; 5.8.0/5.8.1 odb12 compatibility; 5.8.3 autolineage offline fix; 6.0.0 miniprot post-processing, gzipped input, "Allow odb10 and odb12 specific dataset downloads"; 6.1.0 odb12.2 datasets "while retaining compatibility with legacy odb12 and odb10", updated HMMER scoring | confirmed | changelog text (no dates in it) | GL CHANGELOG |
| U3 | Python requirement fits `>=3.10,<3.12` | partly | `setup.py` (5.8.3 and master) has no `python_requires`, uses `distutils` (removed in Python 3.12); bioconda pins only `python >=3.3`. So 3.10/3.11 are fine; 3.12+ from source is doubtful | GL setup.py; BC |
| U4 | bioconda `busco` 6.1.0 run deps: python>=3.3, numpy<2.4, wget, blast>=2.10.1, hmmer>=3.1b2, augustus>=3.3, biopython>=1.79, prodigal, **sepp==4.5.5**, metaeuk>=6.a5d39d9, pandas, bbmap, miniprot, requests, matplotlib-base; 5.8.3: blast>=2.16, augustus>=3.5.0, sepp unpinned, plus r-base/r-ggplot2 | confirmed | meta.yaml | BC |
| U5 | Fit with pixi pins (`augustus>=3.5.0`, `blast>=2.16`, `hmmer>=3.4`, `miniprot 0.18`, `biopython>=1.87`, `pandas`) | partly | no nominal conflict; `numpy<2.4` and `sepp==4.5.5` not checked against the lock; **not solved** (busco absent from pixi.lock, L9) | BC; pixi.toml |
| U6 | Platforms: `busco` is noarch; augustus 3.5.0, metaeuk 7.bba0d80, sepp 4.5.6, miniprot 0.18, blast 2.17.0, bbmap 40.05 have osx-arm64 builds; sepp **4.5.5** on osx-arm64 not checked | partly | files pages | AC |
| U7 | v6.1.0 still supports `--augustus`, `--augustus_species`, `--augustus_parameters`, `--long` in genome mode | confirmed | argparse in run_BUSCO.py; UG "To use the Augustus pipeline add the --augustus option"; UG also: "Augustus gene predictions are the most reliable for downstream analysis" | GL run_BUSCO.py; UG |
| U8 | Augustus mode still writes `augustus_output/{extracted_proteins, gb, gff, predicted_genes_initial_run, predicted_genes_rerun, retraining_parameters/BUSCO_<out>, training_set.db}` | confirmed | UG tree; `augustus.py move_retraining_parameters()` | UG; GL augustus.py |
| U9 | `$AUGUSTUS_CONFIG_PATH/species` must be writable in Augustus mode | confirmed (source) | `os.access(..., os.W_OK)` else "Cannot write to Augustus species folder"; the guide does not say so | GL augustus.py |
| U10 | Gene finders: metaeuk default v5.0-5.6; miniprot default since 5.7.0; `--metaeuk` and `--augustus` opt-in in 6.1.0; UG: miniprot "faster than Augustus and Metaeuk (except on fungi, in which Metaeuk is still faster)" | confirmed | UG; CHANGELOG | UG; GL |
| U11 | No genetic-code / translation-table flag for eukaryotic genome mode; miniprot call hard-coded (`--trans -u -I --outs 0.95 -t N --gff`, no `-T`); no `--miniprot_parameters` | confirmed | run_BUSCO.py argparse; miniprot.py | GL |
| U12 | `--augustus_parameters` is a whitelist that includes `translation_table`; `--metaeuk_parameters` whitelist includes `translation-table`; format `--k=v,--k=v`; other keys dropped with a warning | confirmed | `ACCEPTED_PARAMETERS` in augustus.py / metaeuk.py; UG advises using it "to fix biologically relevant parameters such as the translation table" | GL |
| U13 | Does upstream Augustus mode pass the table to its tblastn step (`-db_gencode`)? | unconfirmed | not checked in the source; the vendored script does (script :1094-1096) | — |
| U14 | miniprot itself has `-T INT` (NCBI tables 1-33 except 7-8, 17-20); MetaEuk/MMseqs2 has `--translation-table` | confirmed / partly | miniprot man page and `main.c`; MMseqs2 `Parameters.cpp` (MetaEuk easy-predict honouring it not verified) | lh3/miniprot; soedinglab/MMseqs2 |
| U15 | `full_table.tsv` genome-mode header: `# Busco id, Status, Sequence, Gene Start, Gene End, Strand, Score, Length, OrthoDB url, Description`; proteins/transcriptome: `# Busco id, Status, Sequence, Score, Length [, OrthoDB url, Description]`; two leading comment lines (version, dataset); `Missing` rows have 2 columns | confirmed | hmmer.py writer; real v6.1.0 file | GL hmmer.py; IUC |
| U16 | Layout: `<out>/run_<lineage>/{busco_sequences/{single_copy,multi_copy,fragmented}_busco_sequences/*.faa,*.gff, full_table.tsv, hmmer_output/, miniprot_output/ or metaeuk_output/ or augustus_output/, missing_busco_list.tsv, short_summary.{txt,json}}`, top-level `short_summary.specific.<dataset>.<out>.txt`, `logs/`; Augustus mode adds `blast_db/`, `blast_output/` | confirmed | UG | UG |
| U17 | Per-gene GFF and protein per BUSCO in `busco_sequences/` (miniprot: no nucleotide FASTA, "can frequently contain internal stop codons"; metaeuk writes `*.fna` too); GFF keeps original coordinates, full_table has trimmed ones | confirmed | UG | UG |
| U18 | v5 vs v6 layout change | partly | no documented change; GFF added to miniprot output in 5.7.0; `<id>:<start>-<stop>` coordinates since 5.4.7; a real v5-vs-v6 diff was not made | GL CHANGELOG |
| U19 | Download: `--download <name|all|eukaryota|...>`, `--download_path` (default `busco_downloads/`), `--offline`, `--list-datasets [odb10|odb12|odb12.2]`, `--auto-lineage`, **`--datasets_version` default odb12.2**, `--download_base_url`; `all_odb10` etc. for legacy sets | confirmed | run_BUSCO.py; UG (UG still mentions the removed `--update-data`) | GL; UG |
| U20 | Offline: with `--offline` `file_versions.tsv` is not fetched ("Unable to verify BUSCO datasets because of offline mode"); dataset must exist as `<download_path>/lineages/<name>` (no date suffix); without `--offline` it fetches the index (retries 10 s/100 s, backup busco-data2.ezlab.org) and auto-updates | confirmed | BuscoDownloadManager.py | GL |
| U21 | Tarball contents odb10 = L5; odb12/odb12.2 same minus `.md5`, `links_to_ODB12[.2].txt`; **`lengths_cutoff` removed for odb12 and "reintroduced" for odb12.2** (UG); eukaryota_odb12 has none, microsporidia_odb12.2 has one | confirmed | extracted tarballs; UG | LIN; UG |
| U22 | file_versions.tsv: 193 odb10 (2024-01-08), 497 odb12 (2026-05-22), 500 odb12.2 (2026-05-13) lineages; fungal sets in all three: fungi, ascomycota, basidiomycota, eurotiomycetes, sordariomycetes, saccharomycetes, dothideomycetes, leotiomycetes, agaricomycetes, tremellomycetes, mucorales, mucoromycota, microsporidia (+ order/family sets); chytridiomycota odb12+ only; capnodiales odb10 only; **no dikarya, pezizomycotina, glomeromycota in any version** | confirmed | FV | FV |
| U23 | Marker counts odb10: eukaryota 255, fungi 758, ascomycota 1706, basidiomycota 1764, mucoromycota 1614, microsporidia 600, saccharomycetes 2137, eurotiomycetes 3546, sordariomycetes 3817, dothideomycetes 3786, leotiomycetes 3234, agaricomycetes 2898, tremellomycetes 4284, mucorales 2449 | confirmed | dataset.cfg (7 sets) + https://busco.ezlab.org/list_of_lineages.html (rest) | LIN; UG list |
| U24 | Marker counts odb12 (OrthoDB 12.1): eukaryota 129, fungi 1122, ascomycota 2826, basidiomycota 2409, saccharomycetes 2319, eurotiomycetes 3665, sordariomycetes 4492, dothideomycetes 4643, microsporidia 425, chytridiomycota 1017, mucoromycota 1560; odb12.2: fungi 1019, ascomycota 2557, basidiomycota 1811, eukaryota 125, microsporidia 514 (others not scanned) | confirmed / partly | dataset.cfg | LIN |
| U25 | Versus odb9 (resources.py:197): fungi 290 -> 758/1122; ascomycota 1315 -> 1706/2826; basidiomycota 1335 -> 1764/2409; microsporidia 518 -> 600/425; eukaryota 303 -> 255/129; eurotiomycetes 4046 -> 3546/3665; dikarya 1312 -> none | confirmed | U23/U24 vs resources.py | — |
| U26 | Id scheme: odb10 and odb12 use `<number>at<taxid>` (`100957at4751`, `364at4751`, `261328at2759`); odb9 `EOG09...` appears in no current file | confirmed | scores_cutoff / links files; real full_table (`774at6656`) | LIN; IUC |
| U27 | A published odb9 -> odb10/odb12 id mapping | unconfirmed (none found) | OrthoDB v10 paper (PMC6323947) says nothing about OG id continuity; `links_to_ODB10.txt` only maps id -> description/url | PMC6323947 |
| U28 | Sizes (MB): odb10 fungi 126, ascomycota 204, basidiomycota 134, eukaryota 66, eurotiomycetes 192, sordariomycetes 201, saccharomycetes 107, microsporidia 26; odb12 fungi 120, ascomycota 296, eukaryota 28; **odb12.2 fungi 949, ascomycota 1592, basidiomycota 442, eukaryota 231** (odb12.2 ~5-8x odb12) | confirmed | Content-Length | LIN |
| U29 | macOS arm64 | partly | UG: "We do not support macOS but BUSCO should work on it, although Augustus seems to cause troubles on some BSD-derived systems, including macOS"; arm64 never mentioned; installable from bioconda in principle (U6); not tried | UG; AC |
| U30 | Docker `ezlabgva/busco:v6.1.0_cv3` ~859 MB compressed (v5.8.2_cv1 1.05 GB) | confirmed | Docker Hub tags | https://hub.docker.com/r/ezlabgva/busco/tags |
| U31 | Papers: v5 Manni et al. 2021 MBE 38:4647 (doi:10.1093/molbev/msab199); v6/odb12 citation per UG is Tegenfeldt et al. 2025 NAR 53:D516 (doi:10.1093/nar/gkae987); no separate v6 paper found | confirmed / partly | UG; journal pages | OUP |
| U32 | "BUSCO v5.8.0 made miniprot the default" | **corrected**: 5.7.0 | CHANGELOG | GL |
| U33 | docs/BUSCO_REVIEW.md claim "v6 drops Augustus and uses odb12, so the Augustus-trained-by-BUSCO workflow disappears" | **wrong** | U7, U8: Augustus mode exists in 6.1.0 with retraining output; odb10 datasets still downloadable and supported (U2, U19) | GL; UG |
| U34 | Upstream per-job timeouts for Augustus/miniprot jobs | unconfirmed | not searched in the source; nothing in the user guide | — |

### 1.2 Rows verified locally (section 7 has the commands)

| # | Claim | Status | Evidence |
|---|---|---|---|
| L1 | `https://busco-data.ezlab.org/v5/data/` now redirects (301) to `https://busco-data.s3.amazonaws.com/`; `file_versions.tsv` has 1,248 rows | confirmed | curl -L; `setupDB.download_busco_v5` (setupDB.py:102-131) does not follow redirects by itself and never downloads a lineage (it builds URLs and stops; dead code) |
| L2 | No `dikarya`, `pezizomycotina`, `saccharomycetales`, `sordariomyceta`, `saccharomyceta`, `glomeromycota` dataset exists in odb10 or odb12 | confirmed | `file_versions.tsv` rows for fungi: `fungi`, `ascomycota`, `basidiomycota`, `agaricomycetes`, `chytridiomycota` (odb12 only), `dothideomycetes`, `eurotiomycetes`, `leotiomycetes`, `microsporidia`, `mucorales`, `mucoromycota`, `saccharomycetes`, `sordariomycetes`, `tremellomycetes`, `eukaryota`; each in `_odb10` (2024-01-08), `_odb12` (2026-05-22) and `_odb12.2` (2026-05-13) |
| L3 | Tarball sizes: fungi_odb10 126 MB, fungi_odb12 120 MB, fungi_odb12.2 949 MB, eukaryota_odb10 66 MB, ascomycota_odb10 204 MB, ascomycota_odb12 296 MB, eurotiomycetes_odb10 192 MB, microsporidia_odb10 26 MB, mucoromycota_odb10 74 MB; unpacked fungi_odb10 278 MB, fungi_odb12 337 MB (odb9 dikarya: 63 MB tarball, 301 MB unpacked) | confirmed | `curl -IL` Content-Length; `du -sh` |
| L4 | Marker counts: fungi_odb10 758 (549 species), fungi_odb12 1,122 (144 species); odb9 fungi 290, dikarya 1,312 (`resources.py:197`, `$FUNANNOTATE_DB/dikarya/dataset.cfg`) | confirmed | `dataset.cfg` of each downloaded dataset |
| L5 | odb10 dataset layout = odb9 layout as read by the vendored script: `hmms/`, `prfl/`, `scores_cutoff`, `lengths_cutoff`, `ancestral`, `ancestral_variants`, `dataset.cfg` (+ `refseq_db.faa.gz`, `links_to_ODB10.txt`); HMMER3/f format, `NAME` = file stem (`100957at4751`) | confirmed | `ls`, `head`; the vendored script's checks are `check_dataset` :480-493 (hmms), :2184-2210 (prfl, ancestral_variants), `_load_score` :1940-1950, `_load_length` :1963-1973 |
| L6 | Vendored script, proteins mode, `-l fungi_odb10` on `fungi_odb10/ancestral`: exit 0, `C:94.3%[S:94.3%,D:0.0%],F:1.3%,M:4.4%,n:758`, 6.3 s with 8 cpus, HMMER 3.4 | confirmed | section 7, run 1 |
| L7 | Vendored script on `fungi_odb12`: fails `Impossible to read the lengths in .../lengths_cutoff` (file absent). Second, independent blocker: odb12 HMM `NAME` is `100036` while the file is `100036at4751.hmm` and `scores_cutoff` keys are `364at4751`; the parser keys hits by the hmmsearch query name (:1733 `busco_query = line[3]`), so even with a synthetic `lengths_cutoff` every hit would miss the cutoff dictionary | confirmed | section 7, run 2; `head` of the files |
| L8 | Augustus 3.5.0 (pixi env, GCC build) accepts odb10 and odb12 `.prfl` files with `--proteinprofile` (warnings about insignificant blocks only; same warnings for odb9) | confirmed | section 7, run 3 |
| L9 | The pixi `default` env has `augustus 3.5.0`, `hmmer 3.4`, `miniprot 0.18`, BLAST+ tblastn, Python 3.11.14; it has **no** `busco` and **no** `metaeuk`; `pixi.lock` has no `busco` entry | confirmed | `ls .pixi/envs/default/bin`, `grep -c busco pixi.lock` = 0; `pixi.toml:40-64, 153, 179` |
| L10 | HMMER 3.4 reads the odb9 HMMs (`HMMER3/f [3.1b1]`) and the odb10 HMMs (`HMMER3/f [3.1b2]`) | confirmed | `hmmstat`, run 1 |

## 2. Use-map of BUSCO in 1.10 (file:line)

Classification: **(a)** predict training-set generation, **(b)** annotate functional annotation,
**(c)** setup/database, **(d)** other (CLI, docs, tests, CI, external package).

### 2.1 predict (a)

| Where | What |
|---|---|
| `predict.py:222-228` | `--busco_seed_species` (default `anidulans`) |
| `predict.py:238` | `--busco_db` (default `dikarya`) |
| `predict.py:240-244` | `--augustus_BUSCO_timeout` (default 300 s) |
| `predict.py:301-302, 312-323` | `--min_training_models`, `--busco_fallback/--no_busco_fallback` |
| `predict.py:465-470` | output-folder name conflict list contains `busco`, `busco_proteins` |
| `predict.py:494-499` | lineage dir `$FUNANNOTATE_DB/<busco_db>` must exist (hint `setup -i busco -b`) |
| `predict.py:795-799` | seed species copied from `$FUNANNOTATE_DB/trained_species` into the run-local Augustus config (`LOCALAUGUSTUS`) |
| `predict.py:975-991, 1035-1036, 1060-1061` | `RunModes[*] = "busco"` when no PASA/pretrained; seed species validated; source text has the typo `BUCSCO` (3 places) |
| `predict.py:1123-1142` | `RunBusco` flag; `training_mode_initial` rows in `training_decisions.tsv` |
| `predict.py:1154-1166` | `augustus --proteinprofile` functional gate (`lib.checkAugustusFunc`, library.py:1131-1175, uses `config/busco_test.fa`) — BUSCO training refused when it fails |
| `predict.py:2051-2252` | `run_busco_training()`: script choice :2055-2062; lineage path :2063; cache key `run_<sp>/full_table_<sp>.tsv` :2065-2069, 2083-2088; seed fallback to `anidulans` with warning :2070-2080; tblastn thread note :2093-2099; command :2106-2126 (`-m genome --lineage <dir> --species <seed> -f --local_augustus <dir>` + `lib.busco_translation_args(args.table)` + `--augustus_timeout`); exit code checked :2136-2142; table checked :2145-2149; `parseBUSCO2genome` :2152-2155; per-marker `augustus_output/gffs/<id>.gff` converted with EVM's `augustus_GTF_to_EVM_GFF3.pl` (`Converter2`, :1401) :2164-2179; `fix_busco_naming` :2187-2194; **second run in proteins mode** on the predicted proteins :2196-2234; `getCompleteBuscos` + `filterGFF3` :2235-2240; decision row `busco_training` :2243-2245 |
| `predict.py:2264-2306` | PASA gate / identity gate -> switch to BUSCO (`lib.resolve_pasa_gate`) |
| `predict.py:2390-2434` | `--min_training_models` fallback to BUSCO (`lib.resolve_min_training_models`) |
| `predict.py:2452-2464` | `gff2gbSmallDNA.pl <FinalTrainingModels> genome 600` + `lib.trainAugustus` |
| `predict.py:2754-2780, 2817-2843` | SNAP and GlimmerHMM train from the same `FinalTrainingModels` (= `busco.final.gff3` in BUSCO mode) |
| `library.py:958-966` | `tblastn_single_thread` (string compare; see tests/test_tblastn_threads.py) |
| `library.py:969-979` | `busco_translation_args(table)`: `-a --translation_table=N --db_gencode N` |
| `library.py:1492-1540` | `fix_busco_naming`: 7-column table (:1505-1508), replaces `Augustus%20prediction` (:1532) |
| `library.py:6972-7008` | `dupBUSCO2gff`: reads `hmmer_output/` and `augustus_output/predicted_genes/` of the vendored layout (ploidy>1 only) |
| `library.py:7011-7026` | `getCompleteBuscos`: unpacks exactly 5 columns (:7023) — proteins-mode table |
| `library.py:7040-7090` | `parseBUSCO2genome`: 7-column genome table (`id status contig start end score length`, :7052-7056), best-scoring duplicate for ploidy>1 (:7061-7077), bed +/-100 bp |
| `library.py:10542-10590` | pure decision functions (`pasa_gate_applies`, `resolve_pasa_gate`, `resolve_min_training_models`) |
| `library.py:1442-1489` | `genemark2busco`, `evidence2busco`: **no callers** (dead) |
| `aux_scripts/funannotate-BUSCO2.py` | genome mode: tblastn :1094-1096 (`-db_gencode` added 2026-10-05, commit 16032e0); Augustus per marker with `prfl/<id>.prfl` :1255, 1264; retraining `new_species.pl`/`etraining`/`optimize_augustus.pl` :1402-1418; rerun of missing/fragmented :1478-1486; hmmsearch :1505; parser :1707-1770; `full_table` writer :1872-1920 (7 columns genome, 5 columns proteins); per-job timeout + reap guard :343-359; worker-queue locks :1141-1170, 2105-2170; `_load_length` :1963-1973; arg parser :3031-3125 (`--augustus_timeout` :3082, `-a` :3087, `--db_gencode` :3093, `--long` :3106); restart command keeps `--db_gencode` :3365-3366; `dataset.cfg` parse :3171-3195 |
| `aux_scripts/funannotate-BUSCO2-py2.py` | selected only when `sys.version_info <= (3, 0)` (predict.py:2055-2062, library.py:6939-6942); **does not accept** `--augustus_timeout` or `--db_gencode` (0 matches) and lacks the deadlock fixes; `setup.py:21` requires Python `>=3.10.0, <3.12`, so it is unreachable |

Vendored-script history relevant to (a) (`git log -- funannotate/aux_scripts/funannotate-BUSCO2.py`):
33558df 2026-05-01 `-a` alias; 1387721 05-18 py2->3 checks; 61fe434 06-06 per-job Augustus timeout and
timing profile (PR #1173, CHANGES.md:375-397); eb8da6c 06-07 timeout-safe rerun extraction; f1e9c73
06-26 two deadlocks fixed; 79b7f0f 08-06 crash/silent-failure fixes; 16032e0 10-05 `-db_gencode`;
c3741f9 10-24 fasta-pearson parsing + escape fixes. The py2 copy received only 1387721 and c3741f9.

### 2.2 annotate (b)

| Where | What |
|---|---|
| `annotate.py:592` | `--busco_db` (default `dikarya`) |
| `annotate.py:746-751` | lineage dir must exist |
| `annotate.py:1593-1608` | `lib.runBUSCO(Proteins, $FUNANNOTATE_DB/<busco_db>, cpus, annotate_misc, annotations.busco.txt)`; cached by file presence |
| `library.py:6937-6969` | `runBUSCO`: runs the vendored script **by shebang** (`cmd[0]` is the script path, no `sys.executable`, :6943-6944) in proteins mode, reads `run_busco/full_table_busco.tsv` (:6961), writes `<protein>\tnote\tBUSCO:<id>` for `Complete` and `Duplicated` (:6968-6969); `runSubprocess` exits on non-zero (library.py:772-784) but a zero-hit run is silent |
| `annotate.py:1824-1845` | `annotations.busco.txt` merged into `all.annotations.txt` -> GenBank `note` qualifiers |
| `library.py:2728, 2858-2861` | `gff2dict` stats count `BUSCO:` notes |
| `library.py:9641, 9674, 9724-9726, 9842, 9945-9947` | annotation table writers: `BUSCO` column from `BUSCO:` notes |
| `library.py:6871-6934` | `runBUSCOv5`: **no callers**; its `genome` branch passes `-m proteins` (:6905-6907); hard-codes `run_<lineage>_odb10/full_table.tsv` (:6921-6923); not a working implementation |

### 2.3 setup / database (c)

| Where | What |
|---|---|
| `funannotate/downloads.json` `"busco"` | 28 names -> `https://github.com/nextgenusfs/funannotate/releases/download/db-mirror/<name>_odb9.tar.gz` and the folder to rename (note the mismatches `sordariomycetes` -> `sordariomyceta_odb9.tar-2.gz`, `saccharomycetes` -> `saccharomyceta_odb9`, `protists`/`alveolata_stramenophiles` -> `_ensembl`) |
| `setupDB.py:496-528` | `download_buscos`: download, `tar -zxf`, rename `<folder>` -> `$FUNANNOTATE_DB/<name>` |
| `setupDB.py:592-595, 605-611, 651-661, 677-679, 723-724` | `-i busco`, `-b` choices (28 fixed names), `busco_outgroups` refused, `busco_links` from the JSON (GitHub master copy first, :664-674) |
| `setupDB.py:102-131` | `download_busco_v5`: unfinished, no callers |
| `resources.py:197` | `buscoTree` with odb9 counts, printed by `database.py:23-27` (`--show-buscos`) |
| `scripts/mirror_db_to_github.py` | mirrors the 28 odb9 tarballs to the `db-mirror` release |
| `Docker/README.md:67-71`, `.github/workflows/image-prune-test.yml:94-109`, `image-test.yml:28` | image notes ("all BUSCO lineage databases", 10-15 GB); CI runs `setup -i repeats busco -b dikarya` under a 900 s timeout |
| `plans/IMAGE_DISTRIBUTION_PLAN.md:59-60` | slim image needs `setup -i busco -b dikarya` before `test -t predict` |

### 2.4 Other (d)

| Where | What |
|---|---|
| `funannotate.py:186-193, 310, 363, 429, 445-450` | help text; `test -t busco`; `setup -b`; removed-command message for `outgroups` |
| `test.py:74-86, 224-288, 408, 426` | `runBuscoTest`: downloads `test-busco.tar.gz` (OSF), runs `predict` with protein evidence only (BUSCO-trained), asserts 1,500-1,800 genes, then re-runs with `-p` |
| `install_scripts/install_augustus_gcc_macos.sh:7,116`, `macos_smoke_checks.sh:14` | macOS: clang Augustus crashes in `--proteinprofile`; GCC build needed for BUSCO training |
| `docs/predict.rst:56-61` | stale pipeline description ("Combined GeneMark and BUSCO2 results, feed into Evidence Modeler") — the code no longer runs EVM on BUSCO models (predict.py:2157 comment vs :2164-2240) |
| `docs/training.rst`, `docs/training_data_selection_methods.md`, `docs/commands.rst:359-365, 607, 644, 675`, `docs/tutorials.rst:144-146, 196-202` (`--show_buscos` should be `--show-buscos`), `docs/annotate.rst:85`, `docs/databases.rst:48`, `docs/conda.rst:404-409`, `docs/containers.rst:51` | user docs |
| `tests/test_busco_db_gencode.py` | `busco_translation_args`; source checks that predict and the script wire `--db_gencode`; real tblastn `-db_gencode 4` test |
| `tests/test_tblastn_threads.py` | thread gate + source sync with the script |
| `tests/test_training_source_decision.py` (31), `test_training_gates.py`, `test_training_selection.py` | PASA-vs-BUSCO routing, pure functions |
| `tests/test_augustus_base.py`, `test_test_data_download.py`, `test_release_version.py` | incidental mentions |
| `tests/test_removed_commands.py:40-43` | `test_setup_busco_outgroups_refused` **fails in isolation** (`mock.patch.object(setupDB.lib, "log")`: `lib.log` exists only after `setupLogging`); passes only when another test set the logger first |
| **No** test exercises `parseBUSCO2genome`, `getCompleteBuscos`, `fix_busco_naming`, `dupBUSCO2gff`, `runBUSCO` or the script itself; the only end-to-end coverage is `funannotate test -t busco` (network, Augustus, ~tens of minutes) |
| External: `funannotate-compare` (`/Users/jstajich/projects/funannotate-compare`, HEAD 3f4cf0e) | `outgroups.py:9-20, 37-60`: finds `aux_scripts/funannotate-BUSCO2.py` inside the installed funannotate (`--busco-script`, `$FUNANNOTATE_BUSCO2`); `:154-166` runs it in proteins mode on `$FUNANNOTATE_DB/<lineage>`; `lib.py:1270-1271` reads `BUSCO:` notes from GenBank; `compare.py:1140-1165, 1198-1220` builds the single-copy-ortholog phylogeny from BUSCO ids shared across genomes and the outgroup `<species>.<lineage>_buscos.fa` (sequence names = `EOG...` ids, `tests/test_outgroups.py:13-20`) |

## 3. What can and cannot be replaced by upstream BUSCO

### 3.1 (b) annotate — replaceable

What the vendored proteins mode does: hmmsearch of every lineage HMM against the proteome, score and
length cutoffs, `Complete/Duplicated/Fragmented/Missing` (script :2930-2990, :1707-1920). Upstream
BUSCO proteins mode does the same with the same cutoff files; the only coupling in funannotate is
`library.runBUSCO` (:6937-6969), which reads two columns (`id`, `status`, `protein`) of a 5-column
table. A 10-line reader for the upstream proteins-mode `full_table.tsv` (U15: `# Busco id, Status,
Sequence, Score, Length[, OrthoDB url, Description]` after two comment lines; same first three columns,
so `runBUSCO`'s `col[0..2]` logic already fits) replaces it.
What changes for users is **the id namespace** (section 5.1), not the mechanics.

Cheaper than upstream for 1.10: keep the vendored proteins mode and allow odb10 lineage directories
(section 1.2 L5-L6). No new dependency, no Python constraint, works offline, same code path for
`funannotate-compare outgroups`.

### 3.2 (a) predict — not replaceable in 1.10

Points where an external run cannot be dropped in, with evidence:

1. **Output coupling.** `parseBUSCO2genome` expects 7 columns with score in column 6
   (library.py:7052-7056); `fix_busco_naming` keys models by `contig:start-end` from the same table
   and rewrites `Augustus%20prediction` (:1505-1508, 1532); per-marker GFF files are read from
   `run_<sp>/augustus_output/gffs/<id>.gff` and converted with an EVM script (predict.py:2164-2179);
   ploidy>1 reads `hmmer_output/` and `augustus_output/predicted_genes/` (library.py:6972-7008).
2. **Second pass.** predict re-runs the script in proteins mode on the predicted proteins and keeps
   only models still `Complete` (predict.py:2196-2240). With an external run this pass needs the same
   lineage in the same format.
3. **Translation table.** `--table` flows to Augustus (`-a --translation_table=N`) and to tblastn
   (`--db_gencode N`) through `lib.busco_translation_args` (library.py:969-979; script :1094-1096,
   3093-3095, 3365-3366), with a real tblastn regression test (tests/test_busco_db_gencode.py). Upstream
   BUSCO has no genetic-code option in eukaryotic genome mode; the default miniprot job is hard-coded
   without `-T` (U11); only the opt-in Augustus/MetaEuk modes can receive a table through the
   `--augustus_parameters`/`--metaeuk_parameters` whitelists (U12), and whether upstream's Augustus
   mode passes it to its tblastn step is unconfirmed (U13). The miniprot/metaeuk modes also change the
   evidence type from ab initio Augustus models to protein-to-genome alignments (miniprot proteins
   "can frequently contain internal stop codons", U17), which the methods draft has not evaluated as
   training data.
4. **Run-local Augustus config.** predict copies the seed species into
   `predict_misc/ab_initio_parameters/augustus` and passes `--local_augustus` (predict.py:795-799,
   2122-2123); the script trains `BUSCO_<abrev><rand>` there. Upstream `--augustus` writes into
   `$AUGUSTUS_CONFIG_PATH/species`, which must be writable and is shared between concurrent runs (the
   BFD production set is 21,762 genomes, docs/training_data_selection_methods.md:137).
5. **Timeouts and deadlocks.** `--augustus_BUSCO_timeout` with SIGKILL of the process group, timing
   summary and the two queue deadlock fixes (script :343-359, 1141-1170, 2105-2170; CHANGES.md:375-397)
   exist only in the vendored script.
6. **Seed species.** `--busco_seed_species` with validation (predict.py:982-990) and a logged fallback
   to `anidulans` (predict.py:2070-2080); upstream `--augustus_species` would need the same mapping, and
   the metaeuk/miniprot modes have no seed at all.
7. **Decision log and tests.** `busco_training`, `training_mode_switch`, `min_training_models` and
   `final_training_source` rows (predict.py:2243-2245, 2302-2305, 2392-2397, 2441-2447) and 31+ tests
   assume `busco.final.gff3` is produced by this path with these counts.
8. **Benchmarked behaviour.** The 1,000-complete-model PASA gate and the BUSCO fallback were tuned
   against BUSCO training *as implemented* (docs/training_data_selection_methods.md sections 4-6,
   Tables M2, M5, M6; funannotate-manuscript experiments A-D). Changing the training-set generator
   invalidates those thresholds until re-measured.

What **is** feasible for (a) in 1.10 without touching the engine: odb10 lineage data (more markers:
fungi 758 vs 290; the dataset `species=aspergillus_nidulans` matches the default seed), because the
genome mode reads the same files (`prfl/`, `ancestral_variants`, cutoffs; L5, L8). This needs the
genome-mode validation run listed in section 4.3 (not run here: no test genome on this machine).

## 4. Feasibility, recommendation, phased plan

### 4.1 Summary

| Item | 1.10 | 1.11 | Not planned |
|---|---|---|---|
| Hygiene (section 6) | yes | | |
| odb10 lineage datasets usable by the vendored script (`--busco_db fungi_odb10`, `setup -i busco -b fungi_odb10`) | yes (medium) | | |
| Id-namespace tag so odb9 and odb10 annotations cannot be mixed silently | yes (small) | | |
| `annotate --busco_results <upstream full_table.tsv>` (user-supplied, any BUSCO version) | optional (small) | yes | |
| `annotate` shelling out to upstream `busco` when installed | | yes (medium) | |
| `predict --busco_results <run dir>` import through a normalizer | | yes, after benchmark (large) | |
| `predict` running upstream `busco --augustus`/`--miniprot` itself | | | until benchmark passes and v6 is pinned in pixi |
| odb12 support | | only via upstream BUSCO (L7) | in the vendored script (odb12.2 restores `lengths_cutoff`, U21, but the HMM `NAME` suffix issue of L7 was not tested on odb12.2) |
| Delete `funannotate-BUSCO2.py` | | | blocked by `funannotate-compare outgroups` and by (a) |
| Delete `funannotate-BUSCO2-py2.py` | yes | | |
| Stop hosting odb9 tarballs | | | keep for id continuity; add odb10 pointers instead |

### 4.2 Phase 0 (1.10): hygiene — section 6, all small, no behaviour change.

### 4.3 Phase 1 (1.10): odb10 lineages for the existing engine

Work items:

1. **(small)** `setupDB.download_buscos`: accept any `<lineage>_odb10` name; resolve the URL from
   `file_versions.tsv` (follow redirects: L1), verify md5 from the same file, extract to
   `$FUNANNOTATE_DB/<name>`; keep the 28 odb9 names and the GitHub mirror. Replace the fixed `-b`
   `choices` list with validation against both lists. Write `dataset.cfg` provenance to
   `funannotate-db-info.txt`.
2. **(small)** `resources.buscoTree` / `database --show-buscos`: print the odb9 tree and the odb10
   fungal lineages with counts read from `dataset.cfg` when installed.
3. **(small)** predict/annotate: accept `--busco_db fungi_odb10`; the only code change is the
   directory-exists check (predict.py:494-499, annotate.py:746-751) plus a one-line log of the dataset
   name and `creation_date` into the run log and `training_decisions.tsv`.
4. **(small)** Id-namespace tag: write the dataset name into the annotation note
   (`BUSCO:100957at4751` is self-describing, `EOG...` is odb9; add a `funannotate-db-info` row and a
   header comment in `annotations.busco.txt`) and make the annotation-table writer keep the full id.
   Document the namespace change for `funannotate-compare`.
5. **(medium)** Validation of genome mode on odb10: run `funannotate test -t busco` with
   `--busco_db fungi_odb10` and the standard `dikarya` back to back on one machine; compare
   `busco.final.gff3` counts, `training_decisions.tsv` and the gene-count assertion (1,500-1,800).
   Also one `--ploidy 2` run (the `dupBUSCO2gff` path) and one `--table 4` run (tblastn + Augustus
   with the table).
6. **(small)** `test.py`: `-t busco` takes `--busco_db` so the CI image test can run the odb10 variant.
7. **(small)** docs: `databases.rst`, `conda.rst:404-409`, `commands.rst`, `tutorials.rst:196-202`.

Acceptance tests for Phase 1:
- unit: `download_buscos` resolves names/urls/md5 from a fixture `file_versions.tsv`; rejects odb12
  names with a message that names the reason (L7); the existing odb9 path unchanged.
- unit: the proteins-mode self-test of L6 reduced to a fixture (10 HMMs + 10 ancestral sequences cut
  from fungi_odb10, < 1 MB), run through `lib.runBUSCO`, expecting 10 `Complete` and the exact note
  lines; the same fixture for `getCompleteBuscos` and `parseBUSCO2genome` (synthetic 7-column table).
- integration (CI image, network): `funannotate test -t busco` with `dikarya` and with `fungi_odb10`;
  `funannotate test -t annotate` with `fungi_odb10`.
- regression: `tests/test_busco_db_gencode.py`, `test_tblastn_threads.py`, the training-decision suite
  (99 tests pass today in 0.9 s; one unrelated failure, section 6 item 9).

### 4.4 Phase 2 (1.11): external BUSCO, import first

1. **(small)** `annotate --busco_results FILE`: accept an upstream `full_table.tsv` (proteins mode) or
   the vendored 5-column table; detect by header/column count; write the same note lines.
2. **(medium)** `annotate` runs `busco -m proteins -l <lineage> --offline --download_path <dir>` when
   `busco` is on PATH and `--busco_db` names an odb10/odb12 lineage that is not installed in the
   vendored layout; version-gate on `busco --version` (>= 5.8.0 for odb12 tables, 6.1.0 for odb12.2;
   U2); pin in `pixi.toml` only if the solver accepts `busco` (noarch; `numpy<2.4`, `sepp==4.5.5`,
   `metaeuk`, `bbmap`, `prodigal`, `matplotlib-base`; U4) next to the current `augustus>=3.5`,
   `blast>=2.16`, `hmmer>=3.4`, `miniprot 0.18` pins on both `linux-64` and `osx-arm64` (U5, U6;
   unresolved today because `busco` is absent from the lock, L9). Always run it with
   `--offline --download_path` and a dataset directory named exactly `<name>` (U20); never rely on
   `-l fungi` without a version suffix, which now pulls odb12.2 (U19, U28).
3. **(large)** `predict --busco_results DIR|TSV|GFF`: the normalizer sketched in
   `docs/BUSCO_REVIEW.md` section 5 (detect layout, parse 7/10-column tables, collect per-marker GFF
   from `augustus_output/gffs`, `miniprot_output` or `metaeuk_output`, filter `Complete` (+ best
   `Duplicated` when ploidy>1), validate CDS under `--table`, write `busco.final.gff3`, skip the
   proteins-mode re-validation or run it on the same lineage). Gate on the benchmark below.
4. **(medium)** Benchmark (reuses the manuscript structure: `funannotate-manuscript/
   assessment_pasa2.6_fun1.9/scripts/expD_arm_opt.sh` = arm runner with step A on training
   chromosomes / step B whole genome with `-p`, `predict_scorer.py split|score` for holdout
   gffcompare, proteome BUSCO 5.8.0 fungi_odb10; BUSCO-forced arm = `--min_pasa_complete_models
   1000000000`). New arms, same genomes as experiment D (16 RefSeq genomes) or the 5 of Table M2:
   - `busco_odb9` (today), `busco_odb10` (Phase 1), `ext_aug_odb10` (upstream `--augustus
     --augustus_species anidulans`, imported), `ext_miniprot_odb10` and, if available,
     `ext_odb12` (imported), each 3 repeats where the arm is stochastic;
   - metrics: training models written / surviving CDS validation; holdout locus, exon, intron-chain
     Sn/Pr (`predict_scorer.py`); proteome BUSCO C %; BUSCO-stage wall time and peak RSS at 16 cpus;
     for one table-12 genome (`Candida albicans`) the count of models with internal stops under
     `--table 12`;
   - decision rule as in the methods draft: adopt an external arm as a *default* only if its holdout
     locus F1 is within the experiment-D noise band (bootstrap 95 % interval of the difference
     includes 0, no genome worse by > 1 point) and it is faster; otherwise it stays opt-in.

### 4.5 Phase 3 (later): decide the default; retire the vendored genome mode only when
`funannotate-compare` no longer imports the script (it needs a proteins-mode engine of its own or an
upstream call), the benchmark passed, and a BUSCO version with the Augustus mode is pinned in pixi for
both platforms.

## 5. Risks

### 5.1 Id compatibility

- odb9 ids are `EOG09260001`-style; odb10 and odb12 ids are `<orthogroup>at<taxid>` (`100957at4751`,
  `364at4751`; L5, L7). The id sets are different orthogroup definitions (fungi_odb10 758 vs dikarya
  1,312 markers), so there is no 1:1 mapping (section 1.1 row on mapping); a genome annotated with odb9
  and one annotated with odb10 share zero `BUSCO:` ids.
- Consequences: `funannotate-compare` phylogeny and ortholog tables use shared ids across genomes
  (`compare.py:1140-1165`); the shipped outgroup proteomes are `EOG`-keyed; mixing namespaces gives an
  empty single-copy set without an error. Mitigation: Phase 1 item 4 (namespace tag), a check in
  `funannotate-compare` (out of scope here), and keeping odb9 as the default in 1.10.
- Existing annotated genomes keep their ids; re-annotation with odb10 changes the `BUSCO` column and
  the `note` qualifiers and therefore any downstream diff.

### 5.2 Offline use and image size

- The vendored path is fully offline after `setup`. Upstream BUSCO needs `--offline --download_path`
  with pre-placed `<download_path>/lineages/<name>` directories (U20); without `--offline` it contacts
  busco-data at every run and auto-updates datasets (U2, U20). The datasets are larger (fungi_odb10 126 MB packed, 278 MB
  unpacked; ascomycota_odb10 204 MB; odb12.2 fungi 949 MB; L3) than odb9 (dikarya 63/301 MB). "All
  BUSCO lineages" in the full image (Docker/README.md:67) would grow by several GB if odb10 and odb12
  were both shipped; ship odb9 `dikarya` + `fungi_odb10` only.
- `busco-data.ezlab.org/v5/data` redirects to S3 (L1): any downloader must follow redirects; CI has a
  900 s budget for the BUSCO download (image-prune-test.yml:109).

### 5.3 macOS arm64

- The vendored path depends on a GCC-built Augustus for `--proteinprofile`
  (install_augustus_gcc_macos.sh:7); that is already handled. Upstream BUSCO is "not supported" on
  macOS by its authors (U29); its extra dependencies (metaeuk, sepp, bbmap, prodigal) do have
  osx-arm64 bioconda builds, except that the exact `sepp==4.5.5` pin was not checked (U6), and nobody
  has run BUSCO 6 on arm64 for this review. Until the pixi `osx-arm64` target (pixi.toml:176-179)
  solves with `busco`, the external path must be optional and Linux-first.

### 5.4 User-facing CLI

- `--busco_db` keeps odb9 names (`dikarya`, default) and gains odb10 names; `setup -i busco -b` loses
  the fixed `choices` list (validation message instead). `database --show-buscos` output changes.
- `funannotate test -t busco` gene-count assertion (1,500-1,800) was calibrated with odb9; may need
  re-calibration for odb10 (Phase 1 item 5).
- Error semantics: predict exits on a non-zero BUSCO exit (predict.py:2136-2142) but `annotate`'s
  `runBUSCO` relies on `runSubprocess` (`sys.exit(1)`) and reads the table without checking for zero
  hits; `funannotate-compare` already works around the silent-zero case (`outgroups.py:120-135`).

### 5.5 Upstream churn

- v5 -> v6 changed the default gene finder and datasets (section 1.1); any external integration must
  be version-gated and layout-detected, and `runBUSCOv5` shows how fast a hard-coded layout rots
  (library.py:6921-6923).

## 6. Cheap hygiene fixes (safe now, independent of the decision)

1. Delete `funannotate/aux_scripts/funannotate-BUSCO2-py2.py` and the `sys.version_info` branches
   (predict.py:2055-2062, library.py:6939-6942). Python is `>=3.10` (setup.py:21); the py2 copy lacks
   `--augustus_timeout` and `--db_gencode`, so if it were ever selected predict's command would fail
   with "unrecognized arguments". Also update CHANGES.md:25 (mentions "the two BUSCO2 helper scripts").
2. Delete `library.runBUSCOv5` (:6871-6934, no callers, wrong `-m` in genome mode) and
   `setupDB.download_busco_v5` (:102-131, no callers, never downloads). Keep a note in CHANGES.md.
3. Delete `library.genemark2busco` and `library.evidence2busco` (:1442-1489, no callers).
4. `library.runBUSCO`: run the script with `sys.executable` (as predict.py:2107 does) instead of the
   shebang `#!/usr/bin/env python`; fail with a message when `full_table_busco.tsv` is missing or has
   zero `Complete`/`Duplicated` rows instead of a traceback / silent empty annotation.
5. Fix the `BUCSCO` typo in the training-source text (predict.py:991, 1036, 1061) — it is written into
   `parameters.json`/the run log.
6. `docs/predict.rst:56-61`: remove the "EVM on BUSCO2 models" steps; describe the current path
   (tblastn + Augustus profile, retrain, protein re-validation). `docs/tutorials.rst:200`:
   `--show_buscos` -> `--show-buscos`.
7. `downloads.json`: document (or normalize at setup time) the folder-name mismatches
   (`sordariomycetes` -> `sordariomyceta_odb9`, `saccharomycetes` -> `saccharomyceta_odb9`); `-u`
   checks depend on the folder rename at setupDB.py:528.
8. Record the lineage provenance (`dataset.cfg` name, creation_date, number_of_BUSCOs) in
   `training_decisions.tsv` (`busco_training` row reason) and in `annotate`'s log.
9. `tests/test_removed_commands.py:40-43`: patch `funannotate.library.log` with `create=True` (or call
   `lib.setupLogging` in `setUp`); the test fails when run alone or first.
10. `predict.py:2093-2099`: the message says tblastn runs single-threaded, but the gate
    (`tblastn_single_thread`) is False for all current BLAST+ (library.py:958-966), so the message never
    appears with BLAST+ >= 2.10; either drop the gate or make the log line say what actually happens.
11. `funannotate check`: add the vendored script's own requirements (`new_species.pl`,
    `optimize_augustus.pl`, `etraining`, `hmmsearch`, `tblastn`, `makeblastdb`) if not already listed,
    since the script exits with its own messages otherwise (script :422-436, :525-535).

## 7. Local experiments (all under `$TMPDIR/busco_probe`, nothing installed, no tracked file touched)

Tools: pixi `default` env (`augustus 3.5.0`, `hmmer 3.4`, `tblastn`, Python 3.11.14). The vendored
script writes `run_<name>` into the current directory; the two runs below were started from the repo
dir via `pixi run`, and the resulting `run_odb10self/` and `run_odb12self/` were deleted (git status
clean afterwards).

1. `funannotate-BUSCO2.py -i fungi_odb10/ancestral -m proteins -l fungi_odb10 -o odb10self -c 8 -f`
   -> rc 0, `C:94.3%[S:94.3%,D:0.0%],F:1.3%,M:4.4%,n:758`, 6.3 s (L6).
2. Same with `fungi_odb12` -> rc 0 from the wrapper but `ERROR Impossible to read the lengths in
   .../lengths_cutoff`, `BUSCO analysis failed` (L7). `fungi_odb12` has `hmms/` (1,122), `prfl/`
   (1,122), `scores_cutoff`, `ancestral`, `ancestral_variants`, `no_hits`, `refseq_db.faa.gz`,
   `links_to_ODB12.txt`, `dataset.cfg` (`OrthoDB_version=12.1`, `dataset_version=02`,
   `number_of_BUSCOs=1122`, `number_of_species=144`, `ncbi_taxid=4751`), no `lengths_cutoff`.
3. `augustus --species=anidulans --proteinprofile=<prfl> config/busco_test.fa` with an odb10, an odb12
   and the odb9 `EOG0926004Z` profile: all run; only the odb9 profile predicts a gene in the test
   sequence (it is the marker the test file was made for); odb10/odb12 profiles give "Block ... not
   significant enough" warnings like odb9 (L8).
4. `curl -IL` sizes and `file_versions.tsv` (L1-L3).
