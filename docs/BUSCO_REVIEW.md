# Review: BUSCO use in `funannotate predict` (internal BUSCO2 vs external BUSCO v5/v6)

Branch reviewed: `review/busco-predict-external` (read-only analysis; no code changed).

Legend: statements about the internal code are from reading the files, with file:line.
Statements about BUSCO v5/v6 are from memory and marked **(verify)** where I am not certain.
I did not run any benchmarks.

## 1. Summary

- `predict` uses BUSCO only as a **training-set generator** for Augustus, SNAP and GlimmerHMM
  when no transcript-derived models exist (or PASA fails its gate / yields too few models).
  The output that matters is `predict_misc/busco.final.gff3` (predict.py:2063), consumed as
  `FinalTrainingModels` (predict.py:2290, 2302, 2412; used at 2452, 2762+, 2824+).
- The "BUSCO" is a vendored, stripped-down BUSCO2/3-era script
  (`aux_scripts/funannotate-BUSCO2.py`, 3451 lines, plus a 3340-line py2 copy) running
  tblastn -> `augustus --proteinprofile` per marker, retrain, rerun on missing/fragmented.
  It is **not** BUSCO v5/v6 and shares no code or lineage data with them.
- External BUSCO v5 (with `--augustus`) and v6 behave differently. Several details that matter
  for us: v6 reportedly drops Augustus and uses odb12, so the "Augustus-trained-by-BUSCO"
  workflow disappears **(verify)**. Output layout and `full_table.tsv` columns differ from what
  `parseBUSCO2genome` expects, so external output cannot be dropped in as-is (see section 5).
- The internal BUSCO also drives `funannotate annotate` BUSCO annotation
  (annotate.py:1580 -> `lib.runBUSCO`, library.py:6942) and probably `compare`
  (library.py:10724 `busco_dictFlip`). Removing it is therefore bigger than `predict` alone.
- **Recommendation: Option C (either/or), phased.** Keep the internal path as default and
  frozen. Add a format-normalizer layer plus `--busco_results` (user-supplied external run)
  first; add optional "run external busco for me" later. Do not make external BUSCO a hard
  dependency (B) until a benchmark shows it is at least as good and v6 is stable in bioconda.

## 2. Comparison table

| Aspect | Internal (funannotate-BUSCO2) | External BUSCO v5 (`--augustus`) | External BUSCO v5.x/v6 (metaeuk/miniprot) |
|---|---|---|---|
| Gene finder for marker | tblastn region finding, then `augustus --proteinprofile` with per-marker `.prfl` (BUSCO2.py:1251, 1474) | Same idea (tblastn + augustus), uses odb10 `prfl` **(verify)** | metaeuk (v5 default for eukaryote genome) or miniprot (newer v5 / v6) **(verify versions)** |
| Seed/start parameters | `--busco_seed_species` (default `anidulans`, predict.py:223); silently falls back to `anidulans` (predict.py:2081-2084) | `--augustus_species`; BUSCO retrains species mid-run | None needed; no Augustus species |
| Retraining | etraining on complete single-copy, optional `optimize_augustus.pl` (BUSCO2.py:1378-1420); then rerun missing/fragmented | Same, `--long` enables optimization **(verify)** | No Augustus retraining output (Augustus not used) |
| Translation table | Augustus gets `--translation_table=N` via `-a` (predict.py:2129-2130); tblastn call has **no** `-db_gencode` (BUSCO2.py:1092) | Not exposed AFAIK **(verify)**; `--augustus_parameters` could pass it | Not exposed AFAIK; metaeuk/miniprot have genetic-code options but BUSCO may not forward them **(verify)** |
| Lineage data | `$FUNDB/<name>` BUSCO2-style dirs (`prfl/`, `scores_cutoff`, `lengths_cutoff`, `ancestral`), from a funannotate-hosted JSON of links (setupDB.py:523-548, 698-700); tree list in resources.py:197 | odb10 from busco-data.ezlab.org (`setupDB.download_busco_v5`, setupDB.py:102, appears unused for predict) | odb10 (metaeuk) or odb12 (v6) **(verify)** |
| Timeouts | per-Augustus-job kill (default 300 s, predict.py:239-243; BUSCO2.py:324-372, 2118) | none per job | none |
| Threads | tblastn forced to 1 thread when `version > "2.2.31"` as a string compare, i.e. only BLAST+ 2.2.32-2.9.x; 2.10+ is multi-threaded (BUSCO2.py:1087; predict.py:2097-2103); Augustus jobs threaded | BUSCO threads tblastn/augustus by job **(verify)** | metaeuk/miniprot multi-threaded, typically much faster **(verify)** |
| Classification | Complete / Duplicated / Fragmented / Missing using score + length sigma cutoffs (BUSCO2.py:1738-1766, 1886-1909) | Same concepts, odb10 cutoffs | Same labels, different evidence/scoring |
| `full_table` columns | `busco, status, contig, start, end, score, length` (no strand) (library.py:7059-7064) | `... start, end, strand, score, length, url, desc` **(verify)** | similar to v5 **(verify)** |
| Output layout | `run_<name>/full_table_<name>.tsv`, `augustus_output/{gffs,predicted_genes}`, `hmmer_output` | `<out>/run_<lineage>/full_table.tsv`, `busco_sequences/`, `augustus_output/` **(verify)** | `run_*/metaeuk_output` or `miniprot_output`, `busco_sequences/`, GFF in output **(verify)** |
| Python | py2 + py3 copies (diff ~283 lines) | py3 only | py3 only (v6 requires a recent Python **(verify)**) |
| Maintenance | vendored, frozen, ours | upstream churn | upstream churn, fastest |

## 3. What the internal BUSCO does (verified from code)

Invocation (predict.py:2062-2250):

1. Pick py3 or py2 script (2066-2073). Lineage dir is `$FUNDB/<busco_db>` (2074); default
   `dikarya` (predict.py:237). Check at predict.py:500-504 only tests the dir exists.
2. Seed species: if `--busco_seed_species` is not a valid Augustus species it silently uses
   `anidulans` (2081-2084). Separately, at ~predict.py:997 the error message contains a `{}`
   but no `.format(...)` (cosmetic bug).
3. Command: `-m genome --lineage <db> -c <cpus> --species <seed> -f --local_augustus <config>`;
   if `--table != 1` adds `-a --translation_table=N`; always adds `--augustus_timeout`
   (2110-2131). Output is cached by existence of `full_table_<aug_species>.tsv` (2087-2092).
   Failure of BUSCO is detected only by the absence of the table (2143-2147); exit code is ignored.
4. Inside BUSCO2.py (genome mode): makeblastdb + `tblastn -evalue 1e-3 -outfmt 7` with
   **num_threads=1** for BLAST > 2.2.31 (1087-1097; EVALUE_DEFAULT at 269). Candidate regions
   are chosen per marker (up to `--limit` regions, default capped 1-20, 3005+/3253), keeping
   regions >= 0.7 x best size (2538-2543), with flank = genome size/50 clamped to 5-20 kb (3268-3287).
5. Augustus pass 1 per (marker, region) with the seed species and the marker's
   `.prfl` (1251-1264). Output filtered by score/length cutoffs and HMMER (parse at 1701-1770).
6. Complete single-copy models are converted to GenBank and used to `new_species.pl` +
   `etraining` a `BUSCO_<abrev><rand>` species, optional `optimize_augustus.pl` (1378-1420).
   Missing and fragmented markers are re-predicted with that species (1455-1495).
7. Back in predict.py: bed of models (+/-100 bp) via `lib.parseBUSCO2genome` (library.py:7045),
   per-marker Augustus GFFs converted with `Converter2` (predict.py:2166-2177), renamed
   (`lib.fix_busco_naming`, library.py:1524), then **re-validated** by running internal
   BUSCO in `proteins` mode on the predicted proteins (predict.py:2194-2221), keeping only
   those still Complete (`getCompleteBuscos`, library.py:7016; `filterGFF3`, 2236).
8. Training: `busco.final.gff3` -> `GFF2GB ... 600` -> `trainAugustus` (predict.py:2452-2464),
   with `numTrainingSet` = 10 % (>1000), 20 % (<500), else 100 (2443ff). SNAP and GlimmerHMM
   train from the same file (predict.py:2762, 2824).

Ploidy > 1: `Duplicated` entries are also accepted and the best-scoring copy is picked, using
`dupBUSCO2gff`, which re-opens internal `hmmer_output/` and `augustus_output/predicted_genes/`
files (library.py:6977-7013, 7070-7085). This is tightly coupled to internal layout.

### Divergences from BUSCO v5/v6

- Gene finder and evidence: internal has a seed-species-dependent Augustus; v5/v6 default gene
  finders (metaeuk/miniprot) are protein-to-genome and seed-free **(verify)**.
- Marker data: BUSCO2-era datasets (our own tarballs; lengths/scores cutoffs, ancestral
  sequences) vs odb10/odb12. Counts in resources.py:197 are per old datasets (e.g. dikarya 1312,
  fungi 290); odb10 `fungi` has ~758 and `ascomycota`/`basidiomycota` ~1700/1764 markers
  **(verify numbers)**. More/better markers means more training candidates.
- Classification: re-scoring uses BUSCO2 sigma-length thresholds; v5 uses HMM score and
  length cutoffs per marker from odb10 with minor differences **(verify)**.
- Duplicates: internal selects best copy only for ploidy > 1; otherwise Duplicated are dropped.
  External the same for training purposes if we filter on `Complete`.
- No per-marker translation-table awareness in tblastn (no `-db_gencode`) even though Augustus
  gets the table; with codon reassignment, marker detection by tblastn may miss or truncate
  hits near reassigned codons. This is an unverified concern worth testing.

## 4. Evaluation against the decision criteria

### (1) Accuracy of the training set

- Internal strengths: tuned for fungi, two-pass retraining reduces seed bias; extra protein
  re-validation (predict.py:2194-2240) removes models that no longer score Complete; works with
  `--table` through Augustus.
- Internal weaknesses: seed-species bias remains for hits found in pass 1; markers from old
  datasets; silent seed fallback; non-Dikarya fungi (Mucoromycota, Chytrids, Microsporidia)
  rely on `fungi`/`microsporidia` datasets with few markers (resources.py:197: fungi 290,
  microsporidia 518), giving small training sets that can trip `--min_training_models`
  (predict.py:301, 2384-2418).
- External v5 `--augustus`: similar approach, richer odb10 data, likely more training genes.
  Retains Augustus-based gene structures, so format of GFF is Augustus-like.
- External metaeuk/miniprot GFFs: protein-spliced alignments, not ab-initio predictions.
  Exon/intron boundaries are protein-evidence based and UTRs/start-stop may be imprecise;
  these are usually good for CDS structure but incomplete ORFs at starts are more common
  **(verify empirically)**. Need CDS sanity filters (start Met, stop, in-frame, no internal stops)
  before `GFF2GB`. This is the same filtering `lib.gff2dict` + protein validation already does.
- Translation table: the current branch (`add_translationtable_support`) threads `--table`
  into BUSCO (2129), Augustus (2516), `gff2dict` (3051+). External BUSCO has no confirmed
  table support, so for `--table != 1` an external run's models could carry wrong stops;
  we would need to re-translate with the table and drop those with internal stops, or
  refuse external mode when `--table != 1`.
- Odd lineages: odb12 (v6) widens taxon coverage **(verify)**; internal can only use what we
  package.

### (2) Speed

- Internal: tblastn is pinned to 1 thread only on BLAST+ 2.2.32-2.9.x. The gate
  (BUSCO2.py:1087; predict.py:2097) compares version strings lexicographically, so 2.10+
  (all current releases) runs with all cpus, the opposite of the code comment's intent
  (see tests/test_tblastn_threads.py). On old BLAST+ this is the serial
  bottleneck on big genomes; Augustus per-marker jobs are parallel but numerous
  (markers x regions); two passes plus the second protein-mode BUSCO run (predict.py:2198)
  adds time. `--optimize_augustus` is separate from BUSCO and is not passed (the `long` mode
  exists in BUSCO2.py:1405 but predict does not pass `--long`).
- External metaeuk/miniprot: expected to be minutes per fungal genome and multithreaded
  **(verify)**; v5 `--augustus` has comparable runtime to ours.
- No measurements were made here; see benchmark plan.

### (3) Maintenance burden

- Internal: 6.8k lines vendored, py2/py3 duplicate (~283 diff lines), custom locking/timeouts,
  dataset tarballs hosted by us (setupDB.py:698-700), tblastn quirks. Support window
  `>=3.6,<3.12` (setup.py:21) is easy for internal code, but py2 copy is dead weight
  (predict.py:2066, library.py:6944 branch on `sys.version_info`, yet package requires >=3.6).
- External: requires BUSCO installed (conda-only convenience), version churn between v5.x
  and v6 (CLI flags, dataset names, output layout), output-format coupling
  (`parseBUSCO2genome` depends on 7-column table and internal directories).
  BUSCO v6 Python requirement may conflict with `<3.12` or with other bioconda pins
  **(verify)**. Plausible conda solver conflicts with augustus/hmmer/blast pins **(verify)**.
- Dead/buggy code: `lib.runBUSCOv5` (library.py:6876) is unused in the files I grepped
  (only `runBUSCO` is called, annotate.py:1580), its `genome` branch passes `-m proteins`
  (so does not run genome mode), and it hardcodes `run_<lineage>_odb10` path (6934-6936).
  This is a starting point but not a working implementation.

## 5. Options

### A. Keep internal
Pros: no new deps, translation-table path already wired, stable behavior.
Cons: stale datasets, py2 duplicate, tblastn thread gate that is a string compare, ever-growing maintenance of a fork.
Cheap wins even under A: delete the py2 script; pass `--long`/dataset choices explicitly;
fix error format at ~predict.py:997; log (not silently) the `anidulans` fallback; check the BUSCO
exit code; consider `-db_gencode` for tblastn with `--table`.

### B. Hard dependency on external BUSCO
Pros: removes vendored code; modern datasets. Cons: breaks `annotate` BUSCO annotation unless
also ported (annotate.py:1580), breaks offline lineage packaging (`funannotate setup -b`),
loses translation-table handling, couples to upstream churn. **Not recommended now.**

### C. Either/or, including user-supplied BUSCO output (recommended)

CLI sketch (predict):

- `--busco_mode {internal,external,auto}` (default `internal`).
- `--busco_results PATH` : a BUSCO run directory, `full_table.tsv`, a `.gff`/`.gff3`, or a `.zip`/`.tar.gz`
  of a run dir. Implies external training; skips running any BUSCO.
- `--busco_genome_ok` or reuse `-i`: genome used for the BUSCO run must be the same masked
  assembly (validate by contig names and lengths).
- `--busco_exe PATH` (phase 3): run external `busco` ourselves
  (`-m genome -l <lineage_odb> --augustus|--miniprot|--metaeuk`).
- Keep `--busco_db`, but allow an odb name (e.g. `fungi_odb10`) when external.

Normalizer layer (new module e.g. `funannotate/busco_import.py`, called from library.py):

1. `detect_busco_output(path)` -> unpack zip, find `run_*/full_table*.tsv`, GFF, `busco_sequences/`
   and tool provenance (augustus / metaeuk / miniprot) from dir names or `short_summary*.txt`.
2. `parse_full_table(path)` -> normalized records
   `{busco_id, status, contig, start, end, strand, score, length}`; handle both our 7-column
   and v5/v6 10-column (strand column) layouts. Current `parseBUSCO2genome`
   (library.py:7045-7064) and `getCompleteBuscos` (7016-7033, unpacks exactly 5 fields)
   would misparse the v5 table: `cols[5]` is strand in v5 not score.
3. `busco_models_to_gff3(records, run_dir, genome)` -> `busco_augustus.gff3`-equivalent,
   sourcing from `augustus_output/gffs/*.gff` (v5 augustus mode), `miniprot_output`/`metaeuk_output`
   GFF (others), or the top-level GFF if provided. Filter to `Complete` (and best copy of
   `Duplicated` if ploidy > 1), one model per BUSCO id.
4. Validate: contigs exist in genome; CDS divisible by 3 and translates under `--table` with
   no internal stops; start/stop sanity; minimum count; log counts per status. If
   `--table != 1` and tool is not Augustus-with-table, re-translate and filter, else refuse.
5. Write `predict_misc/busco.final.gff3` and rejoin at `FinalTrainingModels`.

Where predict.py / library.py would change:

- predict.py:220-243, 301-365: new flags; mutual exclusion checks near 471 (`conflict`)
  and 500-504 (lineage dir check must be skipped when external results supplied).
- predict.py:980-1010, 1044-1070, 1132-1170: RunModes stays `"busco"`; augustus
  functional check (1165) should not exit when external results supplied (no
  `--proteinprofile` needed).
- predict.py:2062-2250 `run_busco_training`: split into `run_internal_busco()`,
  `import_external_busco()`, shared tail (`fix_busco_naming`, validation, filter). Skip the
  protein re-validation (2194-2232) or replace with funannotate's own CDS validation when external.
- library.py: new parsers (above); retire/repair `runBUSCOv5` (6876); generalize
  `parseBUSCO2genome` (7045), `getCompleteBuscos` (7016), `dupBUSCO2gff` (6977),
  `fix_busco_naming` (1524, which assumes ID naming `Augustus%20prediction`).
- annotate.py:1580 and library.py:6942: leave internal for now; optionally later accept an
  external v5/v6 `full_table.tsv` for annotate.
- Docs/tests: fixtures of tiny v5 augustus, metaeuk, miniprot outputs.

Phased migration:

1. Phase 0 (cheap hygiene, no behavior change): the fixes listed under A; remove py2 copy if
   py2 is truly unsupported (requires python >=3.6).
2. Phase 1: normalizer + `--busco_results` (user-supplied only); unit tests with fixtures;
   no new dependency. Document that `--table != 1` requires re-translation filtering.
3. Phase 2: benchmark (below). If external-derived training >= internal on accuracy and
   faster, add `--busco_mode external` that shells out to `busco`, with version check and
   dataset resolution (odb10/odb12).
4. Phase 3: decide default. Only deprecate internal after annotate/compare no longer use it
   and v6 is stable in bioconda with compatible Python/Augustus pins.

## 6. Risks and open questions

- Does BUSCO v6 still ship/use Augustus and `--augustus_species`/`--long`? (I believe not.) **(verify)**
  If not, "external for training" means protein-evidence models (miniprot/metaeuk) only.
- Exact file names and columns for `full_table.tsv`, GFF, `busco_sequences/` in v5.8 and v6. **(verify)**
- Whether BUSCO exposes a genetic code for metaeuk/miniprot in genome mode. **(verify)**
- Whether protein-evidence GFFs train Augustus as well as Augustus's own models
  (start/stop completeness, UTR-less models, possible partial genes).
- Lineage mapping: `--busco_db dikarya` (BUSCO2 name) vs `fungi_odb10/odb12`; counts differ.
- Provenance/reproducibility: record BUSCO version, dataset and command in the training
  decision log (`lib.record_training_decision`, predict.py:2241 pattern).
- Hidden dependents: `annotate`, `compare` and `funannotate util`/`species` paths still use the internal tool.
- Python/bioconda constraints for BUSCO v6 vs `>=3.6,<3.12`. **(verify)**

## 7. Benchmark plan

Genomes (6-8): e.g. Aspergillus nidulans/fumigatus (Eurotiomycetes, close to the seed),
Neurospora crassa (Sordariomycetes), Saccharomyces/Candida (Saccharomycetes, CTG clade with
table 12 for Candida), Cryptococcus neoformans (Basidiomycota), Rhizopus or Batrachochytrium
(non-Dikarya), plus one with alternative table (Candida albicans table 12; a Mucor/mito-style
table only if available). Use each genome's reference annotation as held-out truth.

Arms: (A) internal; (C1) external v5 `--augustus`; (C2) external v5 metaeuk; (C3) external v6
miniprot/odb12; for each the same masked genome, same CPUs.

Metrics:
1. Number of Complete single-copy training models, number surviving CDS validation, fraction
   overlapping reference genes exactly (CDS structure).
2. Augustus accuracy on held-out genes (exclude training loci from the reference set; compute
   gene-/exon-/nucleotide-level sensitivity/specificity with `eval_accuracy` or gffcompare).
3. Final `predict` output quality: BUSCO completeness of final proteins, EVM consensus count,
   difference in gene counts vs reference.
4. Wall time and peak RAM per step (BUSCO stage vs total), at 8 and 32 CPUs.
5. For the table-12 genome: count of models with internal stops under `--table 12`, with and
   without re-translation filter.
6. Seed-bias test: run internal with `anidulans` vs a distant seed (e.g. `saccharomyces`) and
   compare training set overlap.

Success criterion: external arm within, say, 1-2 points of F1 on held-out genes and
no worse on table-12, and materially faster.

## 8. Recommendation

Adopt **Option C, phased**: freeze and tidy the internal path (default), add a validated
`--busco_results` import (works with run dirs, zips, `full_table.tsv`, GFF) through a
normalizer that hides v5/v6 layout differences, and defer making external BUSCO executable
or default until the benchmark and a v6 compatibility check. This keeps the translation-table
work intact, removes the biggest risk (silent format mismatch), and gives users with modern
datasets or non-standard lineages an escape hatch, without committing to upstream churn.
