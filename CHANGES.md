# Changes

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
