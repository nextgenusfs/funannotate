# Spike: moving `funannotate compare` into its own package

Status 2026-10-23. Static analysis of the repo at `master`; nothing has been moved. The compare test (`funannotate test -t compare`)
is queued in `image-prune-test.yml` to prove the current behaviour in the image.

## What compare is made of

- `funannotate/compare.py`: 1,662 lines. Imports only `goatools`, `Bio`, `natsort`, `pandas`, `funannotate.library` and `funannotate.resources`.
- `funannotate/library.py` (12,688 lines): compare calls 52 names directly. Following calls inside `library.py` (an `ast` walk, so
  it can over-include a name that is also an attribute elsewhere) gives a closure of **79 names, 2,266 lines**:
  - **51 names, 1,306 lines, used by nothing else** (compare-only): the phylogeny and orthology code (`ortho2phylogeny`,
    `ReciprocalBlast`, `gb2proteinortho`, `drawPhyMLtree`, `align2Codon`, `rundNdSexhaustive`, `rundNdSestimate`, `parsedNdS`),
    plotting (`drawStackedBar`, `donutplot`, `drawHeatmap`, `drawbarplot`, `distance2mds`), annotation parsing
    (`getGBKannotation`, `genomeStats`, `parseGOterms`, `pfam2dict`, `busco_dictFlip` ...) and the whole HTML report
    (`HEADER`, `FOOTER`, `INDEX`, `SUMMARY`, `PHYLOGENY`, ... `CITATION`).
  - **28 names, 960 lines, shared with annotate/update/predict/fix/gff tools**: the real coupling. The big one is
    **`gb_feature_add2dict` (441 lines)**, the GenBank parser that `annotate`, `update`, `remote`, `gbk2parts` also use, with `getID` (97),
    `runSubprocess` (67), `process_handle` (34), `runMultiProgress` (32), `list_columns`, `CheckDependencies`, `RevComp`, `translate`,
    `setupLogging`, `copyDirectory` and about a dozen tiny ones (`SystemInfo`, `get_version`, `which_path`, `colr` ...).
- Other files that are compare-only: `funannotate/stackedBarGraph.py` (308), the `funannotate/html_template/` css and js,
  `docs/compare.rst`, the compare block in `funannotate/test.py`, and the `compare` entries in `funannotate.py`
  (help text, command table, `-t` choices).
- Related but a separate decision: `funannotate/outgroups.py` (121), `setupDB.py` `outgroupsDB` and `-i busco_outgroups`, and the
  `FUNANNOTATE_DB/outgroups` directory. They exist only for compare's ML tree outgroups. `funannotate util contrast` is unrelated.

Python modules imported only by the compare-only helpers: `matplotlib` (and `pylab`), `seaborn`, `sklearn.manifold`,
`sklearn.metrics.pairwise`, `Bio.Phylo.Consensus`. External programs compare checks for: `mafft`, `trimal`, `proteinortho`,
`ete3`, `phyml` (dN/dS), `find_enrichment.py` (goatools), and `iqtree` or `raxmlHPC-PTHREADS`; plus `diamond` for reciprocal hits.

## What it would take out of the main environment (from `pixi.lock`)

Packages reachable only from `ete3`, `iqtree`, `raxml`, `mafft`, `proteinortho`, `trimal`:

| Platform | packages | compressed |
|---|---|---|
| linux-64 | 43 | 209 MB |
| osx-arm64 | 37 | 153 MB |

The big ones: `qt-main` 52 MB, the LLVM libraries (`libllvm21` 44, `libllvm20` 44, `libclang-cpp20.1` 21, `libclang13` 12), `pyqt`.
That is the Qt and LLVM chain that `qdoc` and `libLTO` kept alive in the Linux image, and the chain behind the `ete3`,
`augustus`, `transdecoder` conflicts in the `linux-aarch64` test solve. (`goatools` stays: annotate and the gff tools use it.
`scikit-learn` and `seaborn` look compare-only in `library.py` but need a check across the other modules before they come out.)

## Extraction plan

1. **Isolate in place, no behaviour change.** `funannotate/compare/` subpackage; move the 51 compare-only names and `stackedBarGraph.py`
   and `html_template/` there; keep the 28 shared names in `library.py` and import them from there. Done when
   `funannotate test -t compare` gives identical output (compare the output folder, not just exit status). Needs a code review of
   the boundary (what moved, what was copied, import cycles, behaviour of `HEADER`/`FOOTER`/`CITATION` strings).
2. **Decide the shared-code contract.** Options: (a) `funannotate-compare` depends on `funannotate` and imports
   `funannotate.library` (simple, but pins the two together and keeps `funannotate` big); (b) copy the ~960 shared lines
   (mainly `gb_feature_add2dict`); (c) put a small GenBank-reading module in its own tiny package both use. Recommended: (a) for
   the first release, (c) later if it hurts.
3. **Extract** with `git filter-repo` to keep history; publish as `funannotate-compare` (PyPI, then Bioconda); own CI that runs the
   compare test data (`test-compare.tar.gz`, 1 MB in the `test-data-v1` release).
4. **Cut over in funannotate.** `funannotate compare` becomes a shim that runs the new tool if installed, otherwise prints how to
   install it; drop `-t compare` from `-t all`; remove `ete3`, `iqtree`, `raxml`, `mafft`, `proteinortho`, `trimal` from `pixi.toml`
   and re-lock; decide where `outgroups` and `busco_outgroups` live.
5. **Deprecate the shim** one release later.

## Risks

- `gb_feature_add2dict` is the dependency; if compare must work without `funannotate` installed, it has to be copied or shared.
- Breaking the `funannotate compare`, `outgroups`, `setup -i busco_outgroups` and `test -t compare` command lines; docs and citation text.
- Compare's output (HTML report, trees) must stay byte-for-byte comparable; the extraction must not change plotting dependencies' versions.
- The analysis is static: dynamic lookups (`getattr`, string-built names) would not show.

## Revision 2026-10-24: decisions so far and what they change

Decided (owner):
1. **Audience is people running 5 to 50 genomes** they annotated themselves. So compare cannot simply be deleted in favour of BFD:
   BFD needs Nextflow, an HPC tool stack and its own `samples.csv`/locus-tag conventions, and has no `funannotate annotate` ingest.
2. **Sections to keep: Pfam, GO (with enrichment), CAZy, MEROPS.** Phylogeny is optional and can be handed to PHYling instead of the
   built-in iqtree/raxml route. A workflow layer (Nextflow) is out of scope here, because it would repeat what BFD does.
3. **Input stays GenBank (`.gbk`) for the first version of the separate tool.** `annotations.txt` as an input is a later option.

Checked with an `ast` walk of `library.py` from every `lib.*` call in `compare.py` (32 roots, 56 functions reachable):
`gb_feature_add2dict` IS reachable, through `gb2proteinortho` (which builds the per-genome protein, transcript and GFF files for
ProteinOrtho from the `.gbk`), and `getID` through `getGBKannotation`. `annotationtable` is not reachable. So the shared-parser
coupling in this spike stands, and an agent's claim during the 2026-10-24 review that compare never uses it was wrong.

What follows from the decisions:
- Because the input is `.gbk`, the new tool needs its own copy of the GenBank feature parser (`gb_feature_add2dict`, 441 lines, plus
  `getID`). Contract option (c) of step 2 (copy now, share later) is the likely one; (a) would keep funannotate installed with it.
- If orthology is limited to Pfam/GO/CAZy/MEROPS plus an optional phylogeny, ProteinOrtho is only needed for the ortholog tables.
  Dropping the dN/dS chain (`ete3`, `phyml`, `codeml`), `outgroups` and the RAxML bootstrap path removes the Qt/LLVM chain from the
  funannotate environment (see the size table above).
- The phylogeny hand-off to PHYling needs an interface: either compare writes the per-genome proteomes (`.faa` from
  `gb2proteinortho`) in the layout PHYling reads, or the docs show how to run PHYling on them. Not designed yet.

Still open (owner):
- Are COG, InterPro, secretome/SignalP, transcription factors and secondary-metabolite counts kept or dropped? They were not named.
- Package name and repository; does `funannotate compare` stay as a shim for one release, or become an error that names the new tool?
- Where `funannotate outgroups` and `setup -i busco_outgroups` go if the built-in phylogeny is dropped.
- Whether BFD should later ingest `annotate_results` columns (GO, InterPro) so the Exophiala-style dashboard can read them.
- Owner note 2026-10-24: InterPro will be added to BFD in future, but BFD and compare differ for real (audience, scale, setup), so
  that does not by itself remove InterPro from compare; whether compare keeps its InterPro section is still to be decided.

### Owner answers, second round (2026-10-24)

- **Sections:** all of them stay (COG, InterPro, secretome/SignalP, TF and secondary-metabolite counts, as well as Pfam/GO/CAZy/MEROPS).
- **Repository:** `funannotate-compare`, ideally under `nextgenusfs`; that needs a conversation with Jon Palmer first. Start it as
  `hyphaltip/funannotate-compare` and transfer later (a transfer keeps history, issues and redirects).
- **Removal in 1.10:** remove `funannotate compare` outright; the command stays in the table only to print a warning and how to
  install and run the new package, then exit non-zero. No shim that runs the new tool.
- **Outgroups (checked in the code):** `funannotate outgroups`, `funannotate database --show-outgroups` and `funannotate setup -i
  busco_outgroups` serve only `funannotate compare --outgroup`, which names a species from `$FUNANNOTATE_DB/outgroups` (a
  `<species>_buscos.fa` file) to root the BUSCO-based ML tree (`compare.py` lines 161-180, 1086-1090, 1483-1500). Nothing in predict,
  train, update or annotate reads them. They belong to compare's tree building, so they leave with it or are dropped once the tree
  is done by PHYling.
- **Phylogeny:** yes, hand it to PHYling. Add an `nf_phyling` folder (a small Nextflow workflow for running PHYling on a set of
  proteomes) to the new project. Keep it small so it does not repeat BFD.
- **Outgroups move to `funannotate-compare` (owner, 2026-10-24).** `funannotate outgroups`, the outgroups listing and the
  `busco_outgroups` download become commands and data of the new package, together with `compare --outgroup`. In funannotate 1.10 the
  `outgroups` command, `database --show-outgroups` and `setup -i busco_outgroups` follow the same rule as `compare`: they print a
  warning with instructions and exit non-zero. The `$FUNANNOTATE_DB/outgroups` directory and its database entry are no longer
  created by funannotate; the new package must say where it keeps them (own data directory or an environment variable) and how an
  existing `$FUNANNOTATE_DB/outgroups` is reused. The `outgroups` code (`funannotate/outgroups.py`, 121 lines) runs BUSCO via
  library functions, so those need the same copy-or-import decision as the GenBank parser.

### PHYling hand-off (defined 2026-10-24)

Written by `funannotate-compare phyling` (repo `hyphaltip/funannotate-compare`, `docs/PHYLING.md`): `pep/<genome>.fa.gz`,
`cds/<genome>.fa` (optional) and `samples.tsv`, one file per genome with the compare genome name as the file name, built from the
`<name>.faa` and `<name>.transcripts.fa` files that `gb2proteinortho` already writes. `stajichlab/nf_phyling` already exists and reads
`.fa.gz` per taxon for protein and `.fa` per taxon for CDS, so no `nf_phyling` folder is created; the new repo only documents the call.
To check: transcripts in frame with proteins and any trailing stop codon, taxon naming from file names.

