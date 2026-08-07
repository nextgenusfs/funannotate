# AGENTS.md

Instructions for AI coding agents working in this repository.

## CLI help text and docs must track argparse

funannotate's top-level dispatcher (`funannotate/funannotate.py`) does **not** let each
subcommand's own argparse parser handle `-h`/`--help`. Instead, `funannotate.py` intercepts
`-h`/`--help` for every subcommand (see the dispatch logic that prints `info[cmdName]['help']`)
and prints a hand-written triple-quoted string instead (`trainHelp`, `predictHelp`,
`updateHelp`, `annotateHelp`, `cleanHelp`, `sortHelp`, `maskHelp`, `fixHelp`, `remoteHelp`,
`setupHelp`, `compareHelp`, `outgroupHelp`, `iprscanHelp`, `utilHelp`, and the `util`-scoped
`statsHelp`/`gff2tblHelp`/`gff2protHelp`/`gbk2partsHelp`/`contrastHelp`/`tbl2gbkHelp`/
`bam2gff3Help`/`stringtieHelp`/`quarryHelp`/`gffrenameHelp`/`predictTRNAHelp`/`prot2genomeHelp`).

These strings are **not** auto-generated, so they silently drift out of sync with the real
`parser.add_argument(...)` calls in each module whenever an argument is added, renamed,
removed, or its default/choices/required-ness changes. The docs under `docs/*.rst`
(`commands.rst`, `predict.rst`, `update.rst`, `annotate.rst`, `compare.rst`, `utilities.rst`)
also embed copies of these same help blocks and drift independently.

**Whenever you add, remove, rename, or change the default/choices/`required` flag of an
argument in any of these files**, you must also update the corresponding hand-written Help
string in `funannotate/funannotate.py`, and the matching example block(s) in `docs/*.rst`:

- `funannotate/train.py` ↔ `trainHelp`
- `funannotate/predict.py` ↔ `predictHelp` ↔ `docs/predict.rst`, `docs/commands.rst`
- `funannotate/update.py` ↔ `updateHelp` ↔ `docs/update.rst`, `docs/commands.rst`
- `funannotate/annotate.py` ↔ `annotateHelp` ↔ `docs/annotate.rst`, `docs/commands.rst`
- `funannotate/clean.py` ↔ `cleanHelp`
- `funannotate/sort.py` ↔ `sortHelp`
- `funannotate/mask.py` ↔ `maskHelp`
- `funannotate/fix.py` ↔ `fixHelp`
- `funannotate/remote.py` ↔ `remoteHelp`
- `funannotate/setupDB.py` ↔ `setupHelp`
- `funannotate/compare.py` ↔ `compareHelp` ↔ `docs/compare.rst`
- `funannotate/outgroups.py` ↔ `outgroupHelp`
- `funannotate/aux_scripts/iprscan-local.py` ↔ `iprscanHelp`
- `funannotate/utilities/*.py` ↔ the matching `*Help` variable ↔ `docs/utilities.rst`
- `funannotate/aux_scripts/funannotate-p2g.py` ↔ `prot2genomeHelp` ↔ `docs/utilities.rst`

### Review checklist when touching CLI arguments

1. Read the module's full argparse block (from `argparse.ArgumentParser(` to
   `args = parser.parse_args(`) and list every `add_argument` flag, its `default`, `choices`,
   and whether `required=True` is actually set.
2. Read the corresponding Help string in `funannotate/funannotate.py` and diff it against that
   list: flag missing from the help text, a flag in the help text that no longer exists, a
   flag shown under "Required" that isn't `required=True` (or vice versa), a stale default or
   choices list, or a flag that argparse now marks deprecated but the help text doesn't (and
   vice versa).
3. Fix the Help string in place, preserving the existing section layout (Required / grouped
   input sections / Optional / ENV Vars) and the `.format(package_name, __version__)`
   interpolation. Keep column alignment consistent with the surrounding lines. Watch for
   literal tab characters in TSV-format descriptions (e.g. `GeneID\tName\tProduct`) — inside a
   non-raw Python triple-quoted string these must be written as `\\t`, otherwise they print as
   real tabs and misalign the help output.
4. Grep `docs/*.rst` for the same subcommand's usage block (`Usage:       funannotate <cmd>
   <arguments>`) and apply the same fix there. These files are plain text quotes of the CLI
   help, not auto-generated, so they need the identical manual edit.
5. Verify by importing the module and printing the string, e.g.:
   ```
   python -c "import funannotate.funannotate as f; print(f.predictHelp)"
   ```
   Confirm no stray literal tabs, correct alignment, and that every real flag appears.
6. Prefer running `funannotate <cmd> --help` (or the equivalent for `util` subcommands) over
   guessing indentation by eye when in doubt.

### Keep examples in docs runnable

When editing `docs/*.rst`, any `funannotate <cmd> ...` example command must use flags that
actually exist in that command's current argparse block (not renamed/removed ones), and must
not claim a flag is required when it is optional (or the reverse). If unsure whether an example
is still valid, check it against the module's `add_argument` calls, not against older prose.
