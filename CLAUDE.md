# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

funannotate is a eukaryotic genome annotation pipeline (built primarily for fungi). It is a Python CLI tool installed as a package, dispatching subcommands to individual module scripts.

## Installation and Development Setup

```bash
# Install in editable/development mode
pip install -e .

# Required environment variable for all database-dependent commands
export FUNANNOTATE_DB=/path/to/funannotate/db
```

Python version requirement: `>=3.6.0, <3.12`

## Common Commands

```bash
# Run the CLI
funannotate <command> <arguments>

# Typical annotation workflow order:
funannotate clean -i genome.fasta -o genome_clean.fasta
funannotate sort -i genome_clean.fasta -o genome_sorted.fasta
funannotate mask -i genome_sorted.fasta -o genome_masked.fasta
funannotate train -i genome_masked.fasta -o output/ --left R1.fastq.gz --right R2.fastq.gz -s "Species name"
funannotate predict -i genome_masked.fasta -o output/ -s "Species name"
funannotate annotate -i output/ -s "Species name"

# Check dependencies
funannotate check --show-versions

# Database management
funannotate database
funannotate setup -d /path/to/db
```

## Linting

The project uses [trunk](https://trunk.io) for linting. Configured linters include `black`, `isort`, `ruff`, `shellcheck`, `shfmt`, and `markdownlint`.

```bash
trunk check         # lint staged/changed files
trunk fmt           # auto-format
```

## Architecture

### Entry Point and Dispatch

`funannotate/funannotate.py` is the single entry point (registered as the `funannotate` console script). It maintains an `info` dict mapping subcommand names to module paths and dispatches by importing the module and calling its `main(args)` function. For commands needing multiprocessing (e.g., `prot2genome`), it runs them as subprocesses instead.

### Module Layout

- `funannotate/library.py` — shared utility library used by nearly every module: logging setup, subprocess runners, FASTA/GFF3 parsers, file format converters
- `funannotate/resources.py` — static data/constants (BUSCO lineage tree, etc.)
- `funannotate/interlap.py` — interval overlap data structure used for genomic coordinate queries
- `funannotate/predict.py` — gene prediction (Augustus, GeneMark, SNAP, GlimmerHMM, CodingQuarry, EVM consensus)
- `funannotate/train.py` — RNA-seq training pipeline (Trinity + PASA)
- `funannotate/update.py` — PASA-based gene model refinement using RNA-seq
- `funannotate/annotate.py` — functional annotation (BUSCO, Pfam, dbCAN, MEROPS, EggNog, InterProScan, SignalP, Phobius)
- `funannotate/compare.py` — comparative genomics across multiple annotated genomes
- `funannotate/setupDB.py` / `funannotate/database.py` — database setup and management
- `funannotate/aux_scripts/` — helper scripts called as subprocesses (e.g., `funannotate-p2g.py` for protein-to-genome alignment, `funannotate-runEVM.py`, `hmmer_parallel.py`)
- `funannotate/utilities/` — format conversion utilities exposed under `funannotate util` (GFF3, TBL, GBK, BAM conversions)

### Database Dependency

All annotation-related commands require `$FUNANNOTATE_DB` pointing to a directory set up via `funannotate setup`. The database contains BUSCO models, Pfam HMMs, MEROPS, dbCAN, UniProt, and other reference data.

### Key External Tools

The pipeline wraps many external bioinformatics tools: Augustus, GeneMark, SNAP, GlimmerHMM, Trinity, PASA, HISAT2, minimap2, samtools, BLAST/diamond, HMMER, InterProScan, RepeatMasker, tbl2asn, SignalP, Phobius, EvidenceModeler (EVM).

### Current Branch

`add_translationtable_support` — adding translation table support to gene prediction/annotation.
