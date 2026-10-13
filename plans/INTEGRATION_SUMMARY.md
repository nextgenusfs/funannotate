# Rust EVM and PASA Integration Summary

## Overview
Successfully integrated Rust-enabled versions of EvidenceModeler and PASA as drop-in replacements for conda-installed versions in funannotate.

## Changes Made

### 1. Rust EvidenceModeler (EVM)
- **Removed**: conda package `evidencemodeler==1.1.1` from pixi.toml
- **Added**: Build script `scripts/pixi_install_rust_evm.sh` to compile and install Rust EVM on activation
- **Source**: https://github.com/hyphaltip/EVidenceModeler_rust
- **Environment**: `FUNANNOTATE_EVM_ENGINE=rust` (set in pixi.toml)
- **Binaries installed**: `evidence_modeler`, `partition_evm_inputs`, `recombine_evm_outputs`, `convert_EVM_outputs_to_GFF3`, `gff3_file_to_proteins`

### 2. Rust-Enabled PASA
- **Removed**: conda package `pasa==2.4.1` from pixi.toml
- **Added**: Build script `scripts/pixi_install_rust_pasa.sh` to compile and install PASA from local checkout
- **Source**: ~/projects/funannotate/PASApipeline (https://github.com/hyphaltip/PASApipeline)
- **Environment**: `PASAHOME=~/projects/funannotate/PASApipeline` (set in pixi.toml)
- **Rust components**: `pasa_rust`, `slclust_rust`, `cdbyank_rust`, `faidx_rust`
- **C++ components**: `pasa`, `slclust`, `cdbfasta`, `cdbyank`, `seqclean`, `mdust`, `psx`, `trimpoly`

## Configuration

### pixi.toml Changes
```toml
[activation]
scripts = ["scripts/pixi_setup_symlinks.sh", "scripts/pixi_install_salmon.sh", 
           "scripts/pixi_install_rust_evm.sh", "scripts/pixi_install_rust_pasa.sh"]

[activation.env]
TRINITY_HOME = "${CONDA_PREFIX}/bin"
TRINITYHOME = "${CONDA_PREFIX}/bin"
PASAHOME = "${HOME}/projects/funannotate/PASApipeline"
FUNANNOTATE_EVM_ENGINE = "rust"

[dependencies]
# Removed: evidencemodeler = "==1.1.1"
# Removed: pasa = "==2.4.1"
rust = ">=1.80"  # Required for building Rust components
```

## funannotate Integration

### predict.py
- Detects Rust EVM via `FUNANNOTATE_EVM_ENGINE` env var
- Validates `evidence_modeler` binary in PATH when engine is "rust"
- No `EVM_HOME` required for Rust engine

### train.py & update.py
- Accepts `--PASAHOME` argument to override env var
- Looks for `Launch_PASA_pipeline.pl` in `PASAHOME` root (v2.3.0+)
- Uses Rust PASA automatically when locally installed

### aux_scripts/funannotate-runEVM.py
- Auto-detects Rust engine via environment variable
- Maps Rust binaries to script interfaces
- Falls back to Perl EVM if Rust not available

## Build Process

### On pixi activation:
1. `pixi_install_rust_evm.sh`: Clones and builds Rust EVM via cargo → binaries in PATH
2. `pixi_install_rust_pasa.sh`: Builds from local checkout via make
   - Rust components: `cargo build --release` in pasa_rust/
   - C++ components: Makefile build system for pasa_cpp/, plugins, seqclean
   - All binaries installed to `${PASAHOME}/bin/`

## Testing

### Integration Tests (Passed ✓)
- PASA installation structure verified
- Rust binaries present and executable
- funannotate commands accept PASA arguments
- Launch_PASA_pipeline.pl found at expected location

### How to Test Locally
```bash
# Set environment
export PASAHOME=~/projects/funannotate/PASApipeline

# funannotate train
funannotate train -i genome.fasta -o output/ \
  --left R1.fastq.gz --right R2.fastq.gz \
  -s "Species name" --PASAHOME "$PASAHOME"

# funannotate update  
funannotate update -i output/ \
  --left R1.fastq.gz --right R2.fastq.gz \
  --PASAHOME "$PASAHOME"
```

## Benefits

- **Performance**: Rust implementations are optimized for speed
- **Compatibility**: Drop-in replacements maintain same interfaces
- **Local Control**: Building from source allows for customizations
- **Modern Tools**: Leverages Rust ecosystem for assembler and clustering

## Dependencies

- `rust >=1.80` (for building Rust components)
- Standard build tools: `gcc`, `g++`, `make` (for C++ components)
- All other funannotate dependencies unchanged

## Commits
- `23d52d3` Integrate Rust-enabled PASA from local checkout
- `587f5f1` Fix PASA C++ build process - add proper linking step
- Previous: `92ff9fe` pixi: install Rust EVM from hyphaltip/EVidenceModeler_rust
- Previous: `87cd2a4` fix Rust EVM CLI arg names and remove redundant -o flag
