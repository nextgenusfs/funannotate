# Funannotate - Conda Recipe

This directory contains the conda recipe for building and distributing funannotate with Rust-optimized PASA and EVidenceModeler components.

## Package Overview

**Name**: funannotate  
**Version**: 1.9.0  
**Build**: rust_1.9.0_0  
**Architecture**: Linux 64-bit only  

Funannotate is an eukaryotic genome annotation pipeline that combines:
- Gene prediction methods (Augustus, GeneMark, SNAP, GlimmerHMM, CodingQuarry)
- Transcript evidence (Trinity RNA-seq, PASA alignment assembly)
- Protein evidence (BLAST, DIAMOND, HMMER searches)
- Consensus modeling (EVidenceModeler)
- Functional annotation (Pfam, InterPro, dbCAN, MEROPS, EggNOG)

This version includes Rust-optimized implementations of key components:
- **pasa-rust** (>=3.0.0) - Fast transcript assembly and alignment
- **evidencemodeler-rust** (>=3.0.0) - High-performance consensus modeling

## Building the Package

### Prerequisites

- conda-build (or mamba build)
- Git
- Internet access (to clone dependencies)

### Local Build

```bash
# Build locally from this directory
cd funannotate-live/conda-recipe
conda build . -c conda-forge -c bioconda

# Faster with mamba
mamba build . -c conda-forge -c bioconda

# Build with specific output directory
conda build . -c conda-forge -c bioconda --croot /path/to/build
```

### Build Time and Size

- **Build time**: 30-60 minutes (most time spent installing dependencies)
- **Package size**: ~500MB compressed, ~2-3GB installed
- **Dependencies**: ~40 bioconda packages + pasa-rust + evidencemodeler-rust

### Testing the Build

```bash
# Install from local build
conda install --use-local funannotate

# Verify installation
funannotate --version
funannotate check --show-versions
```

## Publishing to Bioconda

To publish this package to [bioconda](https://bioconda.github.io/):

### Prerequisites

1. Fork [bioconda/bioconda-recipes](https://github.com/bioconda/bioconda-recipes)
2. Set up development environment:
   ```bash
   git clone https://github.com/YOUR_USERNAME/bioconda-recipes.git
   cd bioconda-recipes
   ```

### Steps

1. Create recipe directory:
   ```bash
   mkdir recipes/funannotate
   cp funannotate-live/conda-recipe/* recipes/funannotate/
   ```

2. Verify recipe:
   ```bash
   conda build recipes/funannotate -c conda-forge -c bioconda
   ```

3. Commit and push:
   ```bash
   git add recipes/funannotate
   git commit -m "Add funannotate 1.9.0-rust recipe"
   git push origin funannotate-1.9-rust
   ```

4. Create pull request on GitHub
   - Bioconda maintainers will review
   - CI will test on multiple platforms
   - Upon approval, package is built and uploaded to conda-forge

## Dependencies

### Build Dependencies
- Python >=3.6,<3.9
- setuptools

### Runtime Dependencies

**Core**:
- Python >=3.6,<3.9
- BioPython <1.80

**Bioconda Packages**:
- trinity >=2.8.5
- **pasa-rust >=3.0.0** ← Rust-optimized
- **evidencemodeler-rust >=3.0.0** ← Rust-optimized
- codingquarry ==2.0
- proteinortho >=6.0
- augustus, GeneMark, snap, glimmerhmm (predictors)
- blast, diamond (protein searches)
- trinity, samtools, hisat2 (RNA-seq)
- hmmer, exonerate (sequence comparison)
- goatools (GO enrichment)
- And 30+ more specialized tools

**System-level** (installed separately):
- perl (with modules: YAML, DBI, DBD-SQLite, DB_File)
- sqlite3 / mysql-client (for database)

## Installation via Conda

Once published to bioconda:

```bash
# Install from bioconda channel
conda install -c bioconda funannotate

# Or with specific version
conda install -c bioconda funannotate==1.9.0

# With Mamba (faster)
mamba install -c bioconda funannotate
```

## Environment Variables

When installed via conda, funannotate automatically sets:

- `FUNANNOTATE_EVM_ENGINE=rust` — Use Rust implementations
- `PASAHOME` — Set by pasa-rust activation scripts
- `EVM_HOME` — Set by evidencemodeler-rust if installed separately

Additional variables you may need to set:

```bash
# Database location (download or mount)
export FUNANNOTATE_DB=/path/to/databases

# Optional: Trinity configuration
export TRINITY_HOME=$(conda info --base)/bin

# Optional: Augustus configuration
export AUGUSTUS_CONFIG_PATH=$(conda info --base)/config/augustus
```

## Quick Start

After installation:

```bash
# Check dependencies
funannotate check --show-versions

# Download databases (first time only)
funannotate setup -i all --wget -d /path/to/databases

# Annotate a genome
funannotate predict \
  -i genome.fasta \
  -o output_directory \
  -s "Species name" \
  --rna_bam alignments.bam

# Annotate with functional information
funannotate annotate \
  -i output_directory \
  -s "Species name" \
  -d /path/to/databases
```

## Conda Activation/Deactivation

When you activate an environment with funannotate:

```bash
conda activate myenv

# Automatically sets:
# - PASAHOME → $CONDA_PREFIX/opt/pasa-rust-3.0/src
# - FUNANNOTATE_EVM_ENGINE=rust
# - Adds Rust tool binaries to PATH
```

When you deactivate:
```bash
conda deactivate

# PASAHOME and other env vars are unset
```

## Troubleshooting

### Build Fails with Dependency Issues

```bash
# Update conda/mamba and channels
conda update -n base conda
conda config --add channels conda-forge
conda config --add channels bioconda

# Try rebuild
mamba build . -c conda-forge -c bioconda --force-rebuild
```

### "pasa-rust not found" Error

```bash
# Ensure bioconda channel is configured
conda config --show channels

# If pasa-rust is not available, build it first
cd ../../../PASA_rust/conda-recipe
mamba build . -c conda-forge
# Then rebuild funannotate
```

### Installation Complains About Python Version

```bash
# Funannotate requires Python <3.9
conda create -n funannotate python=3.8 -c bioconda funannotate
```

### "ModuleNotFoundError: No module named 'funannotate'"

```bash
# Verify installation
python -c "import funannotate; print(funannotate.__file__)"

# If not found, check environment
conda list | grep funannotate
which funannotate
```

## Recipe Structure

```
conda-recipe/
├── meta.yaml           # Package metadata, dependencies, versions
├── build.sh            # Build/install script
└── README.md           # This file
```

### meta.yaml Sections

- **package**: Name and version
- **source**: Where to get source (git repo + branch)
- **build**: Build number, entry points
- **requirements**: Build, host, and runtime dependencies
- **test**: Import tests and command-line verification
- **about**: Homepage, license, summary, documentation
- **extra**: Maintainers, citations

### build.sh

- Runs during conda build process
- Installs Python package using setuptools
- Verifies installation with import test

## Development and Iteration

### Local Development Build

While developing funannotate:

```bash
# Build from local modified source
cd funannotate-live
mamba build conda-recipe -c conda-forge -c bioconda --force-rebuild

# Install locally
conda install --use-local funannotate

# Test changes
funannotate predict -i test_genome.fasta -o test_output/
```

### Updating Version

Edit `meta.yaml`:
```yaml
package:
  name: funannotate
  version: 1.9.1  # ← Update this
```

Edit `build:` section:
```yaml
build:
  number: 1      # ← Increment if same version, rebuild to fix
  string: rust_1.9.1_1
```

### Testing Before Publishing

```bash
# Full test suite
mamba build . -c conda-forge -c bioconda --test

# Test in fresh environment
conda create -n test-funannotate --use-local funannotate
conda activate test-funannotate
funannotate check
```

## CI/CD Integration

For automated conda builds (e.g., GitHub Actions):

```yaml
# .github/workflows/conda-build.yml
name: Conda Build
on: [push, pull_request]
jobs:
  build:
    runs-on: ubuntu-latest
    steps:
      - uses: actions/checkout@v2
      - uses: conda-incubator/setup-miniconda@v2
      - run: |
          mamba install -c conda-forge conda-build
          cd funannotate-live/conda-recipe
          mamba build . -c conda-forge -c bioconda
```

## Related Recipes

This recipe depends on:
- **pasa-rust** — Separate conda recipe in PASA_rust/conda-recipe/
- **evidencemodeler-rust** — Separate conda recipe in EVidenceModeler_rust/conda-recipe/

All three should be published together for best compatibility.

## References

- [Funannotate GitHub](https://github.com/nextgenusfs/funannotate)
- [Bioconda Project](https://bioconda.github.io/)
- [Conda-build Documentation](https://docs.conda.io/projects/conda-build/)
- [Conda Packaging Guide](https://docs.conda.io/projects/conda-build/en/latest/user_guide/)
- [Funannotate Publication](https://doi.org/10.12688/f1000research.8066.2)
