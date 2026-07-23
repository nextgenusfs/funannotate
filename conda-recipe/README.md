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

This version includes Rust-optimized implementations of key components,
built from source by `build.sh` (not bioconda packages -- there is no
`pasa-rust`/`evidencemodeler-rust` package to depend on):
- **PASApipeline** (`hyphaltip/PASApipeline`, `rust_optimize` branch) - Fast transcript assembly and alignment
- **EVidenceModeler** (`hyphaltip/EVidenceModeler`, `rust_optimize` branch) - High-performance consensus modeling

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
- **Dependencies**: ~50 bioconda packages, plus PASApipeline/EVidenceModeler (rust_optimize) built from source during the build

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
- codingquarry ==2.0
- proteinortho >=6.0
- augustus, glimmerhmm (predictors; GeneMark/SNAP are license-gated, not
  packaged -- see the manual-dependency install guide)
- blast, diamond (protein searches)
- trinity, samtools, hisat2 (RNA-seq)
- hmmer, exonerate (sequence comparison)
- goatools (GO enrichment)
- And 30+ more specialized tools

**Built from source by `build.sh`** (not bioconda packages):
- PASApipeline (`hyphaltip/PASApipeline`, `rust_optimize` branch)
- EVidenceModeler (`hyphaltip/EVidenceModeler`, `rust_optimize` branch)

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
- `PASAHOME` — Points at `$CONDA_PREFIX/opt/pasa-rust/src` (built by `build.sh`)
- `EVM_HOME` — Points at `$CONDA_PREFIX/opt/evm-rust/src` (built by `build.sh`)

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
# - PASAHOME → $CONDA_PREFIX/opt/pasa-rust/src
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

### PASApipeline/EVidenceModeler build failures during `build.sh`

There is no separate `pasa-rust` package to install -- `build.sh` clones and
builds PASApipeline/EVidenceModeler (`rust_optimize` branch, pinned commits)
from source directly. If that step fails:

```bash
# Confirm network access to GitHub is available during the build (bioconda's
# sandboxed CI does NOT allow this -- this recipe is local-build-only)
git ls-remote https://github.com/hyphaltip/PASApipeline.git rust_optimize
git ls-remote https://github.com/hyphaltip/EVidenceModeler.git rust_optimize

# Confirm the compiler/build toolchain is present in the build env
# (c-compiler, cxx-compiler, make, cmake, rust/cargo, sqlite -- see
# requirements: build/host in meta.yaml)
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

## Related Repositories

`build.sh` builds these directly from source (pinned commits, see
`scripts/pixi_install_{pasa,evm}_rust.sh`) rather than depending on a
prebuilt package:
- [PASApipeline](https://github.com/hyphaltip/PASApipeline) (`rust_optimize` branch) — has its own standalone `conda-recipe/` upstream, but funannotate's recipe does not depend on it as a package
- [EVidenceModeler](https://github.com/hyphaltip/EVidenceModeler) (`rust_optimize` branch)

Note: the CI example above (`mamba build . -c conda-forge -c bioconda`) will
only work on a runner with outbound network access during `build.sh`, since
these are cloned from GitHub at build time -- this is not bioconda-CI
compatible as-is (see the note at the top of `meta.yaml`'s `requirements:`
section).

All three should be published together for best compatibility.

## References

- [Funannotate GitHub](https://github.com/nextgenusfs/funannotate)
- [Bioconda Project](https://bioconda.github.io/)
- [Conda-build Documentation](https://docs.conda.io/projects/conda-build/)
- [Conda Packaging Guide](https://docs.conda.io/projects/conda-build/en/latest/user_guide/)
- [Funannotate Publication](https://doi.org/10.12688/f1000research.8066.2)
