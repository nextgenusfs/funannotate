# Multi-stage build: Base environment with all dependencies for funannotate
FROM continuumio/miniconda3 AS build

# Update mamba and conda-pack
RUN conda update -n base -c defaults --yes conda && \
    conda install -c conda-forge -n base --yes mamba conda-pack

# Install funannotate core dependencies (bioconda packages)
# Note: Using newer versions where available, legacy conda packages for older tools
#
# c-compiler/cxx-compiler/make/cmake/autoconf/automake/libtool/zlib/bzip2/
# libffi/sqlite are the build toolchain required by the source-built Rust
# forks (trinityrnaseq/PASApipeline/EVidenceModeler @ rust_optimize) invoked
# below via scripts/pixi_install_*.sh -- mirrors pixi.toml so pixi and
# Docker builds can't drift apart again.
RUN mamba create -c conda-forge -c bioconda \
    -n funannotate --yes \
    "python>=3.6,<3.9" \
    "biopython<1.80" \
    xlrd==1.2.0 \
    "trinity>=2.8.5,<3" \
    "codingquarry==2.0" \
    "proteinortho>=6.0" \
    goatools \
    matplotlib-base \
    natsort \
    numpy \
    pigz \
    pandas \
    psutil \
    requests \
    "scikit-learn<1.0.0" \
    scipy \
    seaborn \
    "blast>=2.12" \
    tantan \
    bedtools \
    hmmer \
    exonerate \
    "diamond>=2.0.5" \
    tbl2asn \
    blat \
    "trnascan-se>=2.0" \
    ucsc-pslcdnafilter \
    trimmomatic \
    raxml \
    iqtree \
    trimal \
    "mafft>=7" \
    hisat2 \
    "kallisto==0.46.1" \
    minimap2 \
    stringtie \
    "salmon>=0.9" \
    "samtools>=1.9" \
    glimmerhmm \
    bamtools \
    perl perl-yaml perl-file-which perl-local-lib perl-dbd-mysql \
    perl-clone perl-hash-merge perl-soap-lite perl-json \
    perl-logger-simple perl-scalar-util-numeric perl-math-utils perl-mce \
    perl-text-soundex perl-parallel-forkmanager perl-db-file perl-perl4-corelibs \
    ete3 \
    distro \
    rustc cargo \
    c-compiler cxx-compiler make cmake autoconf automake libtool \
    zlib bzip2 libffi sqlite \
    && conda clean -a -y

# Install funannotate Python package
SHELL ["conda", "run", "-n", "funannotate", "/bin/bash", "-c"]
RUN python -m pip install git+https://github.com/nextgenusfs/funannotate.git@target_1.9/rust_EVM_trinity_PASA

# Package with conda-pack and unpack to /venv *before* installing PASA_rust/EVM_rust,
# since those need to install directly into /venv/opt/... -- creating /venv/opt/...
# before this step would make the `mkdir /venv` below fail with "File exists".
RUN conda-pack --ignore-missing-files -n funannotate -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar && \
    /venv/bin/conda-unpack

# Build trinity/PASA/EVM Rust components using the same install scripts pixi
# uses (scripts/pixi_install_{trinity,pasa,evm}_rust.sh), so Docker and pixi
# builds can't drift apart. Each script is COPYed immediately before the RUN
# that uses it, so editing one script doesn't invalidate the cache for the
# others' (expensive) cargo/make layers. Scripts read $CONDA_PREFIX to
# determine the install prefix; point it at /venv (the unpacked env used at
# container runtime). Commit SHAs are pinned inside each script. The scripts
# use `return` for early exits (they're written to be sourced, matching
# pixi's [activation] mechanism), so they must be `source`d here too --
# `bash script.sh` would let `return` fall through instead of exiting.
WORKDIR /tmp

COPY scripts/pixi_install_trinity_rust.sh /tmp/pixi_install_trinity_rust.sh
RUN CONDA_PREFIX=/venv bash -c "source /tmp/pixi_install_trinity_rust.sh"

COPY scripts/pixi_install_evm_rust.sh /tmp/pixi_install_evm_rust.sh
RUN CONDA_PREFIX=/venv bash -c "source /tmp/pixi_install_evm_rust.sh"

COPY scripts/pixi_install_pasa_rust.sh /tmp/pixi_install_pasa_rust.sh
RUN CONDA_PREFIX=/venv bash -c "source /tmp/pixi_install_pasa_rust.sh"

# ============================================================================
# Final runtime stage
FROM debian:bullseye AS runtime

LABEL maintainer="Jason Stajich <jason.stajich@ucr.edu>" \
      description="funannotate with Rust-optimized PASA and EVidenceModeler"

# Copy conda environment from build stage
COPY --from=build /venv /venv

# Install runtime dependencies
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    snap \
    augustus \
    augustus-data \
    locales \
    locales-all \
    libgl1 \
    libssl-dev \
    libsqlite3-0 \
    mysql-client \
    procps \
    perl \
    && rm -rf /var/lib/apt/lists/* && \
    ln -s /usr/bin/snap-hmm /usr/bin/snap && \
    rm "/venv/bin/fasta" && \
    ln -s "/venv/bin/fasta36" "/venv/bin/fasta"

# Setup environment variables
# EVM binaries land directly in /venv/bin (pixi_install_evm_rust.sh installs
# there, matching pixi's $CONDA_PREFIX/bin convention); PASA binaries land in
# opt/pasa-rust/bin (PASApipeline's own install.sh convention, sibling of
# opt/pasa-rust/src which is what PASAHOME points at). EVM_HOME points at the
# single evm-rust/src tree, which already contains EvmUtils/ from the same
# checkout used to build the Rust binaries -- no separate Perl EVM clone.
ENV PATH="/venv/bin:/venv/opt/pasa-rust/bin:$PATH" \
    AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config" \
    EVM_HOME="/venv/opt/evm-rust/src" \
    PASAHOME="/venv/opt/pasa-rust/src" \
    FUNANNOTATE_EVM_ENGINE="rust" \
    TRINITYHOME="/venv/bin" \
    TRINITY_HOME="/venv/bin" \
    QUARRY_PATH="/venv/opt/codingquarry-2.0/QuarryFiles" \
    ZOE="/usr/share/snap" \
    USER="funannotate" \
    FUNANNOTATE_DB="/opt/databases"

# Verify installations
RUN evidence_modeler --version && \
    Launch_PASA_pipeline.pl --help | head -5 && \
    test -f "$EVM_HOME/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl" && \
    funannotate --version || true

# Create non-root user
RUN useradd -m -u 1000 funannotate
USER funannotate
WORKDIR /work

SHELL ["/bin/bash", "-c"]
CMD ["funannotate", "--help"]
