# Multi-stage build: Base environment with all dependencies for funannotate
FROM continuumio/miniconda3 AS build

# Update conda-pack; install git for cloning external repos
RUN conda update -n base -c defaults --yes conda && \
    conda install -c conda-forge -n base --yes conda-pack git

# Install funannotate core dependencies (bioconda packages)
#
# This list mirrors pixi.toml's [dependencies] table exactly (channel-pinned
# entries use bioconda::pkg to match pixi's per-dependency `channel = ...`
# overrides) so pixi and Docker builds can't drift apart. A few extra
# packages beyond pixi.toml's list are pinned here deliberately: scipy,
# scikit-learn<1.0.0, pandas, psutil, requests, seaborn, and packaging come
# from funannotate's own setup.py install_requires (unpinned there), and are
# pre-installed here via conda so the `pip install` below is satisfied by
# these known-good pinned versions instead of resolving its own (potentially
# incompatible, e.g. scikit-learn>=1.0) versions from PyPI.
#
# Trinity is intentionally NOT installed as a conda package (nor is it in
# pixi.toml) -- it is built from source below via
# install_scripts/pixi_install_trinity.sh, same as pixi's install-trinity
# task. c-compiler/cxx-compiler/make/cmake/autoconf/automake/libtool/zlib/
# bzip2/libffi/sqlite are the build toolchain required by that and the other
# source-built Rust forks (PASApipeline/EVidenceModeler @ rust_optimize).
RUN conda create -c conda-forge -c bioconda \
    -n funannotate --yes \
    "python>=3.6,<3.9" \
    "biopython<1.80" \
    xlrd==1.2.0 \
    "goatools>=1.4.12,<2" \
    "matplotlib-base>=3.7.3,<4" \
    "natsort>=8.4.0,<9" \
    "numpy>=1.24.4,<2" \
    "pigz>=2.8,<3" \
    pandas \
    psutil \
    requests \
    "scikit-learn<1.0.0" \
    scipy \
    seaborn \
    packaging \
    "ete3>=3.1.3,<4" \
    "distro>=1.9.0,<2" \
    "rust>=1.80" \
    c-compiler cxx-compiler "make>=4" "cmake>=3.20" "autoconf>=2.69" "automake>=1.16" \
    "libtool>=2.4" "zlib>=1.2" "bzip2>=1.0" "libffi>=3.5.2,<4" sqlite \
    coreutils libgomp && \
    conda install -c conda-forge -c bioconda -n funannotate --yes \
    "codingquarry==2.0" \
    "augustus>=3.5.0" \
    "proteinortho>=6.3.6,<7" \
    "blast>=2.16.0,<3" \
    "tantan>=51,<52" \
    "bedtools>=2.31.1,<3" \
    "hmmer>=3.4,<4" \
    "exonerate>=2.4.0,<3" \
    "diamond>=2.0.5" \
    "tbl2asn>=25.7,<26" \
    "blat>=35,<36" \
    "trnascan-se>=2.0" \
    "ucsc-pslcdnafilter>=482,<483" \
    "glimmerhmm>=3.0.4,<4" \
    "bamtools>=2.5.3,<3" \
    "repeatmasker>=4.2.3,<5" \
    "repeatmodeler>=1.0.8,<2" \
    "repeatscout>=1.0.7,<2" \
    "recon>=1.8,<2" \
    "fasta3>=36.3.8,<37" \
    "bioconda::snap>=2017_3_1,<2018" \
    "bioconda::gmap>=2025.7.31,<2026" \
    "bioconda::eggnog-mapper>=2.1.13,<3" \
    "kmer-jellyfish>=2.3.0,<3" \
    ucsc-blat lighttpd "pblat>=2.5,<3" "transdecoder>=6.0.0,<7" && \
    conda install -c conda-forge -c bioconda -n funannotate --yes \
    "trimmomatic>=0.40,<0.41" \
    "fastp>=0.23,<1" \
    "raxml>=8.2.13,<9" \
    "iqtree>=3.1.1,<4" \
    "trimal>=1.5.1,<2" \
    "mafft>=7" \
    "hisat2>=2.2.2,<3" \
    "kallisto==0.46.1" \
    "minimap2>=2.30,<3" \
    "miniprot>=0.18,<0.19" \
    "stringtie>=3.0.3,<4" \
    "salmon>=1.0" \
    "samtools>=1.9" \
    "openjdk>=17" "bowtie2>=2.3.0,<3" && \
    conda install -c conda-forge -c bioconda -n funannotate --yes \
    perl "perl-yaml>=1.30,<2" "perl-file-which>=1.24,<2" "perl-local-lib>=2.29,<3" \
    "perl-dbd-mysql>=5.13,<6" "perl-clone>=0.46,<0.47" "perl-hash-merge>=0.302,<0.303" \
    "perl-soap-lite>=1.27,<2" "perl-json>=4.11,<5" "perl-logger-simple>=2.0,<3" \
    "perl-scalar-util-numeric>=0.40,<0.41" "perl-math-utils>=1.14,<2" "perl-mce>=1.902,<2" \
    "perl-text-soundex>=3.5,<4" "perl-parallel-forkmanager>=2.4,<3" "perl-db-file>=1.855,<2" \
    "perl-dbd-sqlite>=1.78,<2" perl-carp perl-uri perl-dbi \
    r-base r-cluster r-gplots r-fastcluster r-argparse r-ape r-phangorn r-tidyverse r-sm r-vioplot \
    bioconductor-qvalue bioconductor-ctc bioconductor-edger bioconductor-goseq \
    "bioconductor-go.db" bioconductor-dexseq && \
    conda clean -a -y

# Install funannotate Python package
SHELL ["conda", "run", "-n", "funannotate", "/bin/bash", "-c"]
RUN python -m pip install --no-cache-dir git+https://github.com/nextgenusfs/funannotate.git@target_1.9/rust_EVM_trinity_PASA

# Package with conda-pack and unpack to /venv *before* installing PASA_rust/EVM_rust,
# since those need to install directly into /venv/opt/... -- creating /venv/opt/...
# before this step would make the `mkdir /venv` below fail with "File exists".
RUN conda-pack --ignore-missing-files -n funannotate -o /tmp/env.tar && \
    mkdir /venv && cd /venv && tar xf /tmp/env.tar && \
    rm /tmp/env.tar && \
    /venv/bin/conda-unpack

# Build trinity/PASA/EVM Rust components using the same install scripts pixi
# uses (install_scripts/pixi_install_{trinity,evm,pasa}.sh -- mirrors
# pixi's `pixi run install-externals` task), so Docker and pixi builds can't
# drift apart. Each script is COPYed immediately before the RUN that uses it,
# so editing one script doesn't invalidate the cache for the others'
# (expensive) cargo/make layers. Scripts read $CONDA_PREFIX to determine the
# install prefix; point it at /venv (the unpacked env used at container
# runtime). Commit SHAs are pinned inside each script.
WORKDIR /tmp

COPY install_scripts/pixi_install_trinity.sh /tmp/pixi_install_trinity.sh
RUN CONDA_PREFIX=/venv bash /tmp/pixi_install_trinity.sh && \
    test -x /venv/bin/Trinity && \
    test -x /venv/bin/sam_to_read_coords

COPY install_scripts/pixi_install_evm.sh /tmp/pixi_install_evm.sh
RUN CONDA_PREFIX=/venv bash /tmp/pixi_install_evm.sh && \
    test -x /venv/bin/evidence_modeler && \
    test -f /venv/opt/evm/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl

COPY install_scripts/pixi_install_pasa.sh /tmp/pixi_install_pasa.sh
RUN CONDA_PREFIX=/venv bash /tmp/pixi_install_pasa.sh && \
    test -x /venv/opt/pasa/src/Launch_PASA_pipeline.pl && \
    test -d /venv/opt/pasa/src

# The conda-packaged bowtie2 2.5.5 falls back to slower baseline x86-64 code
# at runtime instead of using AVX2 (x86-64-v3) instructions; build from
# source with the local toolchain and shadow the conda binaries in-place
# (see install_scripts/pixi_install_bowtie2.sh for details).
COPY install_scripts/pixi_install_bowtie2.sh /tmp/pixi_install_bowtie2.sh
RUN CONDA_PREFIX=/venv bash /tmp/pixi_install_bowtie2.sh && \
    test -x /venv/bin/bowtie2 && \
    test -x /venv/bin/bowtie2-align-s-v256 && \
    test -x /venv/bin/bowtie2-align-l-v256

# ============================================================================
# Final runtime stage
# Debian bookworm (12) ships glibc 2.36, but the conda env built in the build
# stage contains packages that require __glibc >= 2.39 (grep conda-meta/*.json)
# and the source-built Trinity/PASA/EVM externals compiled via the conda-forge
# c/cxx toolchain require GLIBC_2.38. On a libc that old, every one of those
# compiled externals dies at exec with "GLIBC_2.38 not found" -- bamsifter,
# ParaFly, Inchworm and Chrysalis all fail, so genome-guided Trinity cannot run
# at all. Use Debian trixie (Debian 13), which ships glibc 2.41 and satisfies
# both requirements. Bump this together with a conda-env dependency bump if a
# future conda env raises its __glibc floor past 2.41.
FROM debian:trixie AS runtime

LABEL maintainer="Jason Stajich <jason.stajich@ucr.edu>" \
      description="Funannotate with Rust-optimized PASA and EVidenceModeler"

LABEL org.opencontainers.image.description="Eukaryotic Genome Annotation Pipeline"
# Copy conda environment from build stage
COPY --from=build /venv /venv

# Install runtime dependencies. Augustus and SNAP now come from the conda env
# (matching pixi.toml's `augustus`/`snap` dependencies) instead of apt --
# both ship their own binaries and data/config dirs directly in /venv, so the
# old apt augustus/augustus-data/snap packages and the snap-hmm->snap symlink
# are no longer needed.
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
    locales \
    locales-all \
    libgl1 \
    libssl-dev \
    libsqlite3-0 \
    default-mysql-client \
    procps \
    perl \
    && rm -rf /var/lib/apt/lists/* && \
    ln -sf /venv/bin/fasta36 /venv/bin/fasta

# Set locale for bioinformatics tools that require it
ENV LANG=C.UTF-8 \
    LC_ALL=C.UTF-8 \
    LANGUAGE=en_US:en

# Setup environment variables -- mirrors pixi.toml's [activation.env] plus
# the env vars bioconda's augustus/snap packages normally set via their own
# conda activate.d scripts (which conda-pack/manual ENV don't run), so
# Docker and pixi environments resolve externals the same way. EVM binaries
# land directly in /venv/bin (pixi_install_evm.sh installs there, matching
# pixi's $CONDA_PREFIX/bin convention); PASA binaries land in opt/pasa/bin
# (PASApipeline's own install.sh convention, sibling of opt/pasa/src which is
# what PASAHOME points at). EVM_HOME points at opt/evm, which already
# contains EvmUtils/ from the same checkout used to build the Rust binaries --
# no separate Perl EVM clone.
ENV PATH="/venv/bin:/venv/opt/pasa/bin:$PATH" \
    AUGUSTUS_CONFIG_PATH="/venv/config" \
    AUGUSTUS_SCRIPTS_PATH="/venv/bin" \
    AUGUSTUS_BIN_PATH="/venv/bin" \
    EVM_HOME="/venv/opt/evm" \
    PASAHOME="/venv/opt/pasa/src" \
    FUNANNOTATE_EVM_ENGINE="rust" \
    TRINITYHOME="/venv/opt/trinityrnaseq" \
    TRINITY_HOME="/venv/opt/trinityrnaseq" \
    QUARRY_PATH="/venv/opt/codingquarry-2.0/QuarryFiles" \
    ZOE="/venv/share/snap" \
    USER="funannotate" \
    FUNANNOTATE_DB="/opt/databases" \
    LD_LIBRARY_PATH="/venv/lib"

# Create non-root user
RUN useradd -m -u 1000 funannotate && \
    mkdir -p /work && chown funannotate:funannotate /work

# Verify installations as the runtime user. Use bash with pipefail so the
# Launch_PASA_pipeline.pl|head pipeline reports a real failure, and drop the
# old trailing `|| true` so a broken externals install actually fails the build.
USER funannotate
WORKDIR /work
SHELL ["/bin/bash", "-c"]
RUN set -euo pipefail && \
    test -x /venv/bin/evidence_modeler && \
    /venv/opt/pasa/src/Launch_PASA_pipeline.pl 2>&1 | head -5 && \
    test -f "$EVM_HOME/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl" && \
    funannotate --version

CMD ["funannotate", "--help"]
