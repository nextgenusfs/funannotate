# Multi-stage build: Base environment with all dependencies for funannotate
FROM continuumio/miniconda3 AS build

# Update mamba and conda-pack
RUN conda update -n base -c defaults --yes conda && \
    conda install -c conda-forge -n base --yes mamba conda-pack

# Install funannotate core dependencies (bioconda packages)
# Note: Using newer versions where available, legacy conda packages for older tools
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

# Build EVidenceModeler_rust from source, installing straight into /venv/opt/evm-rust
WORKDIR /tmp
RUN git clone https://github.com/hyphaltip/EVidenceModeler_rust.git && \
    cd EVidenceModeler_rust && \
    (CONDA_PREFIX=/venv bash scripts/install.sh --install-prefix /venv/opt/evm-rust || \
     (cargo build --release && \
      mkdir -p /venv/opt/evm-rust/bin && \
      cp target/release/evidence_modeler target/release/EVidenceModeler \
         target/release/partition_evm_inputs target/release/recombine_evm_outputs \
         target/release/convert_EVM_outputs_to_GFF3 target/release/gff3_file_to_proteins \
         /venv/opt/evm-rust/bin/))

# EVidenceModeler_rust never ported EvmUtils/misc/*.pl (augustus/genemark -> EVM GFF3
# converters) -- funannotate/predict.py calls those Perl scripts unconditionally,
# regardless of engine, so a real Perl EVM checkout is required even in rust mode.
RUN git clone --depth 1 https://github.com/hyphaltip/EVidenceModeler.git /venv/opt/evm-perl && \
    test -f /venv/opt/evm-perl/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl

# Build PASA_rust from source, installing straight into /venv/opt/pasa-rust-3.0
RUN git clone https://github.com/hyphaltip/PASA_rust.git && \
    cd PASA_rust && \
    (CONDA_PREFIX=/venv bash scripts/install.sh --install-prefix /venv/opt/pasa-rust-3.0 || \
     (make rust && \
      mkdir -p /venv/opt/pasa-rust-3.0/bin /venv/opt/pasa-rust-3.0/src && \
      cp pasa_rust/target/release/* /venv/opt/pasa-rust-3.0/bin/ && \
      cp -r PerlLib PyLib pasa_conf schema scripts Launch_PASA_pipeline.pl /venv/opt/pasa-rust-3.0/src/))

# funannotate expects PASA binaries (seqclean, etc.) under $PASAHOME/bin, but the
# PASA_rust installer puts them in pasa-rust-3.0/bin, a sibling of pasa-rust-3.0/src
# (which is what PASAHOME points at) -- symlink so both layouts resolve.
RUN test -e /venv/opt/pasa-rust-3.0/src/bin || ln -s ../bin /venv/opt/pasa-rust-3.0/src/bin

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
ENV PATH="/venv/bin:/venv/opt/evm-rust/bin:/venv/opt/pasa-rust-3.0/bin:$PATH" \
    AUGUSTUS_CONFIG_PATH="/usr/share/augustus/config" \
    EVM_HOME="/venv/opt/evm-perl" \
    PASAHOME="/venv/opt/pasa-rust-3.0/src" \
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
