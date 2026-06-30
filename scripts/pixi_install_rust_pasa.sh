#!/usr/bin/env bash
# Build and install Rust-enabled PASA from GitHub
# Clones from https://github.com/hyphaltip/PASApipeline (explore_optimize_AI branch)
# and installs into .pixi/env/default/opt/pasa-rust-3.0
#
# This script is idempotent: if PASA is already built, it exits immediately.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

PASA_INSTALL_PREFIX="${CONDA_PREFIX}/opt/pasa-rust-3.0"
PASA_SRC="${PASA_INSTALL_PREFIX}/src"
PASA_BIN="${PASA_INSTALL_PREFIX}/bin"
PASA_REPO="https://github.com/hyphaltip/PASApipeline"
PASA_BRANCH="explore_optimize_AI"

# Check if PASA is already built
if [ -x "${PASA_BIN}/pasa_rust" ] && [ -x "${PASA_SRC}/Launch_PASA_pipeline.pl" ]; then
    exit 0
fi

# Clone PASA if not already present
if [ ! -d "${PASA_SRC}" ]; then
    echo "[pixi_install_rust_pasa] Cloning PASA from ${PASA_REPO} (${PASA_BRANCH} branch)..."
    mkdir -p "${PASA_INSTALL_PREFIX}"
    git clone --branch "${PASA_BRANCH}" "${PASA_REPO}" "${PASA_SRC}"
fi

echo "[pixi_install_rust_pasa] Building Rust-enabled PASA from ${PASA_SRC}..."

# Build Rust components
cd "${PASA_SRC}"
mkdir -p "${PASA_BIN}"
make rust

echo "[pixi_install_rust_pasa] Building C++ components..."

# Build remaining C++ components using make targets from top-level Makefile
# but skip Rust components which are already built
cd "${PASA_SRC}"

# Build C++ pasa assembler
cd pasa_cpp && make && cp pasa "${PASA_BIN}/"
cd ..

# Build slclust
cd pasa-plugins/slclust && make && cp src/slclust "${PASA_BIN}/"
cd ../../

# Build cdbfasta and cdbyank
cd pasa-plugins/cdbtools/cdbfasta && make && cp cdbfasta "${PASA_BIN}/" && cp cdbyank "${PASA_BIN}/"
cd ../../../

# Build seqclean utilities
cd pasa-plugins/seqclean/mdust && make && cp mdust "${PASA_BIN}/"
cd ../../../

cd pasa-plugins/seqclean/psx && make && cp psx "${PASA_BIN}/"
cd ../../../

cd pasa-plugins/seqclean/trimpoly && make && cp trimpoly "${PASA_BIN}/"
cd ../../../

# Copy seqclean utilities if they exist
if [ -f "${PASA_SRC}/pasa-plugins/seqclean/seqclean/seqclean" ]; then
    cp "${PASA_SRC}/pasa-plugins/seqclean/seqclean/seqclean" "${PASA_BIN}/"
fi
if [ -f "${PASA_SRC}/pasa-plugins/seqclean/seqclean/cln2qual" ]; then
    cp "${PASA_SRC}/pasa-plugins/seqclean/seqclean/cln2qual" "${PASA_BIN}/"
fi
if [ -f "${PASA_SRC}/pasa-plugins/seqclean/seqclean/bin/seqclean.psx" ]; then
    cp "${PASA_SRC}/pasa-plugins/seqclean/seqclean/bin/seqclean.psx" "${PASA_BIN}/"
fi

echo "[pixi_install_rust_pasa] PASA built successfully at ${PASA_SRC}"
