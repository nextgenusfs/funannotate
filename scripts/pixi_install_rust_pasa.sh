#!/usr/bin/env bash
# Build and install Rust-enabled PASA from local checkout
# at ~/projects/funannotate/PASApipeline
#
# This script is idempotent: if PASA is already built (bin/pasa_rust exists),
# it exits immediately.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

PASA_SRC="${HOME}/projects/funannotate/PASApipeline"
PASA_BIN="${PASA_SRC}/bin"

# Check if PASA is already built
if [ -x "${PASA_BIN}/pasa_rust" ] && [ -x "${PASA_SRC}/Launch_PASA_pipeline.pl" ]; then
    exit 0
fi

# Verify source directory exists
if [ ! -d "${PASA_SRC}" ]; then
    echo "[pixi_install_rust_pasa] ERROR: PASA source not found at ${PASA_SRC}"
    echo "[pixi_install_rust_pasa] Please clone https://github.com/hyphaltip/PASApipeline to ~/projects/funannotate/"
    exit 1
fi

echo "[pixi_install_rust_pasa] Building Rust-enabled PASA from ${PASA_SRC}..."

# Build Rust components
cd "${PASA_SRC}"
mkdir -p bin
make rust

echo "[pixi_install_rust_pasa] Building C++ components..."

# Build remaining C++ components using make targets from top-level Makefile
# but skip Rust components which are already built
cd "${PASA_SRC}"

# Build C++ pasa assembler
cd pasa_cpp && make && cp pasa ../bin/.
cd ..

# Build slclust
cd pasa-plugins/slclust && make && cp src/slclust ../../bin/.
cd ../../

# Build cdbfasta and cdbyank
cd pasa-plugins/cdbtools/cdbfasta && make && cp cdbfasta ../../../bin/. && cp cdbyank ../../../bin/.
cd ../../../

# Build seqclean utilities
cd pasa-plugins/seqclean/mdust && make && cp mdust ../../../bin
cd ../../../

cd pasa-plugins/seqclean/psx && make && cp psx ../../../bin
cd ../../../

cd pasa-plugins/seqclean/trimpoly && make && cp trimpoly ../../../bin
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
