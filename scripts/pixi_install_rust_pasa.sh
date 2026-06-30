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

# Build C++ components
cd "${PASA_SRC}/pasa_cpp"
make
cp pasa ../bin/.

cd "${PASA_SRC}/pasa-plugins/slclust"
make
cp src/slclust ../../bin/.

cd "${PASA_SRC}/pasa-plugins/cdbtools/cdbfasta"
make
cp cdbfasta ../../../bin/.
cp cdbyank ../../../bin/.

cd "${PASA_SRC}/pasa-plugins/seqclean/mdust"
make
cp mdust ../../../bin

cd "${PASA_SRC}/pasa-plugins/seqclean/psx"
make
cp psx ../../../bin

cd "${PASA_SRC}/pasa-plugins/seqclean/trimpoly"
make
cp trimpoly ../../../bin

# Copy seqclean utilities
cp "${PASA_SRC}/pasa-plugins/seqclean/seqclean/seqclean" "${PASA_BIN}/"
cp "${PASA_SRC}/pasa-plugins/seqclean/seqclean/cln2qual" "${PASA_BIN}/"
cp "${PASA_SRC}/pasa-plugins/seqclean/seqclean/bin/seqclean.psx" "${PASA_BIN}/"

echo "[pixi_install_rust_pasa] PASA built successfully at ${PASA_SRC}"
