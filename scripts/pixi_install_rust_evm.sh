#!/usr/bin/env bash
# Install the Rust EVidenceModeler from hyphaltip/EVidenceModeler_rust
# into the active pixi/conda environment.
#
# This script is idempotent: if the Rust EVM binaries are already present
# in $CONDA_PREFIX/bin, it exits immediately.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

EVM_REPO="https://github.com/hyphaltip/EVidenceModeler_rust.git"
EVM_BRANCH="main"
SRC_DIR="${CONDA_PREFIX}/opt/evm-rust-src"
BIN_DIR="${CONDA_PREFIX}/bin"

# Rust EVM produces these 6 binaries:
BINS=(
    EVidenceModeler
    evidence_modeler
    partition_evm_inputs
    recombine_evm_outputs
    convert_EVM_outputs_to_GFF3
    gff3_file_to_proteins
)

# Skip if already installed
if [ -x "${BIN_DIR}/evidence_modeler" ]; then
    exit 0
fi

echo "[pixi_install_rust_evm] Cloning ${EVM_REPO} (branch ${EVM_BRANCH})..."
rm -rf "${SRC_DIR}"
git clone --depth 1 --branch "${EVM_BRANCH}" "${EVM_REPO}" "${SRC_DIR}"

echo "[pixi_install_rust_evm] Building release binaries..."
cargo build --release --manifest-path "${SRC_DIR}/Cargo.toml"

echo "[pixi_install_rust_evm] Installing binaries into ${BIN_DIR}..."
for bin in "${BINS[@]}"; do
    cp "${SRC_DIR}/target/release/${bin}" "${BIN_DIR}/${bin}"
    chmod +x "${BIN_DIR}/${bin}"
done

echo "[pixi_install_rust_evm] Done. Rust EVM binaries installed."
