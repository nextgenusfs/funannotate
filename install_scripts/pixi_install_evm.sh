#!/usr/bin/env bash
# Install Rust-enabled EVidenceModeler into the active pixi/conda environment.
#
# The Rust rewrite now lives in the main EVidenceModeler repo (the standalone
# EVidenceModeler_rust repo is retired) under the `rust_optimize` branch (also
# its default branch): https://github.com/hyphaltip/EVidenceModeler
#
# This script first checks for a local checkout (../EVidenceModeler), then
# falls back to cloning from GitHub if not available. The Rust workspace lives
# under evm/ (Cargo.toml at evm/Cargo.toml); there is no scripts/install.sh in
# this repo, so binaries are built with `cargo build --release` directly and
# copied into $CONDA_PREFIX/bin. The Perl EvmUtils/ tree ships in the same
# repo (needed unconditionally by funannotate/predict.py for the
# augustus/genemark -> EVM GFF3 converters, regardless of engine -- see
# scripts/benchmark_predict_evm.sh) and gets copied to opt/evm/src, which
# is what EVM_HOME should point at.
#
# This script is idempotent: if the Rust EVM binaries and EvmUtils are already
# present, it exits immediately.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

EVM_INSTALL_PREFIX="${CONDA_PREFIX}/opt/evm"
EVM_SRC="${EVM_INSTALL_PREFIX}"

EVM_BINARIES=(
    evidence_modeler
    EVidenceModeler
    convert_EVM_outputs_to_GFF3
    partition_evm_inputs
    recombine_evm_outputs
    gff3_file_to_proteins
    augustus_to_evm_gff3
    gff3_gene_prediction_file_validator
)

# Skip if already installed
if [ -x "${CONDA_PREFIX}/bin/evidence_modeler" ] && \
   [ -f "${EVM_SRC}/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl" ]; then
    echo "[pixi_install_evm_rust] Rust EVidenceModeler already installed"
    exit 0
fi

# Try to use local EVidenceModeler checkout first
EVM_LOCAL="../EVidenceModeler"
EVM_REPO="https://github.com/hyphaltip/EVidenceModeler"
# Tracks the rust_optimize branch by default. Override with EVM_RUST_COMMIT to
# pin a specific commit for reproducibility; see
# https://github.com/hyphaltip/EVidenceModeler/commits/rust_optimize
EVM_COMMIT="${EVM_RUST_COMMIT:-rust_optimize}"
EVM_CLONE_DIR=""
EVM_SOURCE_DIR=""

if [ -d "${EVM_LOCAL}/evm" ] && [ -f "${EVM_LOCAL}/evm/Cargo.toml" ]; then
    echo "[pixi_install_evm] Using local EVidenceModeler checkout from ${EVM_LOCAL}..."
    EVM_SOURCE_DIR="${EVM_LOCAL}"
else
    echo "[pixi_install_evm] Local EVidenceModeler not found or missing evm/, cloning from ${EVM_REPO} (commit ${EVM_COMMIT})..."
    EVM_CLONE_DIR=$(mktemp -d)
    trap "rm -rf '${EVM_CLONE_DIR}'" EXIT

    git init -q "${EVM_CLONE_DIR}"
    git -C "${EVM_CLONE_DIR}" fetch --depth 1 "${EVM_REPO}" "${EVM_COMMIT}"
    git -C "${EVM_CLONE_DIR}" checkout -q FETCH_HEAD
    EVM_SOURCE_DIR="${EVM_CLONE_DIR}"

    if [ ! -f "${EVM_SOURCE_DIR}/evm/Cargo.toml" ]; then
        echo "[pixi_install_evm] ERROR: evm/Cargo.toml not found in cloned ${EVM_REPO}@${EVM_COMMIT}" >&2
        exit 1
    fi
fi

echo "[pixi_install_evm] Building EVidenceModeler from ${EVM_SOURCE_DIR}..."

SAVED_CONDA_PREFIX="$CONDA_PREFIX"
SAVED_PATH="$PATH"

(
  cd "${EVM_SOURCE_DIR}/evm"
  cargo clean 2>/dev/null || true
  cargo build --release
)

export CONDA_PREFIX="$SAVED_CONDA_PREFIX"
export PATH="$SAVED_PATH"

# Install binaries to conda prefix (already on PATH, matches Trinity/PASA convention)
RELEASE_DIR="${EVM_SOURCE_DIR}/evm/target/release"
for binary in "${EVM_BINARIES[@]}"; do
    if [ -x "${RELEASE_DIR}/${binary}" ]; then
        cp "${RELEASE_DIR}/${binary}" "${CONDA_PREFIX}/bin/"
        echo "[pixi_install_evm] Installed ${binary}"
    else
        echo "[pixi_install_evm] WARNING: ${binary} not found after build"
    fi
done

# Install the Perl EvmUtils/PerlLib tree -- this is EVM_HOME's target
mkdir -p "${EVM_SRC}"
cp -r "${EVM_SOURCE_DIR}/EvmUtils" "${EVM_SRC}/"
if [ -d "${EVM_SOURCE_DIR}/PerlLib" ]; then
    cp -r "${EVM_SOURCE_DIR}/PerlLib" "${EVM_SRC}/"
fi
echo "[pixi_install_evm] Installed EvmUtils to ${EVM_SRC}"

# Clean up the temporary clone, if one was made
if [ -n "${EVM_CLONE_DIR}" ]; then
    rm -rf "${EVM_CLONE_DIR}"
    trap - EXIT
fi

echo "[pixi_install_evm_rust] Done. evidence_modeler is on PATH; EVM_HOME should be ${EVM_SRC}"
