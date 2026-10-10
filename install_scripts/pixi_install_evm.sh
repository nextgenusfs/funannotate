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

# Print which commit a local checkout is at, and warn when it is not the ref this script
# would otherwise clone. A stale checkout next to this repository was once built silently
# (an old EVM branch gave different command-line flags, so predict failed on every partition).
# Set PIXI_INSTALL_STRICT=1 to make a mismatch an error.
report_local_checkout() {
    local name="$1" dir="$2" ref="$3" head want branch dirty
    head=$(git -C "${dir}" rev-parse HEAD 2> /dev/null || echo "unknown")
    branch=$(git -C "${dir}" rev-parse --abbrev-ref HEAD 2> /dev/null || echo "unknown")
    dirty=$(git -C "${dir}" status --porcelain --untracked-files=no 2> /dev/null | grep -vc '^ *[AM]* *plugins/\|target/' || true)
    echo "[${name}] local checkout ${dir}: ${branch} @ ${head:0:12} ($(git -C "${dir}" describe --tags --always 2> /dev/null || echo no-tag)), ${dirty} modified tracked files"
    want=$(git -C "${dir}" rev-parse --verify -q "${ref}^{commit}" 2> /dev/null || git -C "${dir}" rev-parse --verify -q "origin/${ref}^{commit}" 2> /dev/null || true)
    if [ -z "${want}" ]; then
        echo "[${name}] WARNING: the requested ref '${ref}' does not exist in this checkout; building its HEAD (${branch} @ ${head:0:12}) as is." >&2
    elif [ "${want}" != "${head}" ]; then
        echo "[${name}] WARNING: this checkout's HEAD (${head:0:12}) is not the requested ref '${ref}' (${want:0:12})." >&2
    else
        return 0
    fi
    echo "[${name}] Override the ref with the *_RUST_COMMIT variable, or move the checkout. Set PIXI_INSTALL_STRICT=1 to stop on this." >&2
    if [ "${PIXI_INSTALL_STRICT:-0}" = "1" ]; then
        exit 1
    fi
}

if [ -d "${EVM_LOCAL}/evm" ] && [ -f "${EVM_LOCAL}/evm/Cargo.toml" ]; then
    echo "[pixi_install_evm] Using local EVidenceModeler checkout from ${EVM_LOCAL}..."
    report_local_checkout pixi_install_evm "${EVM_LOCAL}" "${EVM_COMMIT}"
    EVM_SOURCE_DIR="${EVM_LOCAL}"
else
    echo "[pixi_install_evm] Local EVidenceModeler not found or missing evm/, cloning from ${EVM_REPO} (commit ${EVM_COMMIT})..."
    EVM_CLONE_DIR=$(mktemp -d)
    trap "rm -rf '${EVM_CLONE_DIR}'" EXIT

    git clone --recursive --depth 1 --jobs=4 -b "${EVM_COMMIT}" "${EVM_REPO}" "${EVM_CLONE_DIR}"
    # Ensure all submodules are initialized and updated (belt-and-suspenders)
    git -C "${EVM_CLONE_DIR}" submodule update --init --recursive 2>/dev/null || true
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
