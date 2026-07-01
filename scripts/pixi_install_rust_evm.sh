#!/usr/bin/env bash
# Install Rust EVidenceModeler into the active pixi/conda environment.
#
# This script first checks for a local EVidenceModeler_rust checkout (../EVidenceModeler_rust),
# then falls back to cloning from GitHub if not available. Delegates all build logic
# to EVidenceModeler_rust/scripts/install.sh for maintainability.
#
# This script is idempotent: if the Rust EVM binaries are already present, it exits immediately.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

EVM_INSTALL_PREFIX="${CONDA_PREFIX}/opt/evm-rust"
BIN_DIR="${EVM_INSTALL_PREFIX}/bin"

# Skip if already installed
if [ -x "${BIN_DIR}/evidence_modeler" ]; then
    return 0
fi

# Try to use local EVidenceModeler_rust checkout first
EVM_LOCAL="../EVidenceModeler_rust"
if [ -d "${EVM_LOCAL}" ] && [ -f "${EVM_LOCAL}/scripts/install.sh" ]; then
    echo "[pixi_install_rust_evm] Using local EVidenceModeler_rust from ${EVM_LOCAL}..."
    "${EVM_LOCAL}/scripts/install.sh" --install-prefix "${EVM_INSTALL_PREFIX}"
    return 0
fi

# Fall back to cloning from GitHub
EVM_REPO="https://github.com/hyphaltip/EVidenceModeler_rust.git"
EVM_BRANCH="main"
SRC_DIR="${EVM_INSTALL_PREFIX}/src"

echo "[pixi_install_rust_evm] Local EVidenceModeler_rust not found, cloning from ${EVM_REPO} (${EVM_BRANCH} branch)..."
mkdir -p "${EVM_INSTALL_PREFIX}"

# Clean up failed partial installation
if [ -d "${SRC_DIR}" ]; then
    echo "[pixi_install_rust_evm] Cleaning up incomplete installation from ${SRC_DIR}..."
    rm -rf "${SRC_DIR}"
fi

git clone --depth 1 --branch "${EVM_BRANCH}" "${EVM_REPO}" "${SRC_DIR}"

# Run the install script from the cloned repo
"${SRC_DIR}/scripts/install.sh" --install-prefix "${EVM_INSTALL_PREFIX}"
