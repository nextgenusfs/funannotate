#!/usr/bin/env bash
# Install Rust-enabled PASA into the active pixi/conda environment.
#
# This script first checks for a local PASA_rust checkout (../PASA_rust), then
# falls back to cloning from GitHub if not available. Delegates all build logic
# to PASA_rust/scripts/install.sh for maintainability.
#
# This script is idempotent: if PASA is already built, it exits immediately.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

PASA_INSTALL_PREFIX="${CONDA_PREFIX}/opt/pasa-rust-3.0"
PASA_SRC="${PASA_INSTALL_PREFIX}/src"
PASA_BIN="${PASA_INSTALL_PREFIX}/bin"

# Check if PASA is already built
if [ -x "${PASA_BIN}/pasa_rust" ] && [ -x "${PASA_SRC}/Launch_PASA_pipeline.pl" ]; then
    exit 0
fi

# Try to use local PASA_rust checkout first
PASA_LOCAL="../PASA_rust"
if [ -d "${PASA_LOCAL}" ] && [ -f "${PASA_LOCAL}/scripts/install.sh" ]; then
    echo "[pixi_install_rust_pasa] Using local PASA_rust from ${PASA_LOCAL}..."
    "${PASA_LOCAL}/scripts/install.sh" --install-prefix "${PASA_INSTALL_PREFIX}"
    exit 0
fi

# Fall back to cloning from GitHub
PASA_REPO="https://github.com/hyphaltip/PASA_rust"
PASA_BRANCH="main"

echo "[pixi_install_rust_pasa] Local PASA_rust not found, cloning from ${PASA_REPO} (${PASA_BRANCH} branch)..."
mkdir -p "${PASA_INSTALL_PREFIX}"

# Clean up failed partial installation
if [ -d "${PASA_SRC}" ]; then
    echo "[pixi_install_rust_pasa] Cleaning up incomplete installation from ${PASA_SRC}..."
    rm -rf "${PASA_SRC}"
fi

git clone --branch "${PASA_BRANCH}" "${PASA_REPO}" "${PASA_SRC}"

# Run the install script from the cloned repo
"${PASA_SRC}/scripts/install.sh" --install-prefix "${PASA_INSTALL_PREFIX}"
