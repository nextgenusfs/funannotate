#!/usr/bin/env bash
# Install Rust-enabled PASA into the active pixi/conda environment.
#
# This script first checks for a local PASApipeline checkout (../PASApipeline,
# rust_optimize branch), then falls back to cloning from GitHub if not available.
# Delegates all build logic to PASApipeline/scripts/install.sh for maintainability.
#
# This script is idempotent: if PASA is already built, it exits immediately.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

PASA_INSTALL_PREFIX="${CONDA_PREFIX}/opt/pasa-rust"
PASA_SRC="${PASA_INSTALL_PREFIX}/src"
PASA_BIN="${PASA_INSTALL_PREFIX}/bin"

# funannotate expects PASA binaries (seqclean, etc.) under $PASAHOME/bin, but
# PASA_rust's install.sh installs them into <prefix>/bin, a sibling of
# <prefix>/src (which is what PASAHOME points at). Symlink so both layouts work.
link_pasa_bin() {
    if [ -d "${PASA_SRC}" ] && [ ! -e "${PASA_SRC}/bin" ]; then
        ln -s ../bin "${PASA_SRC}/bin"
    fi
}

# Check if PASA is already built
if [ -x "${PASA_BIN}/pasa_rust" ] 2>/dev/null && [ -x "${PASA_SRC}/Launch_PASA_pipeline.pl" ] 2>/dev/null; then
    link_pasa_bin
    return 0
fi

# Also check if PASA is already built in a partial state
if [ -f "${PASA_SRC}/Launch_PASA_pipeline.pl" ] 2>/dev/null; then
    echo "[pixi_install_rust_pasa] PASA source found at ${PASA_SRC}, skipping installation"
    link_pasa_bin
    return 0
fi

# Try to use local PASApipeline checkout first
PASA_LOCAL="../PASApipeline"
if [ -d "${PASA_LOCAL}" ] && [ -f "${PASA_LOCAL}/scripts/install.sh" ]; then
    echo "[pixi_install_rust_pasa] Using local PASApipeline from ${PASA_LOCAL}..."
    "${PASA_LOCAL}/scripts/install.sh" --install-prefix "${PASA_INSTALL_PREFIX}"
    link_pasa_bin
    return 0
fi

# Fall back to cloning from GitHub, pinned to a specific commit for
# reproducibility. Override with PASA_RUST_COMMIT to bump deliberately;
# check https://github.com/hyphaltip/PASApipeline/commits/rust_optimize
# for the latest commit before bumping.
PASA_REPO="https://github.com/hyphaltip/PASApipeline"
PASA_COMMIT="${PASA_RUST_COMMIT:-ae40fddd08bfccbbefb33b94ded43f182750eb10}"

echo "[pixi_install_rust_pasa] Local PASApipeline not found, cloning from ${PASA_REPO} (commit ${PASA_COMMIT})..."
mkdir -p "${PASA_INSTALL_PREFIX}"

# Clean up failed partial installation
if [ -d "${PASA_SRC}" ]; then
    echo "[pixi_install_rust_pasa] Cleaning up incomplete installation from ${PASA_SRC}..."
    rm -rf "${PASA_SRC}"
fi

# Clone into a temporary directory first to avoid "same file" cp errors
PASA_TEMP=$(mktemp -d)
trap "rm -rf '${PASA_TEMP}'" EXIT

git init -q "${PASA_TEMP}"
git -C "${PASA_TEMP}" fetch --depth 1 "${PASA_REPO}" "${PASA_COMMIT}"
git -C "${PASA_TEMP}" checkout -q FETCH_HEAD

# Run the install script with temp directory as the source
# This installs PASA properly to PASA_INSTALL_PREFIX/bin and PASA_INSTALL_PREFIX/src
(cd "${PASA_TEMP}" && ./scripts/install.sh --install-prefix "${PASA_INSTALL_PREFIX}")
link_pasa_bin
