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

PASA_INSTALL_PREFIX="${CONDA_PREFIX}/opt/pasa"
PASA_SRC="${PASA_INSTALL_PREFIX}/src"
PASA_BIN="${PASA_INSTALL_PREFIX}/bin"

# The src/bin -> ../bin symlink that makes $PASAHOME ($PASA_SRC) self-contained is
# created by PASApipeline/scripts/install.sh itself, so this wrapper never needs to
# touch the layout -- it just delegates the build.

# Check if PASA is already fully built. Require all four rust binaries under the
# names PASA's PerlLib probes for (cdbyank_rust / faidx_rust, not the doubled
# cdbyank_rust_rust), the launcher, and the src/bin symlink; otherwise fall
# through and let install.sh repair a stale/misnamed install in place.
if [ -x "${PASA_BIN}/pasa_rust" ] && [ -x "${PASA_BIN}/slclust_rust" ] \
    && [ -x "${PASA_BIN}/cdbyank_rust" ] && [ -x "${PASA_BIN}/faidx_rust" ] \
    && [ -x "${PASA_SRC}/Launch_PASA_pipeline.pl" ] && [ -e "${PASA_SRC}/bin" ]; then
    exit 0
fi

# Try to use local PASApipeline checkout first
PASA_LOCAL="../PASApipeline"
if [ -d "${PASA_LOCAL}" ] && [ -f "${PASA_LOCAL}/scripts/install.sh" ]; then
    echo "[pixi_install] Using local PASApipeline from ${PASA_LOCAL}..."
    "${PASA_LOCAL}/scripts/install.sh" --install-prefix "${PASA_INSTALL_PREFIX}"
    exit 0
fi

# Fall back to cloning from GitHub. Tracks the rust_optimize branch by default;
# override with PASA_RUST_COMMIT to pin a specific commit for reproducibility;
# see https://github.com/hyphaltip/PASApipeline/commits/rust_optimize
PASA_REPO="https://github.com/hyphaltip/PASApipeline"
PASA_COMMIT="${PASA_RUST_COMMIT:-rust_optimize}"

echo "[pixi_install] Local PASApipeline not found, cloning from ${PASA_REPO} (commit ${PASA_COMMIT})..."
mkdir -p "${PASA_INSTALL_PREFIX}"

# Clean up failed partial installation
if [ -d "${PASA_SRC}" ]; then
    echo "[pixi_install] Cleaning up incomplete installation from ${PASA_SRC}..."
    rm -rf "${PASA_SRC}"
fi

# Clone into a temporary directory first to avoid "same file" cp errors
PASA_TEMP=$(mktemp -d)
trap "rm -rf '${PASA_TEMP}'" EXIT

git clone --recursive --depth 1 --jobs=4 -b "${PASA_COMMIT}" "${PASA_REPO}" "${PASA_TEMP}"

# Ensure all submodules are initialized and updated (belt-and-suspenders)
git -C "${PASA_TEMP}" submodule update --init --recursive 2>/dev/null || true

# Run the install script with temp directory as the source
# This installs PASA properly to PASA_INSTALL_PREFIX/bin and PASA_INSTALL_PREFIX/src
# (and links src/bin -> ../bin so PASAHOME=$PASA_SRC is self-contained).
(cd "${PASA_TEMP}" && ./scripts/install.sh --install-prefix "${PASA_INSTALL_PREFIX}")
