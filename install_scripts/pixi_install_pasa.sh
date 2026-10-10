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
PASA_COMMIT="${PASA_RUST_COMMIT:-v2.6.1-rc.3}"

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

PASA_LOCAL="../PASApipeline"
if [ -d "${PASA_LOCAL}" ] && [ -f "${PASA_LOCAL}/scripts/install.sh" ]; then
    echo "[pixi_install] Using local PASApipeline from ${PASA_LOCAL}..."
    report_local_checkout pixi_install "${PASA_LOCAL}" "${PASA_COMMIT}"
    "${PASA_LOCAL}/scripts/install.sh" --install-prefix "${PASA_INSTALL_PREFIX}"
    exit 0
fi

# Fall back to cloning from GitHub. Tracks the rust_optimize branch by default;
# override with PASA_RUST_COMMIT to pin a specific commit for reproducibility;
# see https://github.com/hyphaltip/PASApipeline/commits/rust_optimize
# forcing to a specific tag which is on this branch
PASA_REPO="https://github.com/hyphaltip/PASApipeline"

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
