#!/usr/bin/env bash
# Build and install Rust-optimized Trinity utilities into the active pixi/conda environment.
#
# This script builds the Rust Trinity bio-utilities (sam_to_read_coords, extract_reads_per_partition,
# fragment_coverage_writer, define_coverage_partitions) from a Trinity checkout: it first checks for
# a local checkout (../trinityrnaseq), then falls back to cloning from GitHub if not available.
# These Rust binaries replace the slower Perl versions of:
#   - util/support_scripts/extract_reads_per_partition.pl
#   - util/support_scripts/fragment_coverage_writer.pl
#   - util/support_scripts/define_coverage_partitions.pl
#   - SAM-to-coordinate conversion utilities
#
# The Trinity main script and Perl utilities remain available from the conda Trinity package.
#
# This script is idempotent: if the Rust binaries are already present, it exits immediately.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

TRINITY_RUST_BIN_DIR="${CONDA_PREFIX}/bin"

# Skip if already installed
if [ -x "${TRINITY_RUST_BIN_DIR}/sam_to_read_coords" ] && \
   [ -x "${TRINITY_RUST_BIN_DIR}/extract_reads_per_partition" ] && \
   [ -x "${TRINITY_RUST_BIN_DIR}/fragment_coverage_writer" ] && \
   [ -x "${TRINITY_RUST_BIN_DIR}/define_coverage_partitions" ]; then
    echo "[pixi_install_trinity_rust] Rust Trinity utilities already installed"
    return 0
fi

# Try to use local Trinity checkout first
TRINITY_LOCAL="../trinityrnaseq"
TRINITY_REPO="https://github.com/hyphaltip/trinityrnaseq"
TRINITY_BRANCH="devel"
TRINITY_SRC=""
TRINITY_CLONE_DIR=""

if [ -d "${TRINITY_LOCAL}/rust_bio_utils" ] && [ -f "${TRINITY_LOCAL}/rust_bio_utils/Cargo.toml" ]; then
    echo "[pixi_install_trinity_rust] Using local Trinity checkout from ${TRINITY_LOCAL}..."
    TRINITY_SRC="${TRINITY_LOCAL}"
else
    echo "[pixi_install_trinity_rust] Local Trinity checkout not found or missing rust_bio_utils, cloning from ${TRINITY_REPO} (${TRINITY_BRANCH} branch)..."
    TRINITY_CLONE_DIR=$(mktemp -d)
    trap "rm -rf '${TRINITY_CLONE_DIR}'" EXIT

    git clone --depth 1 --branch "${TRINITY_BRANCH}" "${TRINITY_REPO}" "${TRINITY_CLONE_DIR}"
    TRINITY_SRC="${TRINITY_CLONE_DIR}"

    if [ ! -d "${TRINITY_SRC}/rust_bio_utils" ]; then
        echo "[pixi_install_trinity_rust] ERROR: rust_bio_utils not found in cloned ${TRINITY_REPO}@${TRINITY_BRANCH}" >&2
        return 1
    fi

    if [ ! -f "${TRINITY_SRC}/rust_bio_utils/Cargo.toml" ]; then
        echo "[pixi_install_trinity_rust] ERROR: Cargo.toml not found in rust_bio_utils" >&2
        return 1
    fi
fi

echo "[pixi_install_trinity_rust] Building Rust Trinity utilities from ${TRINITY_SRC}..."

# Save current environment in case Trinity's pixi changes it
SAVED_CONDA_PREFIX="$CONDA_PREFIX"
SAVED_PATH="$PATH"

# Build in a subshell to isolate directory changes from affecting parent environment
(
  cd "${TRINITY_SRC}/rust_bio_utils"

  # Clean any previous builds to avoid stale artifacts
  echo "[pixi_install_trinity_rust] Cleaning previous build artifacts..."
  cargo clean 2>/dev/null || true

  # Build the Rust binaries in release mode
  cargo build --release
)

# Restore original environment
export CONDA_PREFIX="$SAVED_CONDA_PREFIX"
export PATH="$SAVED_PATH"

# Install binaries to conda prefix
RELEASE_DIR="${TRINITY_SRC}/rust_bio_utils/target/release"
for binary in sam_to_read_coords extract_reads_per_partition fragment_coverage_writer define_coverage_partitions; do
    if [ -x "${RELEASE_DIR}/${binary}" ]; then
        cp "${RELEASE_DIR}/${binary}" "${TRINITY_RUST_BIN_DIR}/"
        echo "[pixi_install_trinity_rust] Installed ${binary}"
    else
        echo "[pixi_install_trinity_rust] WARNING: ${binary} not found after build"
    fi
done

# Optional: install the shared library
if [ -f "${RELEASE_DIR}/libtrinity_bio.so" ]; then
    cp "${RELEASE_DIR}/libtrinity_bio.so" "${CONDA_PREFIX}/lib/"
    echo "[pixi_install_trinity_rust] Installed libtrinity_bio.so"
fi

# Clean up the temporary clone, if one was made
if [ -n "${TRINITY_CLONE_DIR}" ]; then
    rm -rf "${TRINITY_CLONE_DIR}"
    trap - EXIT
fi

echo "[pixi_install_trinity_rust] Done. Rust Trinity utilities are now in PATH."
