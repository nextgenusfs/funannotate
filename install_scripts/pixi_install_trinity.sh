#!/usr/bin/env bash
# Install Trinity from source into $CONDA_PREFIX/opt/trinityrnaseq with
# Rust-optimized utilities.
#
# This script:
# 1. Clones Trinity (with submodules) into $CONDA_PREFIX/opt/trinityrnaseq
# 2. Builds Trinity via make (C++ Inchworm/Chrysalis + the Rust bio-utilities)
# 3. Runs install.py to symlink the main Trinity executable into $CONDA_PREFIX/bin
# 4. Symlinks the built Rust bio-utilities (sam_to_read_coords,
#    extract_reads_per_partition, fragment_coverage_writer,
#    define_coverage_partitions, ...) into $CONDA_PREFIX/bin so Trinity picks the
#    Rust versions over the slower Perl fallbacks. NOTE: this fork's install.py
#    only links the Trinity script, so step 4 is done here, not by install.py.
#
# This script is idempotent: if Trinity is already installed, it exits.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

TRINITY_INSTALL_DIR="${CONDA_PREFIX}/opt/trinityrnaseq"
TRINITY_REPO="https://github.com/hyphaltip/trinityrnaseq"
# Tracks the rust_optimize branch by default. Override with TRINITY_RUST_COMMIT
# to pin a specific commit for reproducibility; see
# https://github.com/hyphaltip/trinityrnaseq/commits/rust_optimize
TRINITY_COMMIT="${TRINITY_RUST_COMMIT:-rust_optimize}"
ENV_BIN_DIR="${CONDA_PREFIX}/bin"

RUST_UTILS_DIR="${TRINITY_INSTALL_DIR}/rust_bio_utils/target/release"

# install.py only symlinks the main Trinity script. The Rust bio-utilities must
# also be on PATH: Trinity resolves them by bare name and falls back to the
# slower Perl versions when they are absent (see toggle_trinity_rust.sh). Symlink
# every built release binary (skipping cargo's .d dep files) into the env bin.
link_rust_utils() {
    [ -d "${RUST_UTILS_DIR}" ] || return 0
    for util in "${RUST_UTILS_DIR}"/*; do
        [ -f "${util}" ] && [ -x "${util}" ] || continue
        ln -sf "${util}" "${ENV_BIN_DIR}/$(basename "${util}")"
    done
}

# Skip if already installed
if [ -d "${TRINITY_INSTALL_DIR}" ] && [ -f "${TRINITY_INSTALL_DIR}/install.py" ]; then
    echo "[pixi_install_trinity] Trinity already installed at ${TRINITY_INSTALL_DIR}"
    # Ensure symlinks are in place
    if [ -x "${ENV_BIN_DIR}/Trinity" ] && [ -x "${ENV_BIN_DIR}/sam_to_read_coords" ]; then
        echo "[pixi_install_trinity] Trinity executable and Rust utils already symlinked"
        link_rust_utils
        exit 0
    fi
fi

echo "[pixi_install_trinity] Cloning Trinity from ${TRINITY_REPO} (commit ${TRINITY_COMMIT})..."

# Create the install parent dir if it doesn't exist
mkdir -p "${CONDA_PREFIX}/opt"

# Clone Trinity with submodules
if [ ! -d "${TRINITY_INSTALL_DIR}" ]; then
    git clone --recursive "${TRINITY_REPO}" "${TRINITY_INSTALL_DIR}"
    git -C "${TRINITY_INSTALL_DIR}" checkout "${TRINITY_COMMIT}"
else
    echo "[pixi_install_trinity] Trinity directory already exists at ${TRINITY_INSTALL_DIR}"
fi

echo "[pixi_install_trinity] Building Trinity with make..."
# Trinity's bundled Inchworm/Chrysalis declare cmake_minimum_required(VERSION 3.1),
# which CMake >= 4 rejects ("Compatibility with CMake < 3.5 has been removed").
# CMake >= 3.31 honors this env var as the floor policy version, letting the old
# CMakeLists configure without patching upstream. Remove if Trinity bumps its
# cmake_minimum_required past 3.5.
export CMAKE_POLICY_VERSION_MINIMUM="${CMAKE_POLICY_VERSION_MINIMUM:-3.5}"
make -C "${TRINITY_INSTALL_DIR}"

echo "[pixi_install_trinity] Running install.py to set up symlinks and install Rust utilities..."
# install.py handles all binaries and rust_bio tools setup. Its CLI takes an
# 'install' action plus --install-dir (default ~/.local/bin), not a bare path.
python3 "${TRINITY_INSTALL_DIR}/install.py" install --install-dir "${ENV_BIN_DIR}"

# install.py only links the main Trinity script; put the Rust utils on PATH too.
link_rust_utils

echo "[pixi_install_trinity] Trinity installation complete at ${TRINITY_INSTALL_DIR}"
echo "[pixi_install_trinity] Executables symlinked to ${ENV_BIN_DIR}"
