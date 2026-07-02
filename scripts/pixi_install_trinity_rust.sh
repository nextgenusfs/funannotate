#!/usr/bin/env bash
# Install Trinity from source into /opt/trinity_rust with Rust-optimized utilities.
#
# This script:
# 1. Clones Trinity (with submodules) into /opt/trinity_rust
# 2. Builds Trinity via make
# 3. Runs install.py to symlink the Trinity executable and install Rust bio-utilities
#    into /env/default/bin
#
# The install.py script handles:
#   - Symlinking the Trinity main script to /env/default/bin
#   - Building and installing Rust binaries (sam_to_read_coords, extract_reads_per_partition,
#     fragment_coverage_writer, define_coverage_partitions)
#   - Replacing slower Perl versions with Rust equivalents
#
# This script is idempotent: if Trinity is already installed at /opt/trinity_rust, it exits.

set -euo pipefail

TRINITY_INSTALL_DIR="/opt/trinity_rust"
TRINITY_REPO="https://github.com/hyphaltip/trinityrnaseq"
TRINITY_COMMIT="${TRINITY_RUST_COMMIT:-7f04a37c16cbbcf2e694695a34de01d7941b313f}"
ENV_BIN_DIR="/env/default/bin"

# Skip if already installed
if [ -d "${TRINITY_INSTALL_DIR}" ] && [ -f "${TRINITY_INSTALL_DIR}/install.py" ]; then
    echo "[pixi_install_trinity_rust] Trinity already installed at ${TRINITY_INSTALL_DIR}"
    # Ensure symlinks are in place
    if [ -x "${ENV_BIN_DIR}/Trinity" ]; then
        echo "[pixi_install_trinity_rust] Trinity executable already symlinked"
        return 0
    fi
fi

echo "[pixi_install_trinity_rust] Cloning Trinity from ${TRINITY_REPO} (commit ${TRINITY_COMMIT})..."

# Create /opt if it doesn't exist
mkdir -p /opt

# Clone Trinity with submodules
if [ ! -d "${TRINITY_INSTALL_DIR}" ]; then
    git clone --recursive "${TRINITY_REPO}" "${TRINITY_INSTALL_DIR}"
    git -C "${TRINITY_INSTALL_DIR}" checkout "${TRINITY_COMMIT}"
else
    echo "[pixi_install_trinity_rust] Trinity directory already exists at ${TRINITY_INSTALL_DIR}"
fi

echo "[pixi_install_trinity_rust] Building Trinity with make..."
make -C "${TRINITY_INSTALL_DIR}"

echo "[pixi_install_trinity_rust] Running install.py to set up symlinks and install Rust utilities..."
# install.py handles all binaries and rust_bio tools setup
python3 "${TRINITY_INSTALL_DIR}/install.py" "${ENV_BIN_DIR}"

echo "[pixi_install_trinity_rust] Trinity installation complete at ${TRINITY_INSTALL_DIR}"
echo "[pixi_install_trinity_rust] Executables symlinked to ${ENV_BIN_DIR}"
