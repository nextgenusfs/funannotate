#!/bin/bash -l
# Install Salmon if not already available.
# Salmon is specified in pixi.toml (salmon = ">=0.9") and should be installed
# via conda. This script only downloads manually if conda version not found.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

# Check if salmon is already available in PATH (from conda)
if command -v salmon &> /dev/null; then
    echo "[pixi_install_salmon] Salmon already available in PATH:"
    salmon --version 2>&1 | head -1
    exit 0
fi

# Check if salmon was previously installed manually
SALMON_MANUAL="${CONDA_PREFIX}/opt/salmon-1.10.0/salmon-latest_linux_x86_64/bin/salmon"
if [ -x "${SALMON_MANUAL}" ]; then
    echo "[pixi_install_salmon] Found existing manual Salmon installation"
    mkdir -p "${CONDA_PREFIX}/bin"
    rm -f "${CONDA_PREFIX}/bin/salmon"
    ln -s ../opt/salmon-1.10.0/salmon-latest_linux_x86_64/bin/salmon "${CONDA_PREFIX}/bin/salmon"
    exit 0
fi

# Fall back to manual download if conda version not available
echo "[pixi_install_salmon] Salmon not found, downloading v1.10.0..."
mkdir -p src
cd src

# Download Salmon
wget -q https://github.com/COMBINE-lab/salmon/releases/download/v1.10.0/salmon-1.10.0_linux_x86_64.tar.gz \
    -O salmon-1.10.0_linux_x86_64.tar.gz

# Extract
tar xf salmon-1.10.0_linux_x86_64.tar.gz
cd ..

# Install to conda prefix
rsync -a src/salmon-latest_linux_x86_64 "${CONDA_PREFIX}/opt/salmon-1.10.0"

# Create symlink (rsync creates a subdirectory, so adjust path accordingly)
mkdir -p "${CONDA_PREFIX}/bin"
rm -f "${CONDA_PREFIX}/bin/salmon"
ln -s ../opt/salmon-1.10.0/salmon-latest_linux_x86_64/bin/salmon "${CONDA_PREFIX}/bin/salmon"

echo "[pixi_install_salmon] Salmon v1.10.0 installed"
