#!/usr/bin/env bash
# Build bowtie2 from source to fix AVX2 runtime fallback issues in conda binary.
#
# The conda-packaged bowtie2 2.5.5 fails to launch x86-64-v3 (AVX2) optimized
# code at runtime, falling back to slower baseline x86-64 performance. Building
# from source with the local toolchain resolves this issue.
#
# This script:
# 1. Clones bowtie2 v2.5.5 from GitHub
# 2. Compiles with the pixi toolchain (cxx-compiler, make, etc.)
# 3. Installs binaries to $CONDA_PREFIX/bin, shadowing the conda package
# 4. Is idempotent: safe to run multiple times

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

BT2_VERSION="2.5.5"
BT2_BIN_DIR="${CONDA_PREFIX}/bin"

# Check if already built
if [ -x "${BT2_BIN_DIR}/bowtie2" ] && [ -x "${BT2_BIN_DIR}/bowtie2-align-l" ]; then
    echo "[pixi_install_bowtie2] Already built, checking for runtime issues..."
    if ! "${BT2_BIN_DIR}/bowtie2" --version 2>&1 | grep -q "WARNING"; then
        echo "[pixi_install_bowtie2] No x86-64-v3 warnings detected, installation OK"
        exit 0
    fi
    echo "[pixi_install_bowtie2] Found AVX2 warnings, rebuilding..."
fi

echo "[pixi_install_bowtie2] Building bowtie2 ${BT2_VERSION} from source..."

# Create temporary build directory
BT2_BUILD_DIR=$(mktemp -d)
trap "rm -rf '${BT2_BUILD_DIR}'" EXIT

cd "${BT2_BUILD_DIR}"

echo "[pixi_install_bowtie2] Cloning bowtie2 v${BT2_VERSION}..."
git clone --depth 1 --branch v${BT2_VERSION} \
    https://github.com/BenLangmead/bowtie2.git bowtie2

cd bowtie2

echo "[pixi_install_bowtie2] Compiling with local toolchain..."
# Use make with parallel jobs (use SLURM_CPUS_ON_NODE if available, otherwise 16)
NPROC=${SLURM_CPUS_ON_NODE:-16}
if ! make -j "${NPROC}" > /tmp/bowtie2_build.log 2>&1; then
    echo "[pixi_install_bowtie2] ERROR: Build failed. Last 50 lines of build log:" >&2
    tail -50 /tmp/bowtie2_build.log >&2
    exit 1
fi

echo "[pixi_install_bowtie2] Build successful, installing binaries..."

# Copy all bowtie2 binaries to shadow the conda package
declare -a BINARIES=(
    "bowtie2"
    "bowtie2-align-l"
    "bowtie2-align-s"
    "bowtie2-build"
    "bowtie2-build-l"
    "bowtie2-build-s"
    "bowtie2-inspect"
    "bowtie2-inspect-l"
    "bowtie2-inspect-s"
)

for bin in "${BINARIES[@]}"; do
    if [ -f "$bin" ]; then
        cp "$bin" "${BT2_BIN_DIR}/" || echo "[pixi_install_bowtie2] WARNING: Failed to copy $bin" >&2
        chmod +x "${BT2_BIN_DIR}/$bin"
    fi
done

echo "[pixi_install_bowtie2] Installation complete. Verifying..."

# Test that new binary works and check for warnings
BT2_OUTPUT=$("${BT2_BIN_DIR}/bowtie2" --version 2>&1)
echo "[pixi_install_bowtie2] Version: $BT2_OUTPUT"

if echo "$BT2_OUTPUT" | grep -q "WARNING"; then
    echo "[pixi_install_bowtie2] WARNING: x86-64-v3 fallback still occurring (may be expected)" >&2
else
    echo "[pixi_install_bowtie2] SUCCESS: No AVX2 fallback warnings detected"
fi

exit 0
