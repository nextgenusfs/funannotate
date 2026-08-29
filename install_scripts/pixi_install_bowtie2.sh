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

# Check if already built. Requires both -v256 AVX2 variants (-s and -l) to
# be present, not just -l -- the conda-forge bowtie2 package can ship one
# without the other (that asymmetry is exactly what let this short-circuit
# false-positive: it saw bowtie2-align-l-v256 from the conda package and
# skipped the rebuild, leaving bowtie2-align-s-v256 missing), so checking
# only -l-v256 isn't enough to conclude AVX2 dispatch is fully wired up.
if [ -x "${BT2_BIN_DIR}/bowtie2" ] && [ -x "${BT2_BIN_DIR}/bowtie2-align-l" ] \
    && [ -x "${BT2_BIN_DIR}/bowtie2-align-l-v256" ] \
    && [ -x "${BT2_BIN_DIR}/bowtie2-align-s" ] \
    && [ -x "${BT2_BIN_DIR}/bowtie2-align-s-v256" ]; then
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

# Copy all bowtie2 binaries to shadow the conda package. Includes the
# -v256 AVX2 (x86-64-v3) variants that `make` builds alongside the SSE2
# baseline: bowtie2-align-{s,l} exec into bowtie2-align-{s,l}-v256 at
# runtime (via cpuid check, falling back gracefully) if it's present
# alongside them -- without these, the baseline-only binaries never use
# AVX2 even though the build produced the accelerated variants.
declare -a BINARIES=(
    "bowtie2"
    "bowtie2-align-l"
    "bowtie2-align-l-v256"
    "bowtie2-align-s"
    "bowtie2-align-s-v256"
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
    echo "[pixi_install_bowtie2] ERROR: x86-64-v3 fallback still occurring after build -- the -v256 AVX2 binaries did not build or install correctly" >&2
    exit 1
fi

# The wrapper's warning check above only covers what --version happens to
# print; also confirm the -v256 binaries bowtie2-align-{s,l} exec into at
# runtime were actually built and installed, not just missing from BINARIES.
if [ ! -x "${BT2_BIN_DIR}/bowtie2-align-s-v256" ] || [ ! -x "${BT2_BIN_DIR}/bowtie2-align-l-v256" ]; then
    echo "[pixi_install_bowtie2] ERROR: bowtie2-align-{s,l}-v256 not installed to ${BT2_BIN_DIR}" >&2
    exit 1
fi

echo "[pixi_install_bowtie2] SUCCESS: No AVX2 fallback warnings detected"
exit 0
