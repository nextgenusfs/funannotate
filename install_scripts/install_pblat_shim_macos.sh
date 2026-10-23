#!/usr/bin/env bash
# Make pblat single-threaded in the active conda environment (macOS arm64 only).
#
# Why: bioconda's pblat (2.5.1) crashes with a bus error (SIGBUS/SIGSEGV) whenever it runs
# two or more threads on Apple Silicon. Several builds were tried: the older arm64 build,
# and builds from source with GCC and with clang all crash too, so it is not a compiler
# problem. The Intel build under Rosetta works. Single-threaded, the arm64 build gives the
# same output as the Intel build (byte-identical in a test). PASA runs pblat as
#   pblat ... -threads=$CPU
# in its blat alignment step, so `funannotate train` and `update` fail at PASA's
# alignment step. That step takes seconds, so one thread costs nothing noticeable.
#
# The script renames the real binary to pblat.real and installs a wrapper named pblat that
# removes the -threads=N argument and runs pblat.real with -threads=1. Linux and Intel
# Macs are not affected and nothing is changed there.
#
# Usage (conda environment active):
#   bash install_scripts/install_pblat_shim_macos.sh
# Run it again after `conda install` or an update replaces pblat.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set; activate the funannotate conda environment}"

if [ "$(uname -s)" != "Darwin" ] || [ "$(uname -m)" != "arm64" ]; then
    echo "[install_pblat_shim] not macOS arm64; nothing to do"
    exit 0
fi

BIN="${CONDA_PREFIX}/bin"
MARK="funannotate pblat shim"

if [ ! -e "${BIN}/pblat" ]; then
    echo "[install_pblat_shim] ${BIN}/pblat not found; install it first (conda install pblat)" >&2
    exit 1
fi

# Already the shim: nothing to do.
if grep -q "${MARK}" "${BIN}/pblat" 2> /dev/null; then
    if [ -x "${BIN}/pblat.real" ]; then
        echo "[install_pblat_shim] the shim is already installed"
        exit 0
    fi
    echo "[install_pblat_shim] the shim is installed but pblat.real is missing; reinstall pblat and run again" >&2
    exit 1
fi

# bin/pblat is a real binary: keep it as pblat.real (replacing an older one), then add the shim.
# Move rather than copy over the old file: macOS kills a binary whose contents change under a
# cached code signature.
mv -f "${BIN}/pblat" "${BIN}/pblat.real"

cat > "${BIN}/pblat.shim-tmp" << 'EOF'
#!/bin/bash
# funannotate pblat shim: pblat crashes with 2 or more threads on Apple Silicon.
# Remove every -threads=N argument and run the real pblat with one thread.
REAL="${REAL_PBLAT:-$(dirname "$0")/pblat.real}"
args=()
for a in "$@"; do
    case "$a" in
        -threads=*) ;;
        *) args+=("$a") ;;
    esac
done
exec "$REAL" "${args[@]}" -threads=1
EOF
chmod +x "${BIN}/pblat.shim-tmp"
mv -f "${BIN}/pblat.shim-tmp" "${BIN}/pblat"

# Check. A signal (exit 128 or more: 137 = killed for a bad signature, 138/139 = bus error or
# segmentation fault) is a failure. An ordinary error on this tiny input is not: only the
# crash with several threads matters here.
rc=0
"${BIN}/pblat.real" > /dev/null 2>&1 || rc=$?
if [ "${rc}" -ge 128 ]; then
    echo "[install_pblat_shim] pblat.real ended with signal status ${rc} (137 means a bad code signature: codesign -f -s - ${BIN}/pblat.real)" >&2
    exit 1
fi
TEST="${CONDA_PREFIX}/pblat_shim_test"
printf '>q\nACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT\n' > "${TEST}.fa"
rc=0
"${BIN}/pblat" "${TEST}.fa" "${TEST}.fa" -q=rna -threads=4 -out=pslx "${TEST}.out" > /dev/null 2>&1 || rc=$?
rm -f "${TEST}.fa" "${TEST}.out"
if [ "${rc}" -ge 128 ]; then
    echo "[install_pblat_shim] pblat still ended with signal status ${rc} with -threads=4 through the shim" >&2
    exit 1
fi
echo "[install_pblat_shim] installed: pblat now runs with one thread; the real binary is pblat.real"
