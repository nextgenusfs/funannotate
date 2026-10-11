#!/usr/bin/env bash
# Build Augustus 3.5.0 with GCC and install it into the active conda environment (macOS only).
#
# Why: `augustus --proteinprofile` crashes (SIGSEGV in PP::SubstateModel::clearLowScoring)
# in every Augustus built with Apple clang, including the bioconda osx-64 and osx-arm64
# packages. The same source built with GCC works. funannotate runs this exact test before
# BUSCO training, so `funannotate predict` cannot start with a clang build. Linux uses GCC
# and is not affected.
#
# Needs Homebrew GCC (`brew install gcc`) and git. The conda augustus package supplies the
# config directory ($CONDA_PREFIX/config), so only the augustus binary is replaced. The
# original binary is kept as augustus.clang.
#
# Usage (conda environment active):
#   bash install_scripts/install_augustus_gcc_macos.sh
#
# Environment overrides: GXX (path of the GCC C++ compiler), AUGUSTUS_TAG (default v3.5.0).

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set; activate the funannotate conda environment}"

if [ "$(uname -s)" != "Darwin" ]; then
    echo "[install_augustus_gcc] not macOS; nothing to do (Linux builds use GCC)"
    exit 0
fi

AUGUSTUS_TAG="${AUGUSTUS_TAG:-v3.5.0}"
BIN="${CONDA_PREFIX}/bin"

# Already done: augustus.clang (the kept conda binary) exists, so this script has run before, and the
# binary in place passes funannotate's --proteinprofile test. If a conda update has replaced augustus
# with a clang build, that test fails and the script rebuilds as usual. Set FORCE_REBUILD=1 to
# rebuild anyway.
if [ -z "${FORCE_REBUILD:-}" ] && [ -e "${BIN}/augustus.clang" ] && python3 - << 'PY' > /dev/null 2>&1
import logging
import funannotate.library as lib

lib.log = logging.getLogger("install_augustus_gcc")
raise SystemExit(0 if lib.checkAugustusFunc()[1] else 1)
PY
then
    echo "[install_augustus_gcc] a working GCC build of augustus is already installed; nothing to do (FORCE_REBUILD=1 rebuilds)"
    exit 0
fi

# Find a Homebrew GCC C++ compiler (g++-13, g++-14, ...). Apple's /usr/bin/g++ is clang.
if [ -z "${GXX:-}" ]; then
    GXX=$(ls /opt/homebrew/bin/g++-[0-9]* /usr/local/bin/g++-[0-9]* 2> /dev/null | sort -V | tail -1 || true)
fi
if [ -z "${GXX}" ] || [ ! -x "${GXX}" ]; then
    echo "[install_augustus_gcc] no Homebrew g++ found. Run: brew install gcc" >&2
    exit 1
fi
if ! "${GXX}" --version 2> /dev/null | grep -qi "gcc\|g++"; then
    echo "[install_augustus_gcc] ${GXX} is not GCC (Apple clang crashes); set GXX to a real GCC" >&2
    exit 1
fi
echo "[install_augustus_gcc] using $("${GXX}" --version | head -1)"

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

# GitHub sometimes drops the connection ("HTTP/2 stream was not closed cleanly"); retry.
for attempt in 1 2 3; do
    rm -rf "${WORK}/Augustus"
    if git clone --quiet --depth 1 --branch "${AUGUSTUS_TAG}" https://github.com/Gaius-Augustus/Augustus.git "${WORK}/Augustus"; then
        break
    fi
    if [ "${attempt}" -eq 3 ]; then
        echo "[install_augustus_gcc] git clone failed 3 times" >&2
        exit 1
    fi
    echo "[install_augustus_gcc] git clone failed (attempt ${attempt}); retrying" >&2
    sleep 5
done

# Build only the augustus binary. Optional features (zipped input, comparative
# gene prediction, MySQL and SQLite) are off, as the conda package does not need them here.
# The Makefile copies the result into ../bin, which must exist.
# The active conda env exports clang-specific CFLAGS/CXXFLAGS/CPPFLAGS/LDFLAGS
# (-isystem $CONDA_PREFIX/include, libc++ options); with GCC they hide its own C++
# headers ("fatal error: cmath: No such file or directory"), so build without them.
# Linking libstdc++/libgcc statically keeps the binary independent of Homebrew's
# GCC runtime (otool -L then lists only libSystem).
mkdir -p "${WORK}/Augustus/bin"
env -u CFLAGS -u CXXFLAGS -u CPPFLAGS -u LDFLAGS \
    make -C "${WORK}/Augustus/src" -j 4 \
    CXX="${GXX}" \
    ZIPINPUT=false COMPGENEPRED=false MYSQL=false SQLITE=false \
    LDFLAGS="-static-libstdc++ -static-libgcc" \
    augustus

NEW="${WORK}/Augustus/bin/augustus"
[ -x "${NEW}" ] || NEW="${WORK}/Augustus/src/augustus"
[ -x "${NEW}" ] || { echo "[install_augustus_gcc] build produced no augustus binary" >&2; exit 1; }

# The GCC runtime is linked statically above; warn if that did not take effect.
if otool -L "${NEW}" | grep -q "libstdc++"; then
    echo "[install_augustus_gcc] note: the binary links Homebrew's libstdc++; keep Homebrew GCC installed"
fi

# Keep the conda binary once, then install the GCC build over it.
if [ ! -e "${BIN}/augustus.clang" ]; then
    cp "${BIN}/augustus" "${BIN}/augustus.clang"
fi
# Copy to a temporary name, then rename over the old binary. Writing into the
# existing file with `cp` reuses its inode: macOS then kills the new binary at launch
# (exit 137) because the kernel still holds the old code signature for that file,
# and if conda hard-linked the file from its package cache, the cache copy would be
# overwritten too. A rename gives the new binary its own inode.
cp "${NEW}" "${BIN}/augustus.gcc-tmp"
chmod +x "${BIN}/augustus.gcc-tmp"
mv -f "${BIN}/augustus.gcc-tmp" "${BIN}/augustus"

# Check with the test that funannotate runs before BUSCO training.
echo "[install_augustus_gcc] installed: $("${BIN}/augustus" --version 2>&1 | head -1)"
python3 - << 'PY' || { echo "[install_augustus_gcc] funannotate's --proteinprofile test FAILED" >&2; exit 1; }
import logging
import funannotate.library as lib
# checkAugustusFunc() logs through lib.log, which only exists after setupLogging()
lib.log = logging.getLogger("install_augustus_gcc")
version, ok = lib.checkAugustusFunc()
print("[install_augustus_gcc] funannotate checkAugustusFunc:", version, ok)
raise SystemExit(0 if ok else 1)
PY
echo "[install_augustus_gcc] done"
