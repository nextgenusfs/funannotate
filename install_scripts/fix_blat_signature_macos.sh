#!/usr/bin/env bash
# Re-sign blat in the active conda or pixi environment (macOS arm64 only).
#
# The ucsc-blat package for osx-arm64 ships a binary with an invalid code signature ("code or
# signature have been modified"). macOS kills it at launch (exit status 137), and
# `funannotate check` reports blat as not installed. Sign it again with an ad-hoc signature.
#
# The file is copied, signed and renamed over the original. pixi hard-links packages from its
# cache into the environment, so signing the file in place would change the cached copy too, and
# changing a signed file in place gets it killed at launch.
#
# Usage (environment active):
#   bash install_scripts/fix_blat_signature_macos.sh
# Safe to repeat. Does nothing on Linux, on Intel Macs, or when blat is not installed.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set; activate the funannotate environment}"

if [ "$(uname -s)" != "Darwin" ] || [ "$(uname -m)" != "arm64" ]; then
    echo "[fix_blat_signature] not macOS arm64; nothing to do"
    exit 0
fi

BLAT="${CONDA_PREFIX}/bin/blat"
if [ ! -f "${BLAT}" ]; then
    echo "[fix_blat_signature] ${BLAT} not found; nothing to do"
    exit 0
fi

# Use Apple's codesign: an environment can put a different one first on PATH.
CODESIGN=/usr/bin/codesign

if "${CODESIGN}" -v "${BLAT}" 2> /dev/null; then
    echo "[fix_blat_signature] the signature of blat is valid; nothing to do"
    exit 0
fi

TMP="${BLAT}.sign-tmp"
cp "${BLAT}" "${TMP}"
"${CODESIGN}" -f -s - "${TMP}" 2> /dev/null
mv -f "${TMP}" "${BLAT}"

if "${CODESIGN}" -v "${BLAT}" 2> /dev/null; then
    echo "[fix_blat_signature] blat re-signed (ad-hoc signature)"
else
    echo "[fix_blat_signature] blat still has an invalid signature after signing" >&2
    exit 1
fi
