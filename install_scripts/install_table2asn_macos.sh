#!/usr/bin/env bash
# Install NCBI's table2asn into the active conda or pixi environment (macOS arm64 only).
#
# bioconda has no table2asn (and no tbl2asn) for osx-arm64, but NCBI's own macOS build is a native
# arm64 binary. funannotate prefers table2asn and falls back to tbl2asn ($FUN_TBL2ASN overrides).
#
# Usage (environment active):
#   bash install_scripts/install_table2asn_macos.sh
# Environment override: TABLE2ASN_URL (default: NCBI's mac.table2asn.gz). Safe to repeat. Does nothing
# on Linux or Intel Macs, or when a working table2asn or tbl2asn is already on PATH.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set; activate the funannotate environment}"

if [ "$(uname -s)" != "Darwin" ] || [ "$(uname -m)" != "arm64" ]; then
    echo "[install_table2asn] not macOS arm64; nothing to do"
    exit 0
fi

BIN="${CONDA_PREFIX}/bin"
URL="${TABLE2ASN_URL:-https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/mac.table2asn.gz}"

for tool in table2asn tbl2asn; do
    if command -v "${tool}" > /dev/null 2>&1; then
        echo "[install_table2asn] ${tool} is already on PATH: $(command -v "${tool}"); nothing to do"
        exit 0
    fi
done

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

# curl does not add the macOS quarantine flag; a browser download does.
curl -fsSL --retry 3 -o "${WORK}/table2asn.gz" "${URL}"
gunzip "${WORK}/table2asn.gz"
chmod +x "${WORK}/table2asn"

# It must be a native arm64 Mach-O program that runs.
kind=$(file -b "${WORK}/table2asn")
case "${kind}" in
    *Mach-O*arm64*) ;;
    *)
        echo "[install_table2asn] the download is not an arm64 Mach-O program: ${kind}" >&2
        exit 1
        ;;
esac
version=$("${WORK}/table2asn" -version 2>&1 | head -1 || true)
case "${version}" in
    table2asn:*) ;;
    *)
        echo "[install_table2asn] the downloaded program did not report a version: ${version}" >&2
        exit 1
        ;;
esac

# Copy to a temporary name, then rename, so a program that is running is never rewritten in place.
cp "${WORK}/table2asn" "${BIN}/table2asn.tmp"
mv -f "${BIN}/table2asn.tmp" "${BIN}/table2asn"
echo "[install_table2asn] installed ${version} in ${BIN}/table2asn"
