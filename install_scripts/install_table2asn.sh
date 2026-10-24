#!/usr/bin/env bash
# Install NCBI's table2asn into the active conda or pixi environment (macOS arm64 and Linux x86_64).
#
# bioconda has no table2asn (and no tbl2asn) for osx-arm64, and only the legacy tbl2asn for linux-64. NCBI
# publishes native table2asn builds for both. funannotate prefers table2asn and falls back to tbl2asn
# ($FUN_TBL2ASN overrides).
#
# Usage (environment active):
#   bash install_scripts/install_table2asn.sh
# Environment overrides:
#   TABLE2ASN_URL     download URL (default: NCBI's mac.table2asn.gz or linux64.table2asn.gz)
#   TABLE2ASN_SHA256  expected SHA-256 of the downloaded .gz; the install fails if it differs. NCBI replaces
#                     the file in place, so leave it unset to take the current release; the build log
#                     always records the version and the SHA-256 that was installed.
# Safe to repeat. Does nothing on other platforms (Linux aarch64, Intel Macs), or when table2asn is already
# on PATH. A legacy tbl2asn on PATH does not count: it is the fallback, not the program wanted here.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set; activate the funannotate environment}"

os=$(uname -s)
arch=$(uname -m)
case "${os}-${arch}" in
    Darwin-arm64)
        platform="mac"
        want='Mach-O.*arm64'
        ;;
    Linux-x86_64)
        platform="linux64"
        want='ELF 64-bit.*x86-64'
        ;;
    *)
        echo "[install_table2asn] no NCBI table2asn build for ${os} ${arch}; nothing to do"
        exit 0
        ;;
esac

BIN="${CONDA_PREFIX}/bin"
URL="${TABLE2ASN_URL:-https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/${platform}.table2asn.gz}"

if command -v table2asn > /dev/null 2>&1; then
    echo "[install_table2asn] table2asn is already on PATH: $(command -v table2asn); nothing to do"
    exit 0
fi

sha256() {
    if command -v sha256sum > /dev/null 2>&1; then sha256sum "$1" | cut -d' ' -f1; else shasum -a 256 "$1" | cut -d' ' -f1; fi
}

WORK=$(mktemp -d)
trap 'rm -rf "${WORK}"' EXIT

# curl does not add the macOS quarantine flag; a browser download does.
curl -fsSL --retry 3 -o "${WORK}/table2asn.gz" "${URL}"
digest=$(sha256 "${WORK}/table2asn.gz")
if [ -n "${TABLE2ASN_SHA256:-}" ] && [ "${TABLE2ASN_SHA256}" != "${digest}" ]; then
    echo "[install_table2asn] SHA-256 of ${URL} is ${digest}, expected ${TABLE2ASN_SHA256}" >&2
    exit 1
fi
gunzip "${WORK}/table2asn.gz"
chmod +x "${WORK}/table2asn"

# It must be a native program for this platform that runs.
kind=$(file -b "${WORK}/table2asn")
if ! printf '%s\n' "${kind}" | grep -Eq "${want}"; then
    echo "[install_table2asn] the download is not a native ${os} ${arch} program: ${kind}" >&2
    exit 1
fi
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
echo "[install_table2asn] installed ${version} in ${BIN}/table2asn (download sha256 ${digest})"
