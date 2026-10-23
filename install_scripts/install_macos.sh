#!/usr/bin/env bash
# Install funannotate on an Apple Silicon Mac: clone the release, build the pixi environment from
# the lock and the source-built externals (EVM, PASA, Trinity, the macOS fixes), then run the smoke checks.
#
# A release attaches this script as install-funannotate-macos.sh with FUNANNOTATE_REF set to the
# release tag. It builds on your machine, so nothing is relocated or unsigned. Takes about 30 minutes.
#
#   bash install-funannotate-macos.sh              # clone into ./funannotate-<ref>
#   DEST=~/funannotate bash install-funannotate-macos.sh
#
# Then:  cd <dest> && pixi shell      (or:  pixi run funannotate --help)
#
# Needs Homebrew (for GCC: the Augustus build) and git. pixi is installed if it is not on PATH.

set -euo pipefail

FUNANNOTATE_REF="${FUNANNOTATE_REF:-@FUNANNOTATE_REF@}"
FUNANNOTATE_REPO="${FUNANNOTATE_REPO:-https://github.com/nextgenusfs/funannotate.git}"
PIXI_VERSION="${PIXI_VERSION:-v0.81.0}"
DEST="${DEST:-${PWD}/funannotate-${FUNANNOTATE_REF#v}}"

die() { echo "ERROR: $*" >&2; exit 1; }
step() { echo; echo "==> $*"; }

case "${FUNANNOTATE_REF}" in
    "" | @*) die "FUNANNOTATE_REF is not set. Pass the release tag: FUNANNOTATE_REF=v1.9.1 bash $0" ;;
esac
[ "$(uname -s)" = "Darwin" ] && [ "$(uname -m)" = "arm64" ] \
    || die "this script is for Apple Silicon Macs (found $(uname -s) $(uname -m)); use the container or conda on Linux"
command -v git > /dev/null || die "git is not installed (xcode-select --install)"
command -v curl > /dev/null || die "curl is not installed"
command -v brew > /dev/null || die "Homebrew is not installed (https://brew.sh); it provides the GCC used to build Augustus"
[ ! -e "${DEST}" ] || die "${DEST} already exists; remove it or set DEST to another path"

step "pixi ${PIXI_VERSION}"
if ! command -v pixi > /dev/null; then
    curl -fsSL https://pixi.sh/install.sh | PIXI_VERSION="${PIXI_VERSION}" PIXI_NO_PATH_UPDATE=1 bash
    export PATH="${HOME}/.pixi/bin:${PATH}"
fi
pixi --version

step "Homebrew GCC"
brew list gcc > /dev/null 2>&1 || brew install gcc

step "Cloning funannotate ${FUNANNOTATE_REF} into ${DEST}"
git clone --depth 1 --branch "${FUNANNOTATE_REF}" "${FUNANNOTATE_REPO}" "${DEST}"
cd "${DEST}"

step "Installing the environment from the lock"
pixi install --locked

step "Building the externals and applying the macOS fixes"
pixi run install-externals

step "Smoke checks"
pixi run macos-smoke-checks

cat << DONE

funannotate ${FUNANNOTATE_REF} is installed in ${DEST}.
  cd ${DEST} && pixi shell        # activates the environment (sets PASAHOME, EVM_HOME, TRINITY_HOME, ...)
  pixi run funannotate check --show-versions
You still need a database: funannotate setup -d <dir>, then export FUNANNOTATE_DB=<dir>.
DONE
