#!/usr/bin/env bash
# Install Trinity from source into $CONDA_PREFIX/opt/trinityrnaseq with
# Rust-optimized utilities.
#
# This script:
# 1. Clones Trinity (with submodules) into $CONDA_PREFIX/opt/trinityrnaseq
# 2. Builds Trinity via make (C++ Inchworm/Chrysalis + the Rust bio-utilities)
# 3. Runs install.py to symlink the main Trinity executable into $CONDA_PREFIX/bin
# 4. Symlinks the built Rust bio-utilities (sam_to_read_coords,
#    extract_reads_per_partition, fragment_coverage_writer,
#    define_coverage_partitions, ...) into $CONDA_PREFIX/bin so Trinity picks the
#    Rust versions over the slower Perl fallbacks. NOTE: this fork's install.py
#    only links the Trinity script, so step 4 is done here, not by install.py.
#
# This script is idempotent: if Trinity is already installed, it exits.

set -euo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set — run inside a pixi/conda env}"

TRINITY_INSTALL_DIR="${CONDA_PREFIX}/opt/trinityrnaseq"
TRINITY_REPO="https://github.com/hyphaltip/trinityrnaseq"
# Pinned, not floating: the branch previously tracked here (rust_optimize)
# moved 11 commits between two builds of this image without anyone touching
# this script, and one of those commits (5b6a304, "Replace per-run dynamic
# Butterfly CDS archive with a shipped static one") is exactly what broke
# Butterfly's CDS archive in production -- see the neutralize_broken_butterfly_cds
# comment below. Override with TRINITY_RUST_COMMIT for testing a newer branch
# head, but bumping the default here should be a deliberate, reviewed change;
# see https://github.com/hyphaltip/trinityrnaseq/commits/rust_optimize
TRINITY_COMMIT="${TRINITY_RUST_COMMIT:-5b6a304}"
ENV_BIN_DIR="${CONDA_PREFIX}/bin"

RUST_UTILS_DIR="${TRINITY_INSTALL_DIR}/rust_bio_utils/target/release"

# Delete the pre-built, version-mismatched Butterfly CDS archive that commit
# 5b6a304 checks into the repo at Butterfly/butterfly_cds.jsa (13MB, plus its
# companion Butterfly/butterfly_cds.classlist).
#
# Background: Butterfly.jar launches one JVM per assembled graph component --
# up to hundreds of thousands for a full de-novo assembly -- and Trinity's own
# driver (maybe_apply_bfly_cds() in the `Trinity` script) auto-enables
# `-Xshare:auto -XX:SharedArchiveFile=.../butterfly_cds.jsa` on every one of
# them whenever that file exists. Two independent problems make the checked-in
# file actively harmful rather than merely stale:
#
# 1. Format version mismatch: the archive was dumped with whoever's JDK built
#    it (CDS version 0x12, confirmed by inspecting the file's header), which
#    does not match this image's JDK (openjdk 25.0.1, CDS version 0x13).
#    -Xshare:auto doesn't crash on a mismatch, it just logs
#      [cds] The shared archive file version 0x12 does not match the required version 0x13.
#      [cds] Loading static archive failed.
#    and falls back to running without any archive at all -- including the
#    JVM's own *default* class-data archive, because -XX:SharedArchiveFile
#    replaces rather than layers on top of it. Observed impact: one assembly
#    logged this warning 23,775+ times; another sat 19+ hours at 84% through
#    Butterfly with a 420MB log that was almost entirely this warning.
#    Benchmarked impact of the mismatched archive vs. no archive at all: ~36%
#    slower per JVM launch (285ms -> 387ms on a trivial test graph).
#
# 2. Even a version-matched archive cannot load for this specific jar: its
#    manifest carries `Class-Path: .` (`unzip -p Butterfly.jar
#    META-INF/MANIFEST.MF`), which CDS records as an app classpath entry and
#    then refuses to reuse the archive against ("shared class paths mismatch"
#    / "directory is not empty") because a directory (as opposed to a jar) on
#    the classpath can't be validated as unchanged between dump and load.
#    Confirmed by dumping a fresh archive with this exact image's java (both
#    a dynamic -XX:ArchiveClassesAtExit archive and a correctly-versioned
#    static -Xshare:dump archive built via upstream's own
#    build_butterfly_cds_archive.sh) and reloading it: both fail the same way.
#    So regenerating the archive at build time (an earlier version of this
#    fix attempted exactly that) does not help -- CDS cannot work for
#    Butterfly.jar at all while its manifest carries `Class-Path: .`, unless
#    that manifest entry is stripped upstream in the fork first.
#
# Given (2), regenerating is not a viable build-time fix here. Deleting the
# file is: with it gone, maybe_apply_bfly_cds() has nothing to auto-enable,
# so Butterfly just runs java normally and gets the JVM's own default
# class-data sharing -- verified to be faster than the broken custom-archive
# path (see benchmarks above), at zero risk of a future JDK/archive mismatch.
# If the fork ever fixes the Class-Path manifest issue and ships a working
# archive, this step (and this whole comment) should be removed instead of
# adding a second, competing CDS mechanism next to it.
neutralize_broken_butterfly_cds() {
    local bfly_dir="${TRINITY_INSTALL_DIR}/Butterfly"
    local jsa="${bfly_dir}/butterfly_cds.jsa"
    local classlist="${bfly_dir}/butterfly_cds.classlist"

    if [ -e "${jsa}" ] || [ -e "${classlist}" ]; then
        echo "[pixi_install_trinity] Removing known-broken Butterfly CDS archive (${jsa}) -- see comment in this script for why"
        rm -f "${jsa}" "${classlist}"
    fi
}

# install.py only symlinks the main Trinity script. The Rust bio-utilities must
# also be on PATH: Trinity resolves them by bare name and falls back to the
# slower Perl versions when they are absent (see toggle_trinity_rust.sh). Symlink
# every built release binary (skipping cargo's .d dep files) into the env bin.
link_rust_utils() {
    [ -d "${RUST_UTILS_DIR}" ] || return 0
    for util in "${RUST_UTILS_DIR}"/*; do
        [ -f "${util}" ] && [ -x "${util}" ] || continue
        ln -sf "${util}" "${ENV_BIN_DIR}/$(basename "${util}")"
    done
}

# Skip if already installed
if [ -d "${TRINITY_INSTALL_DIR}" ] && [ -f "${TRINITY_INSTALL_DIR}/install.py" ]; then
    echo "[pixi_install_trinity] Trinity already installed at ${TRINITY_INSTALL_DIR}"
    # Ensure symlinks are in place
    if [ -x "${ENV_BIN_DIR}/Trinity" ] && [ -x "${ENV_BIN_DIR}/sam_to_read_coords" ]; then
        echo "[pixi_install_trinity] Trinity executable and Rust utils already symlinked"
        link_rust_utils
        # Run even on the cached-checkout fast path, in case an existing
        # install predates this fix and still has the broken archive.
        neutralize_broken_butterfly_cds
        exit 0
    fi
fi

echo "[pixi_install_trinity] Cloning Trinity from ${TRINITY_REPO} (commit ${TRINITY_COMMIT})..."

# Create the install parent dir if it doesn't exist
mkdir -p "${CONDA_PREFIX}/opt"

# Clone Trinity with submodules
if [ ! -d "${TRINITY_INSTALL_DIR}" ]; then
    git clone --recursive --jobs=4 -b "${TRINITY_COMMIT}" "${TRINITY_REPO}" "${TRINITY_INSTALL_DIR}"
    # Ensure all submodules are initialized and updated (belt-and-suspenders)
    git -C "${TRINITY_INSTALL_DIR}" submodule update --init --recursive 2>/dev/null || true
else
    echo "[pixi_install_trinity] Trinity directory already exists at ${TRINITY_INSTALL_DIR}"
fi

echo "[pixi_install_trinity] Building Trinity with make..."
# Trinity's bundled Inchworm/Chrysalis declare cmake_minimum_required(VERSION 3.1),
# which CMake >= 4 rejects ("Compatibility with CMake < 3.5 has been removed").
# CMake >= 3.31 honors this env var as the floor policy version, letting the old
# CMakeLists configure without patching upstream. Remove if Trinity bumps its
# cmake_minimum_required past 3.5.
export CMAKE_POLICY_VERSION_MINIMUM="${CMAKE_POLICY_VERSION_MINIMUM:-3.5}"
sed -i '1s/^/SHELL := \/bin\/bash\n/' "${TRINITY_INSTALL_DIR}/Makefile"
# `all:` (built below) depends on butterfly_cds_target, which runs
# build_butterfly_cds_archive.sh to (re)dump Butterfly/butterfly_cds.jsa using
# whatever `java` is on PATH at build time. Skip it: even a freshly-dumped,
# correctly-versioned archive still fails to load at runtime for this jar
# (Butterfly.jar's manifest carries `Class-Path: .`, which fails CDS's
# classpath validation regardless of archive version -- see
# neutralize_broken_butterfly_cds() above), so running this target would just
# spend build time re-creating a still-broken archive. Override it to a
# no-op; the rest of the build (inchworm/chrysalis, Rust utils) is unaffected,
# and the already-broken archive checked in by this same commit is deleted by
# neutralize_broken_butterfly_cds() after `make` finishes.
awk '
  /^butterfly_cds_target:/ { skip=1; next }
  skip && /^[^[:space:]]/ { skip=0 }
  skip { next }
  { print }
' "${TRINITY_INSTALL_DIR}/Makefile" > "${TRINITY_INSTALL_DIR}/Makefile.patched"
echo "butterfly_cds_target:" >> "${TRINITY_INSTALL_DIR}/Makefile.patched"
echo "	@:" >> "${TRINITY_INSTALL_DIR}/Makefile.patched"
mv "${TRINITY_INSTALL_DIR}/Makefile.patched" "${TRINITY_INSTALL_DIR}/Makefile"
make -C "${TRINITY_INSTALL_DIR}"

neutralize_broken_butterfly_cds

echo "[pixi_install_trinity] Running install.py to set up symlinks and install Rust utilities..."
# install.py handles all binaries and rust_bio tools setup. Its CLI takes an
# 'install' action plus --install-dir (default ~/.local/bin), not a bare path.
python3 "${TRINITY_INSTALL_DIR}/install.py" install --install-dir "${ENV_BIN_DIR}"

# install.py only links the main Trinity script; put the Rust utils on PATH too.
link_rust_utils

echo "[pixi_install_trinity] Trinity installation complete at ${TRINITY_INSTALL_DIR}"
echo "[pixi_install_trinity] Executables symlinked to ${ENV_BIN_DIR}"
