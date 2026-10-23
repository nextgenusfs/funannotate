#!/usr/bin/env bash
# Remove what a finished funannotate environment does not run: the Rust toolchain, the C/C++/Fortran
# compilers, cmake, headers, the JDK-only parts of openjdk, documentation and Python caches.
#
# For the Docker build stage (Dockerfile.base), after the source builds (Trinity, PASA, EVM, bowtie2) and
# before the env is copied into the runtime stage. Linux only (it uses ldd). Sizes are from the
# image-size-report workflow on v1.9.1-beta.3 (the env was 7.8 GB: rustlib 390 MB, compilers 670 MB,
# include 334 MB, LLVM/clang libraries ~700 MB, jvm 310 MB).
#
#   bash install_scripts/prune_runtime_env.sh [PREFIX]      # PREFIX defaults to $CONDA_PREFIX
#   PRUNE_TOOLCHAIN=0 ...    keep the Rust toolchain, compilers, cmake and headers (docs/JDK/caches are still pruned)
#   DRY_RUN=1 ...            print what would go, remove nothing
#
# Kept on purpose: mysqld/mariadb and the mysql client tools (PASA's --pasa_db mysql, local runs), R and the
# Bioconductor stack (TransDecoder), both BLAST+ suites, RepeatMasker data, java (trimmomatic, Trinity's Butterfly).
# An LLVM/clang library is removed only if no remaining binary or library depends on it (checked with ldd), and
# the script fails if the prune leaves any binary with a library "not found" that it did not have before.

set -euo pipefail

PREFIX="${1:-${CONDA_PREFIX:-}}"
[ -n "${PREFIX}" ] && [ -d "${PREFIX}/conda-meta" ] || { echo "usage: $0 PREFIX (a conda/pixi environment)" >&2; exit 1; }
[ "$(uname -s)" = "Linux" ] || { echo "prune_runtime_env.sh is Linux only (uses ldd)" >&2; exit 1; }
PRUNE_TOOLCHAIN="${PRUNE_TOOLCHAIN:-1}"
DRY_RUN="${DRY_RUN:-0}"

size() { du -sk "${PREFIX}" | cut -f1; }
rm_paths() { # remove files/dirs matching the globs given (relative to PREFIX), if any
    local p f
    for p in "$@"; do
        for f in "${PREFIX}"/${p}; do
            [ -e "${f}" ] || [ -L "${f}" ] || continue
            if [ "${DRY_RUN}" = "1" ]; then echo "would remove ${f#"${PREFIX}"/} ($(du -sk "${f}" 2> /dev/null | cut -f1) KB)"; else rm -rf "${f}"; fi
        done
    done
}

# Files with a missing shared library: "<file>: <lib> => not found", sorted. Run before and after.
missing_libs() {
    LD_LIBRARY_PATH="${PREFIX}/lib" find "${PREFIX}/bin" "${PREFIX}/lib" "${PREFIX}/libexec" "${PREFIX}/opt" -type f \
        \( -perm -u+x -o -name '*.so*' \) -print0 2> /dev/null \
        | xargs -0 -r -n 50 sh -c 'for f in "$@"; do ldd "$f" 2> /dev/null | grep "not found" | sed "s|^|$f: |"; done' sh \
        | sed 's/[[:space:]]\+/ /g' | sort -u
}

before_kb="$(size)"
echo "[prune_runtime_env] ${PREFIX}: $((before_kb / 1024)) MB"
missing_before="$(missing_libs)"

# Always: documentation, caches, conda leftovers, and what only the JDK (not the JRE) needs. The JDK stays
# in the environment because Trinity's make compiles Butterfly with it; this runs after the build.
rm_paths docs man share/doc share/man share/info share/gir-1.0 __clobbers__
if [ "${DRY_RUN}" != "1" ]; then
    find "${PREFIX}" -name __pycache__ -type d -prune -exec rm -rf {} +
    find "${PREFIX}" -name '*.pyc' -delete
fi
rm_paths lib/jvm/jmods lib/jvm/include lib/jvm/lib/src.zip lib/jvm/bin/javac lib/jvm/bin/javadoc lib/jvm/bin/jlink \
    lib/jvm/bin/jmod lib/jvm/bin/jpackage lib/jvm/bin/jshell lib/jvm/bin/jdb lib/jvm/bin/jconsole lib/jvm/bin/jdeps \
    bin/javac bin/javadoc bin/jlink bin/jmod bin/jpackage bin/jshell bin/jdb bin/jconsole bin/jdeps
# pandoc (158 MB) is only a dependency of r-rmarkdown, which funannotate never renders
rm_paths bin/pandoc

if [ "${PRUNE_TOOLCHAIN}" = "1" ]; then
    # Rust toolchain: the Rust EVM/PASA/Trinity utilities are built binaries; nothing calls cargo or rustc at runtime
    rm_paths lib/rustlib 'lib/librustc_driver-*.so' 'lib/libLLVM*-rust-*' bin/cargo 'bin/cargo-*' bin/rustc bin/rustdoc \
        bin/rustfmt bin/clippy-driver 'bin/rust-*'
    # C/C++/Fortran compilers, the glibc sysroot they use, and cmake
    rm_paths x86_64-conda-linux-gnu lib/gcc libexec/gcc 'bin/x86_64-conda-linux-gnu-*' bin/gcc bin/g++ bin/cc bin/c++ \
        bin/cpp bin/gfortran bin/cmake bin/ctest bin/cpack bin/ccmake 'share/cmake-*' include
    # LLVM/clang libraries: keep any that a remaining binary or library still needs
    cand="$(find "${PREFIX}/lib" -maxdepth 1 \( -name 'libLLVM*.so*' -o -name 'libclang*.so*' \) | sort)"
    if [ -n "${cand}" ]; then
        needed="$(mktemp)"
        # sonames needed by everything that is not itself a candidate
        export PREFIX
        LD_LIBRARY_PATH="${PREFIX}/lib" find "${PREFIX}/bin" "${PREFIX}/lib" "${PREFIX}/libexec" "${PREFIX}/opt" -type f \
            \( -perm -u+x -o -name '*.so*' \) -print0 2> /dev/null \
            | xargs -0 -r -n 50 sh -c 'for f in "$@"; do
                case "$f" in "${PREFIX}"/lib/libLLVM*|"${PREFIX}"/lib/libclang*) continue ;; esac
                ldd "$f" 2> /dev/null | sed -n "s/^[[:space:]]*\([^ ]*\) => .*/\1/p"
              done' sh \
            | sort -u > "${needed}"
        for lib in ${cand}; do
            if grep -qxF "$(basename "${lib}")" "${needed}"; then
                echo "[prune_runtime_env] keeping ${lib#"${PREFIX}"/}: still needed"
            else
                rm_paths "${lib#"${PREFIX}"/}"
            fi
        done
        rm -f "${needed}"
    fi
fi

if [ "${DRY_RUN}" = "1" ]; then exit 0; fi

missing_after="$(missing_libs)"
new_missing="$(comm -13 <(printf '%s\n' "${missing_before}") <(printf '%s\n' "${missing_after}") | sed '/^$/d')"
if [ -n "${new_missing}" ]; then
    echo "[prune_runtime_env] ERROR: the prune left libraries missing:" >&2
    echo "${new_missing}" | head -40 >&2
    exit 1
fi

# Things that must still run
java -version 2>&1 | head -1
mysqld --version
blastn -version | head -1
python -c 'import numpy, pandas, scipy, Bio'
Rscript -e 'invisible(library(stats))'

after_kb="$(size)"
echo "[prune_runtime_env] done: $((before_kb / 1024)) MB -> $((after_kb / 1024)) MB (saved $(((before_kb - after_kb) / 1024)) MB)"
