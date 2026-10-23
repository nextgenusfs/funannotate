#!/usr/bin/env bash
# Check that the macOS fixes took effect in the active environment (macOS arm64).
# Run after `pixi run install-externals`:   pixi run bash install_scripts/macos_smoke_checks.sh
# Exit status 0 only when every check passes. Each check is one line: "ok:" or "FAIL:".

set -uo pipefail

: "${CONDA_PREFIX:?CONDA_PREFIX is not set; run inside the environment (pixi run ...)}"

fail=0
ok() { echo "ok:   $1"; }
bad() { echo "FAIL: $1"; fail=1; }

# 1. Augustus: the --proteinprofile test that funannotate runs before BUSCO training. The clang
#    build in the lock crashes (exit 139); the GCC build from install_augustus_gcc_macos.sh passes.
if python - << 'PY' > /dev/null 2>&1; then ok "augustus --proteinprofile (funannotate checkAugustusFunc)"; else bad "augustus --proteinprofile (funannotate checkAugustusFunc)"; fi
import logging
import funannotate.library as lib
lib.log = logging.getLogger("smoke")
version, passed = lib.checkAugustusFunc()
raise SystemExit(0 if passed else 1)
PY

# 2. perl DB_File: the arm64 conda package crashes; perl's core copy must come first in PERL5LIB.
if perl -e 'use DB_File' > /dev/null 2>&1; then ok "perl loads DB_File"; else bad "perl loads DB_File"; fi

# 3. blat: an invalid code signature gets it killed (exit 137) and reported as not installed.
if /usr/bin/codesign -v "${CONDA_PREFIX}/bin/blat" > /dev/null 2>&1; then ok "blat has a valid signature"; else bad "blat has a valid signature"; fi

# 4. PASA's C++ assembler: needs an rpath for libomp, or it stops with "Library not loaded".
pasa="${CONDA_PREFIX}/opt/pasa/bin/pasa"
rpaths=$(otool -l "${pasa}" 2> /dev/null | grep -c LC_RPATH || true)
out=$("${pasa}" < /dev/null 2>&1 || true)
if [ "${rpaths}" -gt 0 ] && ! printf '%s' "${out}" | grep -q "Library not loaded"; then ok "pasa has an rpath and starts"; else bad "pasa has an rpath and starts"; fi

# 5. pblat: must not crash with several threads (the shim runs it with one).
# The input must have many query sequences: pblat only runs threads in parallel across
# queries, and a single short query does not crash the unshimmed binary at all (so such a
# check passes even without the shim). With 300 queries the real pblat crashed in 10 of 10
# runs on Apple Silicon, and a single thread finishes in under a second.
tmp=$(mktemp -d)
python - "${tmp}" << 'PY' > /dev/null 2>&1 || true
import random
import sys

random.seed(7)
d = sys.argv[1]
genome = ["".join(random.choice("ACGT") for _ in range(60000)) for _ in range(2)]
with open(d + "/g.fa", "w") as out:
    for i, s in enumerate(genome):
        out.write(">g%d\n%s\n" % (i, s))
with open(d + "/q.fa", "w") as out:
    for i in range(300):
        s = random.choice(genome)
        a = random.randrange(0, len(s) - 800)
        out.write(">q%d\n%s\n" % (i, s[a : a + 600]))
PY
rc=0
if [ -s "${tmp}/q.fa" ]; then
    pblat "${tmp}/g.fa" "${tmp}/q.fa" -q=rna -threads=4 -out=pslx "${tmp}/o.pslx" > /dev/null 2>&1 || rc=$?
else
    rc=255
fi
rm -rf "${tmp}"
if [ "${rc}" -lt 128 ]; then ok "pblat -threads=4 does not crash (exit ${rc})"; else bad "pblat -threads=4 crashed (exit ${rc})"; fi

# 6. funannotate check must find the tools that these fixes provide.
report=$(funannotate check 2>&1 || true)
missing=$(printf '%s\n' "${report}" | grep -E "ERROR: (fasta|tbl2asn|blat|DB_File|augustus|PASA|Trinity)" || true)
if [ -z "${missing}" ]; then ok "funannotate check finds fasta, tbl2asn, blat, DB_File, augustus"; else bad "funannotate check reports: ${missing}"; fi

if [ "${fail}" -eq 0 ]; then
    echo "all macOS smoke checks passed"
else
    echo "some macOS smoke checks FAILED" >&2
fi
exit "${fail}"
