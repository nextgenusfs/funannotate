#!/usr/bin/env bash
# Time the Rust EVM GFF3 converter: old build vs new build, on shared storage vs
# node-local $SCRATCH, at several thread counts.
#
# Input is the EVM directory of a finished `funannotate predict` run (the one
# that holds partitions_list.out and one evm.out per contig). The script copies
# the first N unpartitioned evm.out files into two test trees (one under
# --shared-root, one under $SCRATCH), so the real run is not touched.
#
# Usage (inside a SLURM job so $SCRATCH exists and CPUs are reserved):
#   sbatch -p epyc -c 8 --mem 16G -t 02:00:00 --wrap \
#     "tests/benchmark_evm_convert.sh <EVM_dir> <old_bin_dir> <new_bin_dir> <shared_root> [N] [REPS] [OUT_TSV]"
#
#   old_bin_dir / new_bin_dir  directories holding convert_EVM_outputs_to_GFF3
#       built from the commit before and after the optimization changes
#       (cargo build --release with a separate CARGO_TARGET_DIR for each).
#   shared_root  a directory on the shared filesystem (/bigdata) for the test tree.
#   N     contigs to copy (default 20000)
#   REPS  repeats per case (default 3)
#
# Output: a TSV of wall, user and system seconds and files/s per case, plus a
# check that old and new produce the same GFF3 bytes in listing order. Page-cache
# warmth is not controlled: each case is run REPS times in rotating order, and
# the first pass over each tree is a discarded warm-up.
set -euo pipefail

EVM_DIR=${1:?EVM dir}
OLD_BIN=${2:?old bin dir}
NEW_BIN=${3:?new bin dir}
SHARED_ROOT=${4:?shared root dir}
N=${5:-20000}
REPS=${6:-3}
OUT_TSV=${7:-evm_convert_benchmark.tsv}
SCRATCH_DIR=${SCRATCH:?SCRATCH is not set; run inside a SLURM job}
BIN=convert_EVM_outputs_to_GFF3
THREADS_LIST="1 4 8"

for b in "$OLD_BIN/$BIN" "$NEW_BIN/$BIN"; do
    [ -x "$b" ] || { echo "missing $b" >&2; exit 1; }
done
LIST_SRC="$EVM_DIR/partitions_list.out"
[ -s "$LIST_SRC" ] || { echo "missing $LIST_SRC" >&2; exit 1; }

# ---- build the two test trees -------------------------------------------------
build_tree() {
    local root=$1 i=0 acc dir flag
    rm -rf "$root"; mkdir -p "$root"
    : > "$root/partitions_list.out"; : > "$root/gff3_paths.txt"
    while IFS=$'\t' read -r acc dir flag _; do
        [ "$flag" = "N" ] || continue
        [ -s "$dir/evm.out" ] || continue
        i=$((i + 1))
        mkdir -p "$root/c$i"
        cp "$dir/evm.out" "$root/c$i/evm.out"
        printf '%s\t%s\tN\n' "$acc" "$root/c$i" >> "$root/partitions_list.out"
        printf '%s\n' "$root/c$i/evm.out.gff3" >> "$root/gff3_paths.txt"
        [ "$i" -ge "$N" ] && break
    done < "$LIST_SRC"
    echo "tree $root: $i contigs" >&2
}
SHARED_TREE="$SHARED_ROOT/evm_convert_bench"
SCRATCH_TREE="$SCRATCH_DIR/evm_convert_bench"
build_tree "$SHARED_TREE"
build_tree "$SCRATCH_TREE"

clean_outputs() { # $1 tree
    xargs -a "$1/gff3_paths.txt" rm -f; rm -f "$1/all.gff3"
}
listing_md5() { xargs -a "$1/gff3_paths.txt" cat 2>/dev/null | md5sum | cut -d' ' -f1; }

printf 'tree\tcase\trep\twall_s\tuser_s\tsys_s\tfiles_per_s\tmd5\n' > "$OUT_TSV"

run_case() { # tree_name tree case_name bin extra_args...
    local tname=$1 tree=$2 cname=$3 bin=$4 rep=$5; shift 5
    local tfile; tfile=$(mktemp)
    clean_outputs "$tree"
    ( cd "$tree" && /usr/bin/time -f '%e\t%U\t%S' -o "$tfile" \
        "$bin" --partitions partitions_list.out -O evm.out "$@" > /dev/null 2>&1 )
    local wall user sys md5 n
    IFS=$'\t' read -r wall user sys < "$tfile"; rm -f "$tfile"
    n=$(wc -l < "$tree/partitions_list.out")
    if [ -s "$tree/all.gff3" ]; then md5=$(md5sum < "$tree/all.gff3" | cut -d' ' -f1)
    else md5=$(listing_md5 "$tree"); fi
    printf '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n' "$tname" "$cname" "$rep" "$wall" "$user" "$sys" \
        "$(awk -v n="$n" -v w="$wall" 'BEGIN{printf "%.1f", n/(w>0?w:0.01)}')" "$md5" >> "$OUT_TSV"
}

# warm-up pass (discarded): reads every evm.out once
for t in "$SHARED_TREE" "$SCRATCH_TREE"; do
    run_case warmup "$t" warmup "$OLD_BIN/$BIN" 0 >/dev/null 2>&1 || true
done
sed -i '/^warmup\t/d' "$OUT_TSV"

for rep in $(seq 1 "$REPS"); do
    for pair in "shared:$SHARED_TREE" "scratch:$SCRATCH_TREE"; do
        tname=${pair%%:*}; tree=${pair#*:}
        run_case "$tname" "$tree" old "$OLD_BIN/$BIN" "$rep"
        for th in $THREADS_LIST; do
            run_case "$tname" "$tree" "new_t$th" "$NEW_BIN/$BIN" "$rep" --threads "$th"
        done
        run_case "$tname" "$tree" "new_t8_skip_concat" "$NEW_BIN/$BIN" "$rep" \
            --threads 8 --skip_empty --concat_out all.gff3
    done
done

echo "wrote $OUT_TSV" >&2
# Byte comparison: old vs new_t* must match. new_t8_skip_concat drops empty
# contigs and concatenates directly, so it matches only if no contig is empty
# or empty contigs produce no output in the old build either.
awk -F'\t' 'NR>1 && $2!="new_t8_skip_concat"{m[$1]=m[$1] " " $8} END{for(t in m){split(m[t],a," "); u=0; delete s; for(i in a) s[a[i]]=1; for(k in s)u++; print t": distinct md5 among old and new_t*: " u}}' "$OUT_TSV" >&2
