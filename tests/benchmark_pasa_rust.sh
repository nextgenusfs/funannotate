#!/usr/bin/env bash
# Benchmark the PASA stage of `funannotate train` with one engine.
#
# Engines:
#   rust     PASA as shipped: C++ assembler `pasa`, plus slclust_rust,
#            cdbyank_rust and faidx_rust (PASA prefers them when present).
#   perl     the same PASA with all four Rust tools hidden: C++ `slclust`,
#            C++ `cdbyank`, no faidx_rust. ("perl" names the non-Rust baseline,
#            as in the other benchmark scripts.)
#   rustasm  `rust` plus the Rust assembler: PASA_ASSEMBLER=pasa_rust.
#   nocdb    `rust` with cdbyank_rust and faidx_rust hidden (C++ cdbyank is used).
#   noslclust `rust` with slclust_rust hidden (C++ slclust is used).
# nocdb and noslclust attribute any difference between rust and perl to one
# tool. Compare rust with nocdb to measure cdbyank_rust alone.
#
# How the engine is switched. PASA finds its tools by scanning PATH, and
# Launch_PASA_pipeline.pl puts $PASAHOME/bin first. The script copies the PASA
# install to node-local scratch and renames the Rust binaries in the copy for
# the perl leg, so the shared install is never changed. PASA's logs do not
# record which binaries ran, so the proof is the file state: before the run the
# script checks that the perl leg finds no Rust tool and the rust leg finds all
# of them, using the same PATH that PASA will use. A run that fails this check
# is INVALID (exit 4).
#
# Stage. The run reuses an existing Trinity transcript FASTA (--trinity), so
# Trinity does not run. The PASA stage time is taken from two lines in
# funannotate-train.log ("Running PASA alignment step" and "PASA finished"). The
# whole-run time is recorded too. hisat2 and other non-Rust work in `train` is
# the same in all legs.
#
# Usage:
#   TRINITY_FASTA=/path/to/funannotate_train.trinity-GG.fasta \
#   tests/benchmark_pasa_rust.sh <dataset_dir> <benchmark_out_dir> <rust|perl|rustasm|nocdb|noslclust> [cpus] [rep]

set -euo pipefail

DATASET_DIR="${1:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl|rustasm|nocdb|noslclust> [cpus] [rep]}"
BENCHMARK_OUT_DIR="${2:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl|rustasm|nocdb|noslclust> [cpus] [rep]}"
ENGINE="${3:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl|rustasm|nocdb|noslclust> [cpus] [rep]}"
CPUS="${4:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"
REP="${5:-1}"
TRINITY_FASTA="${TRINITY_FASTA:?set TRINITY_FASTA to a finished Trinity-GG transcript FASTA}"
SPECIES="${BENCH_SPECIES:-Cordyceps militaris benchmark}"
MEMORY="${BENCH_MEMORY:-16G}"
: "${CONDA_PREFIX:?CONDA_PREFIX is not set; activate the funannotate environment}"
PASA_INSTALL="${PASA_INSTALL:-$CONDA_PREFIX/opt/pasa}"

case "$ENGINE" in
    rust | perl | rustasm | nocdb | noslclust) ;;
    *) echo "ERROR: engine must be rust, perl, rustasm, nocdb or noslclust, got: $ENGINE" >&2; exit 1 ;;
esac

TEST_NAME="${ENGINE}_rep${REP}"
BENCHMARK_TIME="/usr/bin/time -v -a"
LOG_DIR="${BENCHMARK_OUT_DIR}/logs"
RESULTS_DIR="${BENCHMARK_OUT_DIR}/results"
SCRATCH_BASE="${SCRATCH:?SCRATCH is not set (node-local scratch)}"
RUST_TOOLS=(pasa_rust slclust_rust cdbyank_rust faidx_rust)
# Tools that this engine hides in the private PASA copy.
case "$ENGINE" in
    perl) HIDDEN_TOOLS=(pasa_rust slclust_rust cdbyank_rust faidx_rust) ;;
    nocdb) HIDDEN_TOOLS=(cdbyank_rust faidx_rust) ;;
    noslclust) HIDDEN_TOOLS=(slclust_rust) ;;
    *) HIDDEN_TOOLS=() ;;
esac

log_info() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
log_error() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }

setup_directories() {
    mkdir -p "$LOG_DIR" "$RESULTS_DIR"
    SYNC_TARGET_DIR="$BENCHMARK_OUT_DIR/work_pasa_${TEST_NAME}"
    WORK_DIR="$SCRATCH_BASE/funannotate_pasa_bench_${TEST_NAME}"
    rm -rf "$WORK_DIR"
    mkdir -p "$WORK_DIR" "$SYNC_TARGET_DIR"
    trap 'sync_work_dir_back' EXIT
    trap 'log_error "terminated; syncing partial results"; exit 143' TERM INT
}

sync_work_dir_back() {
    if [ -d "$WORK_DIR/output" ]; then
        rsync -a --exclude='output/training/pasa/*.fasta*' "$WORK_DIR/output" "$SYNC_TARGET_DIR/" \
            || log_error "copy back to $SYNC_TARGET_DIR had errors"
    fi
}

validate_inputs() {
    [ -d "$DATASET_DIR" ] || { log_error "dataset directory not found: $DATASET_DIR"; exit 1; }
    [ -s "$TRINITY_FASTA" ] || { log_error "TRINITY_FASTA is empty or missing: $TRINITY_FASTA"; exit 1; }
    [ -d "$PASA_INSTALL/src" ] && [ -d "$PASA_INSTALL/bin" ] || {
        log_error "PASA install not found at $PASA_INSTALL (need src/ and bin/)"; exit 1; }
    FASTA_FILE=$(find "$DATASET_DIR" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) | head -1)
    LEFT_FASTQ=$(find "$DATASET_DIR" -maxdepth 1 -name "*R1*.fastq.gz" | head -1)
    RIGHT_FASTQ=$(find "$DATASET_DIR" -maxdepth 1 -name "*R2*.fastq.gz" | head -1)
    if [ -z "$FASTA_FILE" ] || [ -z "$LEFT_FASTQ" ] || [ -z "$RIGHT_FASTQ" ]; then
        log_error "need a genome FASTA and R1/R2 .fastq.gz files in $DATASET_DIR"
        exit 1
    fi
}

# Copy the PASA tree to scratch and hide the Rust tools there for the perl leg.
prepare_pasa_copy() {
    PASA_COPY="$WORK_DIR/pasa"
    log_info "copying $PASA_INSTALL to $PASA_COPY"
    cp -a "$PASA_INSTALL" "$PASA_COPY"
    if [ ! -e "$PASA_COPY/src/bin" ]; then
        log_error "$PASA_COPY/src/bin is missing (expected symlink to ../bin)"
        exit 1
    fi
    local t
    for t in ${HIDDEN_TOOLS[@]+"${HIDDEN_TOOLS[@]}"}; do
        if [ -e "$PASA_COPY/bin/$t" ]; then
            mv "$PASA_COPY/bin/$t" "$PASA_COPY/bin/${t}.disabled"
        fi
    done
    export PASAHOME="$PASA_COPY/src"
    if [ "$ENGINE" = "rustasm" ]; then
        export PASA_ASSEMBLER=pasa_rust
    else
        unset PASA_ASSEMBLER
    fi
}

# Check which Rust tools are visible on the PATH that PASA builds (its own
# bin/ first, then the rest of PATH). The visible set must equal the Rust tools
# minus the ones this engine hides.
verify_engine_state() {
    local found=() expected=() t h hidden
    for t in "${RUST_TOOLS[@]}"; do
        if PATH="$PASAHOME/bin:$PATH" command -v "$t" > /dev/null 2>&1; then
            found+=("$t")
        fi
        hidden=0
        for h in ${HIDDEN_TOOLS[@]+"${HIDDEN_TOOLS[@]}"}; do
            [ "$h" = "$t" ] && hidden=1
        done
        [ "$hidden" -eq 0 ] && expected+=("$t")
    done
    echo "rust_tools_found=${found[*]:-none}" >> "$RESULTS_DIR/${TEST_NAME}.metrics"
    echo "rust_tools_expected=${expected[*]:-none}" >> "$RESULTS_DIR/${TEST_NAME}.metrics"
    if [ "${found[*]:-}" != "${expected[*]:-}" ]; then
        log_error "$ENGINE leg sees [${found[*]:-none}] but expected [${expected[*]:-none}]; the run is INVALID"
        echo "status=INVALID" >> "$RESULTS_DIR/${TEST_NAME}.metrics"
        exit 4
    fi
    log_info "engine check ok: visible Rust tools = [${found[*]:-none}]"
}

record_environment() {
    {
        echo "engine=$ENGINE"
        echo "rep=$REP"
        echo "stage=pasa"
        echo "cpus=$CPUS"
        echo "memory=$MEMORY"
        echo "host=$(hostname)"
        echo "cpu_model=$(lscpu | awk -F: '/Model name/ {gsub(/^ +/, "", $2); print $2; exit}')"
        echo "cpu_count_node=$(nproc --all)"
        echo "mem_total_kb=$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
        echo "slurm_job_id=${SLURM_JOB_ID:-none}"
        echo "slurm_partition=${SLURM_JOB_PARTITION:-none}"
        echo "PASAHOME=$PASAHOME"
        echo "PASA_ASSEMBLER=${PASA_ASSEMBLER:-unset}"
        echo "pasa_version=$(grep -m1 '^my \$VERSION' "$PASAHOME/Launch_PASA_pipeline.pl" 2> /dev/null | tr -d '\n')"
        echo "trinity_fasta=$TRINITY_FASTA"
        echo "funannotate_version=$(funannotate --version 2>&1 | head -1)"
        echo "start=$(date -Is)"
    } > "$RESULTS_DIR/${TEST_NAME}.env"
}

# Seconds between two "[MM/DD/YY HH:MM:SS]: ..." funannotate log lines.
stage_seconds() {
    local logfile="$1" start_pat="$2" end_pat="$3" s e
    s=$(grep -m1 -- "$start_pat" "$logfile" | sed -E 's/^\[([0-9]{2})\/([0-9]{2})\/([0-9]{2}) ([0-9:]{8})\].*/20\3-\1-\2 \4/')
    e=$(grep -m1 -- "$end_pat" "$logfile" | sed -E 's/^\[([0-9]{2})\/([0-9]{2})\/([0-9]{2}) ([0-9:]{8})\].*/20\3-\1-\2 \4/')
    if [ -z "$s" ] || [ -z "$e" ]; then
        echo "NA"
        return
    fi
    echo $(($(date -d "$e" +%s) - $(date -d "$s" +%s)))
}

run_benchmark() {
    local out="$WORK_DIR/output"
    local test_log="$LOG_DIR/pasa_${TEST_NAME}.log"
    local test_metrics="$RESULTS_DIR/${TEST_NAME}.metrics"
    mkdir -p "$out"

    local rc=0
    $BENCHMARK_TIME -o "$test_metrics" \
        funannotate train \
            -i "$FASTA_FILE" \
            -o "$out" \
            --left "$LEFT_FASTQ" \
            --right "$RIGHT_FASTQ" \
            --left_norm "$LEFT_FASTQ" \
            --right_norm "$RIGHT_FASTQ" \
            --trinity "$TRINITY_FASTA" \
            --species "$SPECIES" \
            --cpus "$CPUS" \
            --memory "$MEMORY" \
            --jaccard_clip \
            &> "$test_log" || rc=$?

    echo "end=$(date -Is)" >> "$RESULTS_DIR/${TEST_NAME}.env"
    echo "exit_status=$rc" >> "$test_metrics"
    if [ "$rc" -ne 0 ]; then
        log_error "funannotate train exited with $rc; see $test_log"
        echo "status=FAILED" >> "$test_metrics"
        return 1
    fi

    local trainlog="$out/logfiles/funannotate-train.log"
    echo "pasa_stage_seconds=$(stage_seconds "$trainlog" 'Running PASA alignment step' 'PASA finished')" >> "$test_metrics"

    # Output summary, to check that the engines give the same result.
    local gff="$out/training/funannotate_train.pasa.gff3"
    if [ -s "$gff" ]; then
        echo "pasa_gff3_genes=$(grep -c -P '\tgene\t' "$gff" || true)" >> "$test_metrics"
        echo "pasa_gff3_md5=$(sort "$gff" | md5sum | cut -d' ' -f1)" >> "$test_metrics"
    else
        log_error "expected PASA output not found: $gff"
        echo "status=FAILED" >> "$test_metrics"
        return 1
    fi
    echo "status=OK" >> "$test_metrics"
    log_info "OK: $TEST_NAME"
}

main() {
    log_info "PASA $ENGINE benchmark, replicate $REP"
    setup_directories
    validate_inputs
    prepare_pasa_copy
    : > "$RESULTS_DIR/${TEST_NAME}.metrics"
    record_environment
    verify_engine_state
    run_benchmark
}

main "$@"
