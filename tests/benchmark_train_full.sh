#!/usr/bin/env bash
# Benchmark the whole `funannotate train` step with Rust tools and without them.
#
# One run of `funannotate train` (hisat2, genome-guided Trinity, StringTie, PASA, TransDecoder, kallisto
# and the final model choice), start to end, for one engine:
#   rust   as shipped: Trinity's Rust utilities, and PASA's slclust_rust, cdbyank_rust and faidx_rust
#          with the C++ assembler.
#   perl   the original tools everywhere: TRINITY_NO_RUST=1, and all four Rust PASA tools hidden in a
#          private copy of the PASA tree (see tests/benchmark_pasa_rust.sh for how the PASA switch works).
#
# Each run proves its engine. After the run the script counts calls to rust_bio_utils binaries in the logs
# (the rust leg needs at least one, the perl leg none) and, before the run, checks which Rust PASA tools
# PATH shows. A run that fails a check is INVALID (exit 4).
#
# The stage profile (minutes per stage) is read from the console log by tests/benchmark_summarize.py.
#
# Usage:
#   tests/benchmark_train_full.sh <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus] [rep]

set -euo pipefail

DATASET_DIR="${1:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus] [rep]}"
BENCHMARK_OUT_DIR="${2:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus] [rep]}"
ENGINE="${3:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus] [rep]}"
CPUS="${4:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"
REP="${5:-1}"
SPECIES="${BENCH_SPECIES:-Cordyceps militaris benchmark}"
MEMORY="${BENCH_MEMORY:-16G}"
: "${CONDA_PREFIX:?CONDA_PREFIX is not set; activate the funannotate environment}"
PASA_INSTALL="${PASA_INSTALL:-$CONDA_PREFIX/opt/pasa}"

case "$ENGINE" in
    rust | perl) ;;
    *) echo "ERROR: engine must be rust or perl, got: $ENGINE" >&2; exit 1 ;;
esac

TEST_NAME="${ENGINE}_rep${REP}"
BENCHMARK_TIME="/usr/bin/time -v -a"
LOG_DIR="${BENCHMARK_OUT_DIR}/logs"
RESULTS_DIR="${BENCHMARK_OUT_DIR}/results"
SCRATCH_BASE="${SCRATCH:?SCRATCH is not set (node-local scratch)}"
RUST_TOOLS=(pasa_rust slclust_rust cdbyank_rust faidx_rust)
if [ "$ENGINE" = "perl" ]; then
    HIDDEN_TOOLS=(pasa_rust slclust_rust cdbyank_rust faidx_rust)
else
    HIDDEN_TOOLS=()
fi

log_info() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
log_error() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }

# funannotate makes some files in training/ as absolute symlinks into the work directory on node-local
# scratch. After copying back, those links point at a directory that no longer exists. Re-point them.
fix_scratch_links() {
    local l t new
    [ -d "$SYNC_TARGET_DIR/output" ] || return 0
    while IFS= read -r l; do
        t=$(readlink "$l")
        case "$t" in
            "$WORK_DIR"/*)
                new="$SYNC_TARGET_DIR/${t#"$WORK_DIR"/}"
                [ -e "$new" ] && ln -sfn "$new" "$l"
                ;;
        esac
    done < <(find "$SYNC_TARGET_DIR/output" -type l)
}

sync_work_dir_back() {
    if [ -d "$WORK_DIR/output" ]; then
        rsync -a --exclude='output/training/trinity_gg' --exclude='output/training/pasa/*.fasta*' \
            "$WORK_DIR/output" "$SYNC_TARGET_DIR/" || log_error "copy back to $SYNC_TARGET_DIR had errors"
        fix_scratch_links
    fi
}

setup_directories() {
    mkdir -p "$LOG_DIR" "$RESULTS_DIR"
    SYNC_TARGET_DIR="$BENCHMARK_OUT_DIR/work_train_${TEST_NAME}"
    WORK_DIR="$SCRATCH_BASE/funannotate_train_bench_${TEST_NAME}"
    rm -rf "$WORK_DIR"
    mkdir -p "$WORK_DIR" "$SYNC_TARGET_DIR"
    trap 'sync_work_dir_back' EXIT
    trap 'log_error "terminated; syncing partial results"; exit 143' TERM INT
}

validate_inputs() {
    [ -d "$DATASET_DIR" ] || { log_error "dataset directory not found: $DATASET_DIR"; exit 1; }
    FASTA_FILE=$(find "$DATASET_DIR" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) | head -1)
    LEFT_FASTQ=$(find "$DATASET_DIR" -maxdepth 1 -name "*R1*.fastq.gz" | head -1)
    RIGHT_FASTQ=$(find "$DATASET_DIR" -maxdepth 1 -name "*R2*.fastq.gz" | head -1)
    if [ -z "$FASTA_FILE" ] || [ -z "$LEFT_FASTQ" ] || [ -z "$RIGHT_FASTQ" ]; then
        log_error "need a genome FASTA and R1/R2 .fastq.gz files in $DATASET_DIR"
        exit 1
    fi
    [ -d "$PASA_INSTALL/src" ] && [ -d "$PASA_INSTALL/bin" ] || { log_error "PASA install not found at $PASA_INSTALL"; exit 1; }
    local home="${TRINITY_HOME:-${TRINITYHOME:-}}"
    if [ -z "$home" ] || [ ! -d "$home/rust_bio_utils/target/release" ]; then
        log_error "TRINITY_HOME/rust_bio_utils/target/release not found (TRINITY_HOME='$home')"
        exit 1
    fi
}

# Private copy of the PASA tree; the Rust tools are hidden there for the perl leg.
prepare_pasa_copy() {
    PASA_COPY="$WORK_DIR/pasa"
    cp -a "$PASA_INSTALL" "$PASA_COPY"
    [ -e "$PASA_COPY/src/bin" ] || { log_error "$PASA_COPY/src/bin is missing"; exit 1; }
    local t
    for t in ${HIDDEN_TOOLS[@]+"${HIDDEN_TOOLS[@]}"}; do
        if [ -e "$PASA_COPY/bin/$t" ]; then
            mv "$PASA_COPY/bin/$t" "$PASA_COPY/bin/${t}.disabled"
        fi
    done
    export PASAHOME="$PASA_COPY/src"
    unset PASA_ASSEMBLER
}

# The visible Rust PASA tools must equal the Rust tools minus the hidden ones.
verify_pasa_state() {
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
    echo "rust_pasa_tools_found=${found[*]:-none}" >> "$RESULTS_DIR/${TEST_NAME}.metrics"
    if [ "${found[*]:-}" != "${expected[*]:-}" ]; then
        log_error "$ENGINE leg sees PASA tools [${found[*]:-none}], expected [${expected[*]:-none}]; INVALID"
        echo "status=INVALID" >> "$RESULTS_DIR/${TEST_NAME}.metrics"
        exit 4
    fi
}

record_environment() {
    {
        echo "engine=$ENGINE"
        echo "rep=$REP"
        echo "stage=train_full"
        echo "cpus=$CPUS"
        echo "memory=$MEMORY"
        echo "host=$(hostname)"
        echo "cpu_model=$(lscpu | awk -F: '/Model name/ {gsub(/^ +/, "", $2); print $2; exit}')"
        echo "cpu_count_node=$(nproc --all)"
        echo "mem_total_kb=$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
        echo "slurm_job_id=${SLURM_JOB_ID:-none}"
        echo "slurm_partition=${SLURM_JOB_PARTITION:-none}"
        echo "TRINITY_NO_RUST=${TRINITY_NO_RUST:-unset}"
        echo "PASAHOME=$PASAHOME"
        echo "trinity_version=$(Trinity --version 2>&1 | head -1)"
        echo "pasa_version=$(grep -m1 '^my \$VERSION' "$PASAHOME/Launch_PASA_pipeline.pl" 2> /dev/null | tr -d '\n')"
        echo "funannotate_version=$(funannotate --version 2>&1 | head -1)"
        echo "start=$(date -Is)"
    } > "$RESULTS_DIR/${TEST_NAME}.env"
}

count_rust_calls() {
    local out="$1"
    { grep -rhs "rust_bio_utils/target/release" "$out/logfiles" "$out/training" --include='*.log' 2> /dev/null || true; } | wc -l
}

run_benchmark() {
    local out="$WORK_DIR/output"
    local test_log="$LOG_DIR/train_${TEST_NAME}.log"
    local test_metrics="$RESULTS_DIR/${TEST_NAME}.metrics"
    mkdir -p "$out"

    if [ "$ENGINE" = "perl" ]; then
        export TRINITY_NO_RUST=1
    else
        unset TRINITY_NO_RUST
    fi

    local rc=0
    $BENCHMARK_TIME -o "$test_metrics" \
        funannotate train \
            -i "$FASTA_FILE" \
            -o "$out" \
            --left "$LEFT_FASTQ" \
            --right "$RIGHT_FASTQ" \
            --left_norm "$LEFT_FASTQ" \
            --right_norm "$RIGHT_FASTQ" \
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

    local gff="$out/training/funannotate_train.pasa.gff3"
    if [ ! -s "$gff" ]; then
        log_error "expected train output not found: $gff"
        echo "status=FAILED" >> "$test_metrics"
        return 1
    fi
    echo "pasa_gff3_genes=$(grep -c -P '\tgene\t' "$gff" || true)" >> "$test_metrics"
    echo "pasa_gff3_struct_md5=$(grep -v '^#' "$gff" | awk -F'\t' 'NF>=7 {print $1"\t"$3"\t"$4"\t"$5"\t"$7}' | sort | md5sum | cut -d' ' -f1)" >> "$test_metrics"
    local fasta="$out/training/trinity.fasta"
    if [ -s "$fasta" ]; then
        echo "transcripts=$(grep -c '^>' "$fasta")" >> "$test_metrics"
    fi

    # Prove the Trinity engine.
    local calls
    calls=$(count_rust_calls "$out")
    echo "rust_calls_in_logs=$calls" >> "$test_metrics"
    if { [ "$ENGINE" = "rust" ] && [ "$calls" -eq 0 ]; } || { [ "$ENGINE" = "perl" ] && [ "$calls" -ne 0 ]; }; then
        log_error "Trinity engine check failed for $ENGINE (Rust calls in logs: $calls); INVALID"
        echo "status=INVALID" >> "$test_metrics"
        return 4
    fi
    echo "status=OK" >> "$test_metrics"
    log_info "OK: $TEST_NAME"
}

main() {
    log_info "full train $ENGINE benchmark, replicate $REP"
    setup_directories
    validate_inputs
    prepare_pasa_copy
    : > "$RESULTS_DIR/${TEST_NAME}.metrics"
    record_environment
    verify_pasa_state
    run_benchmark
}

main "$@"
