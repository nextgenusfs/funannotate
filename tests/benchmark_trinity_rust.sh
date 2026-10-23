#!/usr/bin/env bash
# Benchmark funannotate's Trinity stage with one engine (rust or perl).
#
# Runs `funannotate train` once and collects wall time and memory for that one
# run. Run it once per engine and replicate, as separate SLURM jobs (see
# benchmark_trinity_rust.sbatch), then merge with
# install_scripts/compare_trinity_rust_results.sh.
#
# How the engine is switched. Trinity's Perl scripts find the Rust utilities in
# the install tree ($TRINITY_HOME/rust_bio_utils/target/release), not on PATH.
# The only switch is the environment variable TRINITY_NO_RUST=1, which the
# fork's find_rust_binary() honors. Hiding binaries on PATH does nothing. An
# earlier version of this script did that, so both of its legs ran the Rust
# utilities and its "comparison" was a replicate pair.
#
# The script proves which engine ran. After the run it counts calls to
# rust_bio_utils binaries in the train and Trinity logs. The rust leg must show
# at least one. The perl leg must show none. A run that fails this check is
# marked INVALID and the script exits 4.
#
# Stage. By default (BENCH_STAGE=trinity) it passes --stop_after_trinity, so the
# timing covers read preparation, hisat2, normalization and genome-guided
# Trinity, and PASA does not run. BENCH_STAGE=full runs all of `funannotate
# train`.
#
# Usage:
#   tests/benchmark_trinity_rust.sh <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus] [rep]
#
# Example:
#   tests/benchmark_trinity_rust.sh \
#       ~/projects/funannotate/trinity_example/Cordyceps_militaris \
#       ~/projects/funannotate/BENCHMARK/trinity_clean rust 32 1

set -euo pipefail

DATASET_DIR="${1:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus] [rep]}"
BENCHMARK_OUT_DIR="${2:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus] [rep]}"
ENGINE="${3:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus] [rep]}"
CPUS="${4:-${SLURM_CPUS_PER_TASK:-$(nproc)}}"
REP="${5:-1}"
BENCH_STAGE="${BENCH_STAGE:-trinity}"
SPECIES="${BENCH_SPECIES:-Cordyceps militaris benchmark}"
MEMORY="${BENCH_MEMORY:-16G}"

case "$ENGINE" in
    rust) ;;
    perl) ;;
    *) echo "ERROR: engine must be 'rust' or 'perl', got: $ENGINE" >&2; exit 1 ;;
esac
case "$BENCH_STAGE" in
    trinity | full) ;;
    *) echo "ERROR: BENCH_STAGE must be 'trinity' or 'full', got: $BENCH_STAGE" >&2; exit 1 ;;
esac

TEST_NAME="${ENGINE}_rep${REP}"
BENCHMARK_TIME="/usr/bin/time -v"
# Small, important files go straight to BENCHMARK_OUT_DIR (shared storage).
# Trinity writes thousands of small files, so the work directory goes on
# node-local $SCRATCH and only the final results are copied back.
LOG_DIR="${BENCHMARK_OUT_DIR}/logs"
RESULTS_DIR="${BENCHMARK_OUT_DIR}/results"
SCRATCH_BASE="${SCRATCH:-$BENCHMARK_OUT_DIR}"

log_info() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] $*"; }
log_error() { echo "[$(date +'%Y-%m-%d %H:%M:%S')] ERROR: $*" >&2; }

setup_directories() {
    mkdir -p "$LOG_DIR" "$RESULTS_DIR"
    SYNC_TARGET_DIR="$BENCHMARK_OUT_DIR/work_${TEST_NAME}"
    WORK_DIR="$SCRATCH_BASE/funannotate_trinity_bench_${TEST_NAME}"
    rm -rf "$WORK_DIR"
    mkdir -p "$WORK_DIR" "$SYNC_TARGET_DIR"
    # Copy results back on normal exit, on error, and when SLURM ends the job.
    trap 'sync_work_dir_back' EXIT
    trap 'log_error "terminated (for example SLURM walltime); syncing partial results"; exit 143' TERM INT
}

# funannotate makes some files in training/ as absolute symlinks into the work
# directory on node-local scratch. After copying back, those links point at a
# directory that no longer exists. Re-point them at the copied files.
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

# Copy logs and final files back. Skip trinity_gg/, the large per-partition
# directory that nothing downstream reads.
sync_work_dir_back() {
    if [ -d "$WORK_DIR/output" ]; then
        rsync -a --exclude='output/training/trinity_gg' "$WORK_DIR/output" "$SYNC_TARGET_DIR/" \
            || log_error "copy back to $SYNC_TARGET_DIR had errors"
        fix_scratch_links
    fi
}

validate_dataset() {
    [ -d "$DATASET_DIR" ] || { log_error "dataset directory not found: $DATASET_DIR"; exit 1; }
    FASTA_FILE=$(find "$DATASET_DIR" -maxdepth 1 \( -name "*.fasta" -o -name "*.fa" \) | head -1)
    LEFT_FASTQ=$(find "$DATASET_DIR" -maxdepth 1 -name "*R1*.fastq.gz" | head -1)
    RIGHT_FASTQ=$(find "$DATASET_DIR" -maxdepth 1 -name "*R2*.fastq.gz" | head -1)
    if [ -z "$FASTA_FILE" ] || [ -z "$LEFT_FASTQ" ] || [ -z "$RIGHT_FASTQ" ]; then
        log_error "need a genome FASTA and R1/R2 .fastq.gz files in $DATASET_DIR"
        exit 1
    fi
    log_info "genome: $FASTA_FILE"
    log_info "reads:  $LEFT_FASTQ $RIGHT_FASTQ"
}

check_dependencies() {
    local missing=0
    for cmd in funannotate Trinity samtools hisat2; do
        if ! command -v "$cmd" > /dev/null 2>&1; then
            log_error "$cmd not in PATH"
            missing=1
        fi
    done
    [ "$missing" -eq 0 ] || exit 1
    # The Rust utilities must exist for the rust leg to mean anything.
    local home="${TRINITY_HOME:-${TRINITYHOME:-}}"
    if [ -z "$home" ] || [ ! -d "$home/rust_bio_utils/target/release" ]; then
        log_error "TRINITY_HOME/rust_bio_utils/target/release not found (TRINITY_HOME='$home')"
        exit 1
    fi
    TRINITY_RUST_DIR="$home/rust_bio_utils/target/release"
}

# Write the hardware and software that this run used, so legs can be checked
# for comparability afterwards.
record_environment() {
    local env_file="$RESULTS_DIR/${TEST_NAME}.env"
    {
        echo "engine=$ENGINE"
        echo "rep=$REP"
        echo "stage=$BENCH_STAGE"
        echo "cpus=$CPUS"
        echo "memory=$MEMORY"
        echo "host=$(hostname)"
        echo "cpu_model=$(lscpu | awk -F: '/Model name/ {gsub(/^ +/, "", $2); print $2; exit}')"
        echo "cpu_count_node=$(nproc --all)"
        echo "mem_total_kb=$(awk '/MemTotal/ {print $2}' /proc/meminfo)"
        echo "slurm_job_id=${SLURM_JOB_ID:-none}"
        echo "slurm_partition=${SLURM_JOB_PARTITION:-none}"
        echo "TRINITY_NO_RUST=${TRINITY_NO_RUST:-unset}"
        echo "trinity_version=$(Trinity --version 2>&1 | head -1)"
        echo "funannotate_version=$(funannotate --version 2>&1 | head -1)"
        echo "start=$(date -Is)"
    } > "$env_file"
}

# Count calls to Rust utilities in the run's logs.
count_rust_calls() {
    local out="$1"
    # grep exits 1 when nothing matches (the correct result for the perl leg),
    # which would end the script under pipefail; "|| true" keeps the count.
    { grep -rhs "rust_bio_utils/target/release" "$out/logfiles" "$out/training" --include='*.log' 2> /dev/null || true; } | wc -l
}

run_benchmark() {
    local out="$WORK_DIR/output"
    local test_log="$LOG_DIR/${TEST_NAME}.log"
    local test_metrics="$RESULTS_DIR/${TEST_NAME}.metrics"
    mkdir -p "$out"

    log_info "run: $TEST_NAME stage=$BENCH_STAGE cpus=$CPUS"
    if [ "$ENGINE" = "perl" ]; then
        export TRINITY_NO_RUST=1
    else
        unset TRINITY_NO_RUST
    fi
    log_info "TRINITY_NO_RUST=${TRINITY_NO_RUST:-unset}"

    local stage_args=()
    if [ "$BENCH_STAGE" = "trinity" ]; then
        stage_args+=(--stop_after_trinity)
    fi

    # The command is run as an if-condition so the script can still write the
    # exit status when funannotate fails.
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
            "${stage_args[@]}" \
            &> "$test_log" || rc=$?

    echo "end=$(date -Is)" >> "$RESULTS_DIR/${TEST_NAME}.env"
    echo "exit_status=$rc" >> "$test_metrics"
    if [ "$rc" -ne 0 ]; then
        log_error "funannotate train exited with $rc; see $test_log"
        echo "status=FAILED" >> "$test_metrics"
        return 1
    fi

    # Output summary, used later to check that the engines give the same result.
    # --stop_after_trinity leaves the transcripts in training/trinity.fasta; the
    # funannotate_train.trinity-GG.fasta name is made only by a full train run.
    local fasta="$out/training/funannotate_train.trinity-GG.fasta"
    if [ "$BENCH_STAGE" = "trinity" ]; then
        fasta="$out/training/trinity.fasta"
    fi
    if [ -s "$fasta" ]; then
        echo "transcripts=$(grep -c '^>' "$fasta")" >> "$test_metrics"
        echo "total_bases=$(grep -v '^>' "$fasta" | tr -d '\n' | wc -c)" >> "$test_metrics"
    else
        log_error "expected Trinity output not found: $fasta"
        echo "status=FAILED" >> "$test_metrics"
        return 1
    fi

    # Prove which engine ran.
    local calls
    calls=$(count_rust_calls "$out")
    echo "rust_calls_in_logs=$calls" >> "$test_metrics"
    if [ "$ENGINE" = "rust" ] && [ "$calls" -eq 0 ]; then
        log_error "rust leg made no Rust utility calls; the run is INVALID"
        echo "status=INVALID" >> "$test_metrics"
        return 4
    fi
    if [ "$ENGINE" = "perl" ] && [ "$calls" -ne 0 ]; then
        log_error "perl leg made $calls Rust utility calls; the run is INVALID"
        echo "status=INVALID" >> "$test_metrics"
        return 4
    fi
    echo "status=OK" >> "$test_metrics"
    log_info "OK: $TEST_NAME rust_calls=$calls"
    return 0
}

main() {
    log_info "Trinity $ENGINE benchmark, replicate $REP, stage $BENCH_STAGE"
    setup_directories
    validate_dataset
    check_dependencies
    record_environment
    run_benchmark
}

main "$@"
