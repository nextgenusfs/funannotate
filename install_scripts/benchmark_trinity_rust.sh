#!/usr/bin/env bash
# Benchmark funannotate's Trinity workflow with a single engine (rust or perl).
#
# This script runs `funannotate train` once, either with the Rust Trinity
# utilities on PATH or with them hidden (falls back to Perl), and collects
# timing/memory metrics for that one run. Run it twice (rust and perl) --
# as two independent SLURM jobs via benchmark_trinity_rust_rust.sbatch /
# benchmark_trinity_rust_perl.sbatch -- so each gets its own walltime budget,
# then use scripts/compare_trinity_rust_results.sh to merge the two into a
# comparison report.
#
# Usage:
#   ./scripts/benchmark_trinity_rust.sh [dataset_dir] [output_dir] <rust|perl> [cpus]
#
# The optional [cpus] argument controls the thread count passed to
# `funannotate train --cpus`. If omitted it defaults to the number of CPUs
# detected on the host (nproc), so a SLURM launcher can pass in
# $SLURM_CPUS_PER_TASK instead of relying on a hardcoded value.
#
# Example:
#   ./scripts/benchmark_trinity_rust.sh \
#       ~/projects/funannotate/trinity_example/Cordyceps_militaris \
#       ~/benchmark_results \
#       rust \
#       48

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

DATASET_DIR="${1:-.}"
BENCHMARK_OUT_DIR="${2:-.}"
ENGINE="${3:?Usage: $0 <dataset_dir> <benchmark_out_dir> <rust|perl> [cpus]}"
# Number of CPUs to hand to `funannotate train --cpus`. Defaults to the host's
# detected core count so it can be driven by the launching environment (e.g.
# $SLURM_CPUS_PER_TASK) rather than a hardcoded value.
CPUS="${4:-$(nproc)}"
case "$ENGINE" in
    rust) TEST_NAME="rust_enabled"; USE_RUST="true" ;;
    perl) TEST_NAME="perl_only"; USE_RUST="false" ;;
    *) echo "ERROR: engine must be 'rust' or 'perl', got: $ENGINE" >&2; exit 1 ;;
esac

BENCHMARK_TIME="/usr/bin/time -v"
# logs/results are small and precious -- always write them directly to
# BENCHMARK_OUT_DIR (typically /bigdata, a network filesystem) so nothing is
# lost if the job is killed. The working directory is the opposite: Trinity-GG
# churns through thousands of small per-partition files (trinity_gg/ alone ran
# ~6GB of small files in a prior run) which are painfully slow on a network
# filesystem, so that goes on node-local $SCRATCH instead and only the final
# artifacts get synced back (see sync_work_dir_back).
LOG_DIR="${BENCHMARK_OUT_DIR}/logs"
RESULTS_DIR="${BENCHMARK_OUT_DIR}/results"
SCRATCH_BASE="${SCRATCH:-$BENCHMARK_OUT_DIR}"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Helper functions
log_info() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $*"
}

log_success() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')] ✓ $*${NC}"
}

log_error() {
    echo -e "${RED}[$(date +'%Y-%m-%d %H:%M:%S')] ✗ $*${NC}" >&2
}

log_warning() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')] ⚠ $*${NC}"
}

# Setup
setup_directories() {
    log_info "Setting up benchmark directories..."
    mkdir -p "$LOG_DIR" "$RESULTS_DIR"

    # Working directory lives on scratch (see comment above); final artifacts
    # get synced back to this durable location on BENCHMARK_OUT_DIR.
    SYNC_TARGET_DIR="$BENCHMARK_OUT_DIR/work_${ENGINE}"
    WORK_DIR="$SCRATCH_BASE/funannotate_trinity_bench_${ENGINE}"
    mkdir -p "$WORK_DIR" "$SYNC_TARGET_DIR"

    log_success "Working directory (scratch): $WORK_DIR"
    log_success "Sync target (durable): $SYNC_TARGET_DIR"

    # Sync back on normal exit, error exit, or SLURM killing the job for
    # hitting its walltime (SIGTERM) -- so a timeout still leaves whatever
    # progress was made instead of losing it all on scratch.
    trap 'sync_work_dir_back' EXIT
    trap 'log_warning "Received termination signal (e.g. SLURM walltime) -- syncing partial results before exit..."; exit 143' TERM INT
}

# Copy the final training/ artifacts back to durable storage, skipping
# Trinity's disposable per-partition working directory (trinity_gg/ -- by far
# the largest and most numerous set of files, and never read by downstream
# `funannotate predict`). Registered as a trap so it runs on success, failure,
# or SLURM killing the job for hitting its walltime.
sync_work_dir_back() {
    if [ -d "$WORK_DIR/output" ]; then
        log_info "Syncing results from scratch back to $SYNC_TARGET_DIR..."
        mkdir -p "$SYNC_TARGET_DIR"
        rsync -a --exclude='output/training/trinity_gg' "$WORK_DIR/output" "$SYNC_TARGET_DIR/" \
            && log_success "Sync complete" \
            || log_warning "Sync back to $SYNC_TARGET_DIR encountered errors (partial results may be present)"
    fi
}

validate_dataset() {
    log_info "Validating input dataset..."

    if [ ! -d "$DATASET_DIR" ]; then
        log_error "Dataset directory not found: $DATASET_DIR"
        exit 1
    fi

    # Find FASTA file
    FASTA_FILE=$(find "$DATASET_DIR" -maxdepth 1 -name "*.fasta" -o -name "*.fa" | head -1)
    if [ -z "$FASTA_FILE" ]; then
        log_error "No FASTA genome file found in $DATASET_DIR"
        exit 1
    fi

    # Find fastq files
    LEFT_FASTQ=$(find "$DATASET_DIR" -maxdepth 1 -name "*R1*.fastq.gz" | head -1)
    RIGHT_FASTQ=$(find "$DATASET_DIR" -maxdepth 1 -name "*R2*.fastq.gz" | head -1)

    if [ -z "$LEFT_FASTQ" ] || [ -z "$RIGHT_FASTQ" ]; then
        log_error "Could not find paired-end FASTQ files in $DATASET_DIR"
        exit 1
    fi

    log_success "Found dataset files:"
    echo "  Genome:      $(basename "$FASTA_FILE") ($(du -h "$FASTA_FILE" | cut -f1))"
    echo "  Left reads:  $(basename "$LEFT_FASTQ") ($(du -h "$LEFT_FASTQ" | cut -f1))"
    echo "  Right reads: $(basename "$RIGHT_FASTQ") ($(du -h "$RIGHT_FASTQ" | cut -f1))"
}

check_dependencies() {
    log_info "Checking dependencies..."

    local missing=0

    for cmd in funannotate Trinity samtools hisat2; do
        if ! command -v "$cmd" &> /dev/null; then
            log_warning "  ✗ $cmd not in PATH"
            missing=$((missing + 1))
        else
            log_info "  ✓ $cmd found"
        fi
    done

    if [ $missing -gt 0 ]; then
        log_error "Some dependencies are missing. Ensure pixi environment is activated."
        exit 1
    fi

    # Check for Rust Trinity tools
    local rust_tools=("sam_to_read_coords" "extract_reads_per_partition" "fragment_coverage_writer" "define_coverage_partitions")
    local rust_available=0

    for tool in "${rust_tools[@]}"; do
        if command -v "$tool" &> /dev/null; then
            rust_available=$((rust_available + 1))
        fi
    done

    if [ $rust_available -eq ${#rust_tools[@]} ]; then
        log_success "All Rust Trinity tools found ($rust_available/${#rust_tools[@]})"
    else
        log_warning "Only $rust_available/${#rust_tools[@]} Rust Trinity tools found"
    fi
}

get_system_info() {
    log_info "System information:"
    echo "  CPU cores: $(nproc)"
    echo "  Memory: $(free -h | grep Mem | awk '{print $2}')"
    echo "  Available: $(free -h | grep Mem | awk '{print $7}')"

    if command -v uname &> /dev/null; then
        echo "  OS: $(uname -s) $(uname -r)"
    fi
}

run_benchmark() {
    local test_name="$1"
    local work_dir="$2"
    local use_rust="$3"

    local test_log="$LOG_DIR/${test_name}.log"
    local test_metrics="$RESULTS_DIR/${test_name}.metrics"

    log_info "Starting $test_name benchmark..."
    log_info "  Output: $work_dir"
    log_info "  Using Rust tools: $use_rust"
    log_info "  CPUs: $CPUS"

    # Setup PATH
    local PATH_TO_USE="$PATH"
    if [ "$use_rust" == "false" ]; then
        # Remove Rust Trinity tools from PATH by creating a temp directory without them
        local temp_bin_dir=$(mktemp -d)
        trap "rm -rf '$temp_bin_dir'" RETURN

        # Copy all bin contents except Rust Trinity tools
        for file in $(dirname "$(command -v funannotate)")/*; do
            local basename=$(basename "$file")
            # Skip Rust Trinity tools
            if [[ ! "$basename" =~ ^(sam_to_read_coords|extract_reads_per_partition|fragment_coverage_writer|define_coverage_partitions)$ ]]; then
                ln -s "$file" "$temp_bin_dir/$basename" 2>/dev/null || true
            fi
        done
        PATH_TO_USE="$temp_bin_dir:$(dirname "$(command -v samtools)"):/usr/bin:/bin"
    fi

    # Create output directory structure for funannotate
    local output_base="$work_dir/output"
    mkdir -p "$output_base"

    local start_time=$(date +%s.%N)
    local start_memory=$(grep MemAvailable /proc/meminfo | awk '{print $2}')

    # Run funannotate train with the test configuration.
    # Thread count comes from $CPUS (see arg parsing above) so it tracks the
    # CPUs the launching environment actually allocated.
    if PATH="$PATH_TO_USE" $BENCHMARK_TIME -o "$test_metrics" \
        funannotate train \
            -i "$FASTA_FILE" \
            -o "$output_base" \
            --left "$LEFT_FASTQ" \
            --right "$RIGHT_FASTQ" \
            --left_norm "$LEFT_FASTQ" \
            --right_norm "$RIGHT_FASTQ" \
            --species "Cordyceps militaris benchmark" \
            --cpus "$CPUS" \
            --memory 16G \
            --jaccard_clip \
            &> "$test_log"; then

        local end_time=$(date +%s.%N)
        local end_memory=$(grep MemAvailable /proc/meminfo | awk '{print $2}')
        local wall_time=$(echo "$end_time - $start_time" | bc)
        local memory_used=$(echo "($start_memory - $end_memory) / 1024" | bc)

        log_success "$test_name completed in $wall_time seconds"

        # Extract Trinity output stats
        if [ -f "$output_base/transcripts.fasta" ]; then
            local num_transcripts=$(grep -c "^>" "$output_base/transcripts.fasta")
            local total_bases=$(grep -v "^>" "$output_base/transcripts.fasta" | wc -c)
            echo "transcripts=$num_transcripts" >> "$test_metrics"
            echo "total_bases=$total_bases" >> "$test_metrics"
            log_info "  Transcripts assembled: $num_transcripts"
            log_info "  Total bases: $total_bases"
        fi

        return 0
    else
        log_error "$test_name failed"
        log_info "See $test_log for details"
        return 1
    fi
}

main() {
    log_info "Trinity ${ENGINE^} Engine Benchmark"
    log_info "===================================="
    echo ""

    setup_directories
    validate_dataset
    check_dependencies
    echo ""
    get_system_info
    echo ""

    log_info "Running benchmark..."
    echo ""

    if run_benchmark "$TEST_NAME" "$WORK_DIR" "$USE_RUST"; then
        log_success "${ENGINE^} benchmark completed"
    else
        log_warning "${ENGINE^} benchmark encountered issues"
        log_info "Results directory: $RESULTS_DIR"
        log_info "Logs directory: $LOG_DIR"
        exit 1
    fi

    echo ""
    log_info "Benchmark complete!"
    log_info "Metrics: $RESULTS_DIR/${TEST_NAME}.metrics"
    log_info "Log: $LOG_DIR/${TEST_NAME}.log"
    echo ""
    log_info "Once both rust and perl runs have completed, run:"
    log_info "  scripts/compare_trinity_rust_results.sh \"$BENCHMARK_OUT_DIR\""
}

# Run main function
main "$@"
