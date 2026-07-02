#!/usr/bin/env bash
# Benchmark funannotate predict with the Rust EvidenceModeler engine vs. the Perl engine
#
# EVM (not PASA) is the component actually exercised live during `funannotate predict`
# -- PASA only contributes a pre-computed GFF3 (funannotate_train.pasa.gff3) generated
# during `funannotate train`, it does not run again during predict. The engine switch
# is controlled entirely by the FUNANNOTATE_EVM_ENGINE env var (rust|perl); see
# funannotate/predict.py around line 398.
#
# This script runs `funannotate predict` once, against a copy of a completed
# `funannotate train` output directory, with a single EVM engine (rust or perl).
# Run it twice (rust and perl) -- as two independent SLURM jobs via
# benchmark_predict_evm_rust.sbatch / benchmark_predict_evm_perl.sbatch -- so
# each gets its own walltime budget, then use
# scripts/compare_predict_evm_results.sh to merge the two into a report.
#
# Prerequisites:
#   - $FUNANNOTATE_DB set and `funannotate setup` already run (predict defaults
#     --protein_evidence to $FUNANNOTATE_DB/uniprot_sprot.fasta)
#   - A Perl EvidenceModeler checkout (EvmUtils/partition_EVM_inputs.pl must
#     exist under it) -- pass via arg 4 or the PERL_EVM_HOME env var. Needed
#     for BOTH engines: predict.py calls Perl EvmUtils/misc converter scripts
#     unconditionally, regardless of which engine builds the consensus models.
#   - evidence_modeler (Rust) available on PATH for the "rust" run (installed by
#     scripts/pixi_install_evm_rust.sh as part of pixi activation).
#
# Usage:
#   ./scripts/benchmark_predict_evm.sh <train_output_dir> <benchmark_out_dir> <rust|perl> [evm_home] [cpus]
#
# The optional [cpus] argument controls the thread count passed to
# `funannotate predict --cpus`. If omitted it defaults to the number of CPUs
# detected on the host (nproc), so a SLURM launcher can pass in
# $SLURM_CPUS_PER_TASK instead of relying on a hardcoded value.
#
# Example:
#   ./scripts/benchmark_predict_evm.sh \
#       /bigdata/stajichlab/jstajich/projects/funannotate/BENCHMARK/trinity_rust/work_rust/output \
#       /bigdata/stajichlab/jstajich/projects/funannotate/BENCHMARK/predict_evm \
#       rust \
#       ~/projects/funannotate/EVidenceModeler \
#       48

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

TRAIN_DIR="${1:?Usage: $0 <train_output_dir> <benchmark_out_dir> <rust|perl> [evm_home]}"
BENCHMARK_OUT_DIR="${2:-.}"
ENGINE="${3:?Usage: $0 <train_output_dir> <benchmark_out_dir> <rust|perl> [evm_home]}"
case "$ENGINE" in
    rust) TEST_NAME="rust_enabled" ;;
    perl) TEST_NAME="perl_only" ;;
    *) echo "ERROR: engine must be 'rust' or 'perl', got: $ENGINE" >&2; exit 1 ;;
esac
PERL_EVM_HOME="${4:-${PERL_EVM_HOME:-}}"
# Number of CPUs to hand to `funannotate predict --cpus`. Defaults to the
# host's detected core count so it can be driven by the launching environment
# (e.g. $SLURM_CPUS_PER_TASK) rather than a hardcoded value.
CPUS="${5:-$(nproc)}"

SPECIES="${PREDICT_SPECIES:-Cordyceps militaris benchmark}"
AUGUSTUS_SPECIES="${PREDICT_AUGUSTUS_SPECIES:-fusarium_graminearum}"
BENCHMARK_TIME="/usr/bin/time -v"
# logs/results are small and precious -- always write them directly to
# BENCHMARK_OUT_DIR (typically /bigdata, a network filesystem) so nothing is
# lost if the job is killed. The working copy (predict churns through many
# per-contig/per-partition Augustus/GeneMark/EVM files) goes on node-local
# $SCRATCH instead, synced back to BENCHMARK_OUT_DIR only once at the end.
LOG_DIR="${BENCHMARK_OUT_DIR}/logs"
RESULTS_DIR="${BENCHMARK_OUT_DIR}/results"
SCRATCH_BASE="${SCRATCH:-$BENCHMARK_OUT_DIR}"

# Color output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

log_info() {
    echo -e "${BLUE}[$(date +'%Y-%m-%d %H:%M:%S')]${NC} $*"
}

log_success() {
    echo -e "${GREEN}[$(date +'%Y-%m-%d %H:%M:%S')] \xe2\x9c\x93 $*${NC}"
}

log_error() {
    echo -e "${RED}[$(date +'%Y-%m-%d %H:%M:%S')] \xe2\x9c\x97 $*${NC}" >&2
}

log_warning() {
    echo -e "${YELLOW}[$(date +'%Y-%m-%d %H:%M:%S')] \xe2\x9a\xa0 $*${NC}"
}

setup_directories() {
    log_info "Setting up benchmark directories..."
    mkdir -p "$LOG_DIR" "$RESULTS_DIR"

    SYNC_TARGET_DIR="$BENCHMARK_OUT_DIR/predict_work_${ENGINE}"
    WORK_DIR="$SCRATCH_BASE/funannotate_predict_bench_${ENGINE}"
    mkdir -p "$SYNC_TARGET_DIR"

    log_success "Working directory (scratch): $WORK_DIR"
    log_success "Sync target (durable): $SYNC_TARGET_DIR"

    # Sync back on normal exit, error exit, or SLURM killing the job for
    # hitting its walltime (SIGTERM) -- so a timeout still leaves whatever
    # progress was made instead of losing it all on scratch.
    trap 'sync_work_dir_back' EXIT
    trap 'log_warning "Received termination signal (e.g. SLURM walltime) -- syncing partial results before exit..."; exit 143' TERM INT
}

# Copy the predict working directory back to durable storage. No excludes
# yet (unlike the Trinity benchmark's trinity_gg/ exclusion) -- add one here
# once a real run shows which subdirectory is the disposable bulk.
sync_work_dir_back() {
    if [ -d "$WORK_DIR" ]; then
        log_info "Syncing results from scratch back to $SYNC_TARGET_DIR..."
        mkdir -p "$SYNC_TARGET_DIR"
        rsync -a "$WORK_DIR/" "$SYNC_TARGET_DIR/" \
            && log_success "Sync complete" \
            || log_warning "Sync back to $SYNC_TARGET_DIR encountered errors (partial results may be present)"
    fi
}

validate_train_dir() {
    log_info "Validating training directory..."

    if [ ! -d "$TRAIN_DIR" ]; then
        log_error "Train output directory not found: $TRAIN_DIR"
        exit 1
    fi

    if [ ! -d "$TRAIN_DIR/training" ]; then
        log_error "$TRAIN_DIR does not look like a funannotate train output (no training/ subfolder)"
        exit 1
    fi

    GENOME_FASTA="$TRAIN_DIR/training/genome.fasta"
    if [ ! -f "$GENOME_FASTA" ]; then
        log_error "Genome fasta not found at $GENOME_FASTA"
        exit 1
    fi

    if [ ! -f "$TRAIN_DIR/training/funannotate_train.pasa.gff3" ]; then
        log_warning "No funannotate_train.pasa.gff3 found -- predict will run without PASA transcript evidence"
    fi

    log_success "Training directory looks valid: $TRAIN_DIR"
}

check_dependencies() {
    log_info "Checking dependencies..."

    local missing=0
    for cmd in funannotate augustus tbl2asn diamond; do
        if ! command -v "$cmd" &> /dev/null; then
            log_warning "  \xe2\x9c\x97 $cmd not in PATH"
            missing=$((missing + 1))
        else
            log_info "  \xe2\x9c\x93 $cmd found"
        fi
    done
    if [ $missing -gt 0 ]; then
        log_error "Some dependencies are missing. Ensure pixi environment is activated."
        exit 1
    fi

    if [ -z "${FUNANNOTATE_DB:-}" ]; then
        log_error "\$FUNANNOTATE_DB is not set. predict needs it (e.g. for --protein_evidence default)."
        exit 1
    fi
    log_info "  \xe2\x9c\x93 FUNANNOTATE_DB=$FUNANNOTATE_DB"

    if [ "$ENGINE" == "rust" ]; then
        if command -v evidence_modeler &> /dev/null; then
            log_success "Rust evidence_modeler found ($(command -v evidence_modeler))"
        else
            log_error "Rust evidence_modeler not found in PATH -- cannot run rust_enabled benchmark"
            exit 1
        fi
    fi

    # Needed for both engines: predict.py calls Perl EvmUtils/misc converter
    # scripts unconditionally, regardless of which engine builds consensus models.
    if [ -z "$PERL_EVM_HOME" ]; then
        log_error "No Perl EVM_HOME given (arg 4, or \$PERL_EVM_HOME)"
        exit 1
    fi
    if [ ! -f "$PERL_EVM_HOME/EvmUtils/partition_EVM_inputs.pl" ]; then
        log_error "Perl EVM_HOME invalid: $PERL_EVM_HOME/EvmUtils/partition_EVM_inputs.pl not found"
        exit 1
    fi
    log_success "Perl EVM_HOME validated: $PERL_EVM_HOME"
}

run_predict() {
    local test_name="$1"
    local work_dir="$2"
    local engine="$3"

    local test_log="$LOG_DIR/${test_name}.log"
    local test_metrics="$RESULTS_DIR/${test_name}.metrics"

    log_info "Starting $test_name benchmark..."
    log_info "  Output: $work_dir"
    log_info "  EVM engine: $engine"
    log_info "  CPUs: $CPUS"

    rm -rf "$work_dir"
    mkdir -p "$(dirname "$work_dir")"
    log_info "  Copying training directory (this can take a minute)..."
    cp -r "$TRAIN_DIR" "$work_dir"

    # NOTE: predict.py calls perl $EVM_HOME/EvmUtils/misc/augustus_GFF3_to_EVM_GFF3.pl
    # (and the GTF variant) unconditionally, regardless of engine -- these Perl helper
    # scripts were never ported into EVidenceModeler_rust (it ships bin/augustus_to_evm_gff3
    # instead, which predict.py doesn't call). So EVM_HOME must point at a real Perl EVM
    # checkout for BOTH the rust and perl runs, or Augustus/GeneMark GFF conversion fails.
    local env_args=(FUNANNOTATE_EVM_ENGINE="$engine" EVM_HOME="$PERL_EVM_HOME")
    local extra_args=(--EVM_HOME "$PERL_EVM_HOME")

    local start_time
    start_time=$(date +%s.%N)

    if env "${env_args[@]}" $BENCHMARK_TIME -o "$test_metrics" \
        funannotate predict \
            -i "$work_dir/training/genome.fasta" \
            -o "$work_dir" \
            -s "$SPECIES" \
            --augustus_species "$AUGUSTUS_SPECIES" \
            --cpus "$CPUS" \
            "${extra_args[@]}" \
            &> "$test_log"; then

        local end_time
        end_time=$(date +%s.%N)
        local wall_time
        wall_time=$(echo "$end_time - $start_time" | bc)
        log_success "$test_name completed in $wall_time seconds"

        local gff3_file
        gff3_file=$(find "$work_dir/predict_results" -maxdepth 1 -name "*.gff3" 2>/dev/null | head -1)
        if [ -n "$gff3_file" ]; then
            local num_genes
            num_genes=$(grep -c $'\tgene\t' "$gff3_file" 2>/dev/null || echo "NA")
            echo "genes=$num_genes" >> "$test_metrics"
            log_info "  Gene models predicted: $num_genes"
        fi

        return 0
    else
        log_error "$test_name failed"
        log_info "See $test_log for details"
        return 1
    fi
}

main() {
    log_info "Funannotate Predict ${ENGINE^} EVM Engine Benchmark"
    log_info "========================================="
    echo ""

    setup_directories
    validate_train_dir
    check_dependencies
    echo ""

    log_info "Running benchmark..."
    echo ""

    if run_predict "$TEST_NAME" "$WORK_DIR" "$ENGINE"; then
        log_success "${ENGINE^} EVM benchmark completed"
    else
        log_warning "${ENGINE^} EVM benchmark encountered issues"
        log_info "Results directory: $RESULTS_DIR"
        exit 1
    fi

    echo ""
    log_info "Benchmark complete!"
    log_info "Metrics: $RESULTS_DIR/${TEST_NAME}.metrics"
    log_info "Log: $LOG_DIR/${TEST_NAME}.log"
    echo ""
    log_info "Once both rust and perl runs have completed, run:"
    log_info "  scripts/compare_predict_evm_results.sh \"$BENCHMARK_OUT_DIR\""
}

main "$@"
