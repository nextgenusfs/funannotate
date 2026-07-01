#!/usr/bin/env bash
# Benchmark funannotate predict with the Rust EvidenceModeler engine vs. the Perl engine
#
# EVM (not PASA) is the component actually exercised live during `funannotate predict`
# -- PASA only contributes a pre-computed GFF3 (funannotate_train.pasa.gff3) generated
# during `funannotate train`, it does not run again during predict. The engine switch
# is controlled entirely by the FUNANNOTATE_EVM_ENGINE env var (rust|perl); see
# funannotate/predict.py around line 398.
#
# This script:
# 1. Copies a completed `funannotate train` output directory (it must contain a
#    training/ subfolder with funannotate_train.pasa.gff3, coordSorted.bam, etc.)
#    into two independent working copies.
# 2. Runs `funannotate predict` on each copy: one with FUNANNOTATE_EVM_ENGINE=rust
#    (the evidence_modeler Rust binary must be on PATH), one with
#    FUNANNOTATE_EVM_ENGINE=perl (requires a working Perl EVM_HOME).
# 3. Collects timing/memory metrics and produces a comparison report.
#
# Prerequisites:
#   - $FUNANNOTATE_DB set and `funannotate setup` already run (predict defaults
#     --protein_evidence to $FUNANNOTATE_DB/uniprot_sprot.fasta)
#   - A Perl EvidenceModeler checkout available for the "perl" run (EvmUtils/
#     partition_EVM_inputs.pl must exist under it) -- pass via --evm-home or
#     the PERL_EVM_HOME env var.
#   - evidence_modeler (Rust) available on PATH for the "rust" run (installed by
#     scripts/pixi_install_rust_evm.sh as part of pixi activation).
#
# Usage:
#   ./scripts/benchmark_predict_evm.sh <train_output_dir> [benchmark_out_dir] [evm_home]
#
# Example:
#   ./scripts/benchmark_predict_evm.sh \
#       ~/bench_results/work_rust/output \
#       ~/predict_bench_results \
#       ~/projects/funannotate/EVidenceModeler

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

TRAIN_DIR="${1:?Usage: $0 <train_output_dir> [benchmark_out_dir] [evm_home]}"
BENCHMARK_OUT_DIR="${2:-.}"
PERL_EVM_HOME="${3:-${PERL_EVM_HOME:-}}"

SPECIES="${PREDICT_SPECIES:-Cordyceps militaris benchmark}"
AUGUSTUS_SPECIES="${PREDICT_AUGUSTUS_SPECIES:-fusarium_graminearum}"
BENCHMARK_TIME="/usr/bin/time -v"
LOG_DIR="${BENCHMARK_OUT_DIR}/logs"
RESULTS_DIR="${BENCHMARK_OUT_DIR}/results"

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

    RUST_WORK_DIR="$BENCHMARK_OUT_DIR/predict_work_rust"
    PERL_WORK_DIR="$BENCHMARK_OUT_DIR/predict_work_perl"

    log_success "Created working directories"
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

    if command -v evidence_modeler &> /dev/null; then
        log_success "Rust evidence_modeler found ($(command -v evidence_modeler))"
    else
        log_error "Rust evidence_modeler not found in PATH -- cannot run rust_enabled benchmark"
        exit 1
    fi

    if [ -z "$PERL_EVM_HOME" ]; then
        log_error "No Perl EVM_HOME given (arg 3, or \$PERL_EVM_HOME) -- cannot run perl_only benchmark"
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
            --cpus 4 \
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

compare_results() {
    log_info "Generating comparison report..."

    local report_file="$RESULTS_DIR/PREDICT_BENCHMARK_REPORT.md"
    local timestamp
    timestamp=$(date +'%Y-%m-%d %H:%M:%S')

    {
        echo "# Funannotate Predict EVM Engine Benchmark Report"
        echo ""
        echo "**Generated**: $timestamp"
        echo ""
        echo "## Dataset"
        echo "- **Train output**: $TRAIN_DIR"
        echo "- **Genome**: $(basename "$GENOME_FASTA")"
        echo ""
        echo "## Benchmark Results"
        echo ""
        echo "### Raw Metrics"
        echo ""
        echo "#### Rust EVM engine"
        echo '```'
        cat "$RESULTS_DIR/rust_enabled.metrics" 2>/dev/null | grep -E "^(Elapsed|Maximum resident|User|System|Minor faults|Major faults|Command|genes)" || echo "Metrics not available"
        echo '```'
        echo ""
        echo "#### Perl EVM engine"
        echo '```'
        cat "$RESULTS_DIR/perl_only.metrics" 2>/dev/null | grep -E "^(Elapsed|Maximum resident|User|System|Minor faults|Major faults|Command|genes)" || echo "Metrics not available"
        echo '```'
        echo ""
        echo "### Performance Comparison"
        echo ""
        if [ -f "$RESULTS_DIR/rust_enabled.metrics" ] && [ -f "$RESULTS_DIR/perl_only.metrics" ]; then
            echo "| Metric | Rust | Perl | Speedup |"
            echo "|--------|------|------|---------|"

            local rust_elapsed perl_elapsed
            rust_elapsed=$(grep "Elapsed" "$RESULTS_DIR/rust_enabled.metrics" | head -1 | awk '{print $NF}' | sed 's/[:,]//g')
            perl_elapsed=$(grep "Elapsed" "$RESULTS_DIR/perl_only.metrics" | head -1 | awk '{print $NF}' | sed 's/[:,]//g')
            if [ -n "$rust_elapsed" ] && [ -n "$perl_elapsed" ]; then
                local speedup
                speedup=$(echo "scale=2; $perl_elapsed / $rust_elapsed" | bc)
                echo "| Elapsed Time | $rust_elapsed | $perl_elapsed | ${speedup}x |"
            fi

            local rust_mem perl_mem
            rust_mem=$(grep "Maximum resident" "$RESULTS_DIR/rust_enabled.metrics" | awk '{print $(NF-1)}')
            perl_mem=$(grep "Maximum resident" "$RESULTS_DIR/perl_only.metrics" | awk '{print $(NF-1)}')
            if [ -n "$rust_mem" ] && [ -n "$perl_mem" ]; then
                echo "| Peak Memory | ${rust_mem}K | ${perl_mem}K | N/A |"
            fi
        fi
        echo ""
        echo "## Logs"
        echo ""
        echo "- Rust run: \`$LOG_DIR/rust_enabled.log\`"
        echo "- Perl run: \`$LOG_DIR/perl_only.log\`"
        echo ""
        echo "## Notes"
        echo ""
        echo "- Both runs used 4 CPUs, augustus_species=$AUGUSTUS_SPECIES (a generic fungal stand-in;"
        echo "  this branch of funannotate train does not do BUSCO-based ab initio training)"
        echo "- Ab initio prediction (Augustus/GeneMark) work is identical between runs; only EVM"
        echo "  consensus building differs, so total wall time includes shared, engine-independent overhead"
    } > "$report_file"

    log_success "Benchmark report generated: $report_file"
}

main() {
    log_info "Funannotate Predict EVM Engine Benchmark"
    log_info "========================================="
    echo ""

    setup_directories
    validate_train_dir
    check_dependencies
    echo ""

    log_info "Running benchmarks..."
    echo ""

    if run_predict "rust_enabled" "$RUST_WORK_DIR" "rust"; then
        log_success "Rust EVM benchmark completed"
    else
        log_warning "Rust EVM benchmark encountered issues"
    fi

    echo ""

    if run_predict "perl_only" "$PERL_WORK_DIR" "perl"; then
        log_success "Perl EVM benchmark completed"
    else
        log_warning "Perl EVM benchmark encountered issues"
    fi

    echo ""
    compare_results

    echo ""
    log_info "Benchmark complete!"
    log_info "Results directory: $RESULTS_DIR"

    echo ""
    echo -e "${GREEN}========== BENCHMARK SUMMARY ==========${NC}"
    cat "$RESULTS_DIR/PREDICT_BENCHMARK_REPORT.md" || true
}

main "$@"
