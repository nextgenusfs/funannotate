#!/usr/bin/env bash
# Benchmark funannotate Trinity workflow with and without Rust optimizations
#
# This script tests Trinity's RNA-seq assembly performance using Rust-optimized
# utilities vs. traditional Perl-based versions. It:
#
# 1. Runs funannotate train with Rust Trinity tools in PATH
# 2. Runs the same workflow with Rust tools hidden (falls back to Perl)
# 3. Collects timing, memory, and CPU metrics
# 4. Generates a comprehensive comparison report
#
# Usage:
#   ./scripts/benchmark_trinity_rust.sh [dataset_dir] [output_dir]
#
# Example:
#   ./scripts/benchmark_trinity_rust.sh \
#       ~/projects/funannotate/trinity_example/Cordyceps_militaris \
#       ~/benchmark_results

set -euo pipefail

# Configuration
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

DATASET_DIR="${1:-.}"
BENCHMARK_OUT_DIR="${2:-.}"
BENCHMARK_TIME="/usr/bin/time -v"
LOG_DIR="${BENCHMARK_OUT_DIR}/logs"
RESULTS_DIR="${BENCHMARK_OUT_DIR}/results"

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

    # Create temporary working directories
    RUST_WORK_DIR="$BENCHMARK_OUT_DIR/work_rust"
    PERL_WORK_DIR="$BENCHMARK_OUT_DIR/work_perl"
    mkdir -p "$RUST_WORK_DIR" "$PERL_WORK_DIR"

    log_success "Created working directories"
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
            ((missing++))
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
            ((rust_available++))
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

    # Run funannotate train with the test configuration
    # Using --cpus 4 to keep runtime reasonable for benchmarking
    if PATH="$PATH_TO_USE" $BENCHMARK_TIME -o "$test_metrics" \
        funannotate train \
            -i "$FASTA_FILE" \
            -o "$output_base" \
            --left_norm "$LEFT_FASTQ" \
            --right_norm "$RIGHT_FASTQ" \
            -s "Cordyceps_militaris_benchmark" \
            --cpus 4 \
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

compare_results() {
    log_info "Generating comparison report..."

    local report_file="$RESULTS_DIR/BENCHMARK_REPORT.md"
    local timestamp=$(date +'%Y-%m-%d %H:%M:%S')

    {
        echo "# Trinity Rust Optimization Benchmark Report"
        echo ""
        echo "**Generated**: $timestamp"
        echo ""
        echo "## Dataset"
        echo "- **Genome**: $(basename "$FASTA_FILE")"
        echo "- **Left reads**: $(basename "$LEFT_FASTQ") ($(du -h "$LEFT_FASTQ" | cut -f1))"
        echo "- **Right reads**: $(basename "$RIGHT_FASTQ") ($(du -h "$RIGHT_FASTQ" | cut -f1))"
        echo ""
        echo "## Benchmark Results"
        echo ""
        echo "### Raw Metrics"
        echo ""
        echo "#### With Rust Optimization"
        echo "\`\`\`"
        cat "$RESULTS_DIR/rust_enabled.metrics" 2>/dev/null | grep -E "^(Elapsed|Maximum resident|User|System|Minor faults|Major faults|Command|transcripts|total_bases)" || echo "Metrics not available"
        echo "\`\`\`"
        echo ""
        echo "#### With Perl (Rust tools disabled)"
        echo "\`\`\`"
        cat "$RESULTS_DIR/perl_only.metrics" 2>/dev/null | grep -E "^(Elapsed|Maximum resident|User|System|Minor faults|Major faults|Command|transcripts|total_bases)" || echo "Metrics not available"
        echo "\`\`\`"
        echo ""
        echo "### Performance Comparison"
        echo ""

        # Try to extract and compare metrics
        if [ -f "$RESULTS_DIR/rust_enabled.metrics" ] && [ -f "$RESULTS_DIR/perl_only.metrics" ]; then
            echo "| Metric | Rust | Perl | Speedup |"
            echo "|--------|------|------|---------|"

            local rust_elapsed=$(grep "Elapsed" "$RESULTS_DIR/rust_enabled.metrics" | head -1 | awk '{print $NF}' | sed 's/[:,]//g')
            local perl_elapsed=$(grep "Elapsed" "$RESULTS_DIR/perl_only.metrics" | head -1 | awk '{print $NF}' | sed 's/[:,]//g')

            if [ -n "$rust_elapsed" ] && [ -n "$perl_elapsed" ]; then
                local speedup=$(echo "scale=2; $perl_elapsed / $rust_elapsed" | bc)
                echo "| Elapsed Time | $rust_elapsed | $perl_elapsed | ${speedup}x |"
            fi

            local rust_mem=$(grep "Maximum resident" "$RESULTS_DIR/rust_enabled.metrics" | awk '{print $(NF-1)}')
            local perl_mem=$(grep "Maximum resident" "$RESULTS_DIR/perl_only.metrics" | awk '{print $(NF-1)}')

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
        echo "## Working Directories"
        echo ""
        echo "- Rust output: \`$RUST_WORK_DIR/output\`"
        echo "- Perl output: \`$PERL_WORK_DIR/output\`"
        echo ""
        echo "## Notes"
        echo ""
        echo "- Both runs used 4 CPUs and 16GB memory limit"
        echo "- Trinity may produce slightly different transcripts due to non-deterministic assembly"
        echo "- Focus on elapsed time differences for performance comparison"

    } > "$report_file"

    log_success "Benchmark report generated: $report_file"
}

main() {
    log_info "Trinity Rust Optimization Benchmark"
    log_info "===================================="
    echo ""

    setup_directories
    validate_dataset
    check_dependencies
    echo ""
    get_system_info
    echo ""

    # Run benchmarks
    log_info "Running benchmarks..."
    echo ""

    if run_benchmark "rust_enabled" "$RUST_WORK_DIR" "true"; then
        log_success "Rust optimization benchmark completed"
    else
        log_warning "Rust benchmark encountered issues"
    fi

    echo ""

    if run_benchmark "perl_only" "$PERL_WORK_DIR" "false"; then
        log_success "Perl-only benchmark completed"
    else
        log_warning "Perl benchmark encountered issues"
    fi

    echo ""

    # Generate report
    compare_results

    echo ""
    log_info "Benchmark complete!"
    log_info "Results directory: $RESULTS_DIR"
    log_info "Logs directory: $LOG_DIR"

    # Print summary
    echo ""
    echo -e "${GREEN}========== BENCHMARK SUMMARY ==========${NC}"
    cat "$RESULTS_DIR/BENCHMARK_REPORT.md" || true
}

# Run main function
main "$@"
