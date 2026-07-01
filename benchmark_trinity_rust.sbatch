#!/bin/bash
#SBATCH --job-name=trinity_rust_bench
#SBATCH --account=stajichlab
#SBATCH --partition=epyc
#SBATCH --nodes=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32GB
#SBATCH --time=0-08:00:00
#SBATCH --output=slurm_logs/%x_%A.out
#SBATCH --error=slurm_logs/%x_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jasonst@ucr.edu

# ============================================================================
# Trinity Rust Optimization Benchmark - SLURM Job Script
# ============================================================================
#
# This script benchmarks funannotate Trinity workflow performance with and
# without Rust-optimized utilities on the UCR HPCC cluster.
#
# Submission:
#   sbatch benchmark_trinity_rust.sbatch
#
# Monitor:
#   squeue -u $USER
#   tail -f slurm_logs/trinity_rust_bench_<job_id>.out
#
# Results:
#   ~/bench_results/results/BENCHMARK_REPORT.md

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

PROJECT_ROOT="/bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live"
DATASET_DIR="$HOME/projects/funannotate/trinity_example/Cordyceps_militaris"
BENCHMARK_OUTPUT_DIR="$HOME/bench_results"

# Ensure log directory exists
mkdir -p slurm_logs

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

echo "========== Trinity Rust Benchmark Job =========="
echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo ""

# Module setup (if needed on your system)
# module load pixi 2>/dev/null || true

# ============================================================================
# PRE-RUN VALIDATION
# ============================================================================

echo "Validating setup..."

if [ ! -d "$PROJECT_ROOT" ]; then
    echo "ERROR: Project root not found: $PROJECT_ROOT"
    exit 1
fi

if [ ! -d "$DATASET_DIR" ]; then
    echo "ERROR: Dataset not found: $DATASET_DIR"
    exit 1
fi

# Check if pixi is available
if ! command -v pixi &> /dev/null; then
    echo "ERROR: pixi not found in PATH"
    echo "Please ensure pixi is installed or loaded via modules"
    exit 1
fi

echo "✓ Project root: $PROJECT_ROOT"
echo "✓ Dataset: $DATASET_DIR"
echo "✓ Output: $BENCHMARK_OUTPUT_DIR"
echo ""

# ============================================================================
# PIXI ENVIRONMENT ACTIVATION
# ============================================================================

echo "Initializing pixi environment..."
cd "$PROJECT_ROOT"

# Use pixi shell to activate environment
# This runs the activation scripts including pixi_install_trinity_rust.sh
eval "$(pixi shell-hook -s bash)"
pixi install

echo "✓ Pixi environment ready"
echo ""

# ============================================================================
# DEPENDENCY VALIDATION
# ============================================================================

echo "Checking dependencies..."

# Check essential tools
for tool in funannotate Trinity samtools hisat2; do
    if command -v "$tool" &> /dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool NOT FOUND"
        exit 1
    fi
done

# Check Rust Trinity tools
echo ""
echo "Checking Rust Trinity tools..."
for tool in sam_to_read_coords extract_reads_per_partition fragment_coverage_writer define_coverage_partitions; do
    if command -v "$tool" &> /dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ⚠ $tool not found (will use Perl versions)"
    fi
done

echo ""

# ============================================================================
# SYSTEM INFORMATION
# ============================================================================

echo "System Information"
echo "=================="
echo "CPUs available: $(nproc)"
echo "Memory total: $(free -h | grep Mem | awk '{print $2}')"
echo "Memory available: $(free -h | grep Mem | awk '{print $7}')"
echo "Disk space (project): $(df -h "$PROJECT_ROOT" | tail -1 | awk '{print $4}')"
echo ""

# ============================================================================
# RUN BENCHMARK
# ============================================================================

echo "Starting benchmark..."
echo "===================="
echo ""

BENCHMARK_START=$(date +%s)

# Run the main benchmark script
if "$PROJECT_ROOT/scripts/benchmark_trinity_rust.sh" \
    "$DATASET_DIR" \
    "$BENCHMARK_OUTPUT_DIR"; then

    BENCHMARK_END=$(date +%s)
    BENCHMARK_DURATION=$((BENCHMARK_END - BENCHMARK_START))

    echo ""
    echo "========== BENCHMARK COMPLETED =========="
    echo "Duration: $((BENCHMARK_DURATION / 60)) minutes"
    echo ""

    # Display results if available
    if [ -f "$BENCHMARK_OUTPUT_DIR/results/BENCHMARK_REPORT.md" ]; then
        echo "Results Summary:"
        echo "================"
        cat "$BENCHMARK_OUTPUT_DIR/results/BENCHMARK_REPORT.md"
    fi

else
    echo ""
    echo "ERROR: Benchmark failed!"
    echo "Check logs: $BENCHMARK_OUTPUT_DIR/logs/"
    exit 1
fi

# ============================================================================
# POST-RUN SUMMARY
# ============================================================================

echo ""
echo "========== JOB SUMMARY =========="
echo "Job ID: $SLURM_JOB_ID"
echo "Status: COMPLETED"
echo "Results: $BENCHMARK_OUTPUT_DIR"
echo "Report: $BENCHMARK_OUTPUT_DIR/results/BENCHMARK_REPORT.md"
echo "Logs: $BENCHMARK_OUTPUT_DIR/logs/"
echo ""
echo "To view results:"
echo "  cat $BENCHMARK_OUTPUT_DIR/results/BENCHMARK_REPORT.md"
echo "  ls -lh $BENCHMARK_OUTPUT_DIR/results/"
echo ""
