#!/bin/bash
#SBATCH --job-name=trinity_bench_perl
#SBATCH --account=stajichlab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=64GB
#SBATCH --time=1-00:00:00
#SBATCH --output=slurm_logs/%x_%A.out
#SBATCH --error=slurm_logs/%x_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jasonst@ucr.edu

# ============================================================================
# Trinity Perl-engine Benchmark - SLURM Job Script
# ============================================================================
#
# Runs funannotate train with the Rust Trinity utilities hidden from PATH
# (falls back to the Perl versions). Submit alongside
# benchmark_trinity_rust_rust.sbatch (separate job, own walltime budget) then
# run install_scripts/compare_trinity_rust_results.sh once both finish. Splitting
# these into two jobs is what fixes the previous single-job runs hitting the
# SLURM --time limit before the perl leg even got through Trinity assembly.
#
# Submission:
#   sbatch benchmark_trinity_rust_perl.sbatch
#
# Monitor:
#   squeue -u $USER
#   tail -f slurm_logs/trinity_bench_perl_<job_id>.out

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

PROJECT_ROOT="/bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live"
DATASET_DIR="$HOME/projects/funannotate/trinity_example/Cordyceps_militaris"
BENCHMARK_OUTPUT_DIR="/bigdata/stajichlab/jstajich/projects/funannotate/BENCHMARK/trinity_rust"
ENGINE="perl"

# Number of CPUs to hand to the benchmark. Prefer the CPUs SLURM actually
# allocated to this job (SLURM_CPUS_PER_TASK), falling back to nproc when run
# outside SLURM. This keeps the thread count in sync with --cpus-per-task above
# instead of relying on a hardcoded value in the benchmark script.
CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# Node-local scratch for the working directory (see tests/benchmark_*.sh);
# normally auto-set by the workspace/scratch module, but default it defensively.
export SCRATCH="${SCRATCH:-/scratch/${USER}/${SLURM_JOB_ID}}"

# Ensure log directory exists
mkdir -p slurm_logs

# ============================================================================
# ENVIRONMENT SETUP
# ============================================================================

echo "========== Trinity Perl-engine Benchmark Job =========="
echo "Job ID: $SLURM_JOB_ID"
echo "Host: $(hostname)"
echo "Date: $(date)"
echo ""

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

if ! command -v pixi &> /dev/null; then
    echo "ERROR: pixi not found in PATH"
    echo "Please ensure pixi is installed or loaded via modules"
    exit 1
fi

echo "✓ Project root: $PROJECT_ROOT"
echo "✓ Dataset: $DATASET_DIR"
echo "✓ Output: $BENCHMARK_OUTPUT_DIR"
echo "✓ Engine: $ENGINE"
echo "✓ CPUs: $CPUS"
echo ""

# ============================================================================
# PIXI ENVIRONMENT ACTIVATION
# ============================================================================

echo "Initializing pixi environment..."
cd "$PROJECT_ROOT"

eval "$(pixi shell-hook -s bash)"
pixi install

export PATH="${CONDA_PREFIX}/bin:${PATH}"

echo "✓ Pixi environment ready"
echo "  CONDA_PREFIX: $CONDA_PREFIX"
echo "  Trinity path: $(command -v Trinity 2>/dev/null || echo 'NOT FOUND')"
echo ""

# ============================================================================
# DEPENDENCY VALIDATION
# ============================================================================

echo "Checking dependencies..."

for tool in funannotate Trinity samtools hisat2; do
    if command -v "$tool" &> /dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ✗ $tool NOT FOUND"
        exit 1
    fi
done

echo ""
echo "Checking Rust Trinity tools..."
for tool in sam_to_read_coords extract_reads_per_partition fragment_coverage_writer define_coverage_partitions; do
    if command -v "$tool" &> /dev/null; then
        echo "  ✓ $tool"
    else
        echo "  ⚠ $tool not found"
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

echo "Starting benchmark (engine=$ENGINE)..."
echo "===================="
echo ""

BENCHMARK_START=$(date +%s)

if "$PROJECT_ROOT/tests/benchmark_trinity_rust.sh" \
    "$DATASET_DIR" \
    "$BENCHMARK_OUTPUT_DIR" \
    "$ENGINE" \
    "$CPUS"; then

    BENCHMARK_END=$(date +%s)
    BENCHMARK_DURATION=$((BENCHMARK_END - BENCHMARK_START))

    echo ""
    echo "========== BENCHMARK COMPLETED =========="
    echo "Duration: $((BENCHMARK_DURATION / 60)) minutes"
    echo ""
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
echo "Engine: $ENGINE"
echo "Results: $BENCHMARK_OUTPUT_DIR"
echo ""
echo "Once benchmark_trinity_rust_rust.sbatch also completes, run:"
echo "  install_scripts/compare_trinity_rust_results.sh $BENCHMARK_OUTPUT_DIR"
echo ""
