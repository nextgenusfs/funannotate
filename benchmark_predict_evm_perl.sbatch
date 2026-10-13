#!/bin/bash
#SBATCH --job-name=predict_evm_bench_perl
#SBATCH --account=stajichlab
#SBATCH --partition=stajichlab
#SBATCH --nodes=1
#SBATCH --cpus-per-task=48
#SBATCH --mem=96GB
#SBATCH --time=0-08:00:00
#SBATCH --output=slurm_logs/%x_%A.out
#SBATCH --error=slurm_logs/%x_%A.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jasonst@ucr.edu

# ============================================================================
# Funannotate Predict Perl-EVM-engine Benchmark - SLURM Job Script
# ============================================================================
#
# Runs `funannotate predict` with the Rust EvidenceModeler engine, reusing a
# previously completed `funannotate train` output directory (see
# benchmark_trinity_rust_rust.sbatch / tests/benchmark_trinity_rust.sh).
# Submit alongside benchmark_predict_evm_perl.sbatch (separate job, own
# walltime budget) then run install_scripts/compare_predict_evm_results.sh once both
# finish.
#
# Submission:
#   sbatch benchmark_predict_evm_perl.sbatch
#
# Monitor:
#   squeue -u $USER
#   tail -f slurm_logs/predict_evm_bench_perl_<job_id>.out

set -euo pipefail

# ============================================================================
# CONFIGURATION
# ============================================================================

PROJECT_ROOT="/bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live"
# Output directory from a completed `funannotate train` run (must contain training/)
TRAIN_DIR="${TRAIN_DIR:-/bigdata/stajichlab/jstajich/projects/funannotate/BENCHMARK/trinity_rust/work_perl/output}"
BENCHMARK_OUTPUT_DIR="/bigdata/stajichlab/jstajich/projects/funannotate/BENCHMARK/predict_evm"
PERL_EVM_HOME="${PERL_EVM_HOME:-/bigdata/stajichlab/jstajich/projects/funannotate/EVidenceModeler}"
ENGINE="perl"

# Number of CPUs to hand to the benchmark. Prefer the CPUs SLURM actually
# allocated to this job (SLURM_CPUS_PER_TASK), falling back to nproc when run
# outside SLURM. This keeps the thread count in sync with --cpus-per-task above
# instead of relying on a hardcoded value in the benchmark script.
CPUS="${SLURM_CPUS_PER_TASK:-$(nproc)}"

# Node-local scratch for the working directory (see tests/benchmark_*.sh);
# normally auto-set by the workspace/scratch module, but default it defensively.
export SCRATCH="${SCRATCH:-/scratch/${USER}/${SLURM_JOB_ID}}"
# Shared lab funannotate database (see `module show funannotate/dev_live`)
export FUNANNOTATE_DB="${FUNANNOTATE_DB:-/srv/projects/db/funannotate/1.8}"

mkdir -p slurm_logs

echo "========== Predict Perl-EVM Benchmark Job =========="
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

if [ ! -d "$TRAIN_DIR" ]; then
    echo "ERROR: Train output directory not found: $TRAIN_DIR"
    echo "Run benchmark_trinity_rust_perl.sbatch first (or set \$TRAIN_DIR) to produce one."
    exit 1
fi

if ! command -v pixi &> /dev/null; then
    echo "ERROR: pixi not found in PATH"
    exit 1
fi

echo "✓ Project root: $PROJECT_ROOT"
echo "✓ Train dir: $TRAIN_DIR"
echo "✓ Perl EVM home: $PERL_EVM_HOME"
echo "✓ Engine: $ENGINE"
echo "✓ Output: $BENCHMARK_OUTPUT_DIR"
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
echo "  evidence_modeler path: $(command -v evidence_modeler 2>/dev/null || echo 'NOT FOUND')"
echo ""

if [ -z "${FUNANNOTATE_DB:-}" ]; then
    echo "ERROR: \$FUNANNOTATE_DB is not set. export it before submitting (funannotate setup required)."
    exit 1
fi

# ============================================================================
# RUN BENCHMARK
# ============================================================================

echo "Starting benchmark (engine=$ENGINE)..."
echo "===================="
echo ""

BENCHMARK_START=$(date +%s)

if "$PROJECT_ROOT/tests/benchmark_predict_evm.sh" \
    "$TRAIN_DIR" \
    "$BENCHMARK_OUTPUT_DIR" \
    "$ENGINE" \
    "$PERL_EVM_HOME" \
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

echo ""
echo "========== JOB SUMMARY =========="
echo "Job ID: $SLURM_JOB_ID"
echo "Status: COMPLETED"
echo "Engine: $ENGINE"
echo "Results: $BENCHMARK_OUTPUT_DIR"
echo ""
echo "Once benchmark_predict_evm_rust.sbatch also completes, run:"
echo "  install_scripts/compare_predict_evm_results.sh $BENCHMARK_OUTPUT_DIR"
echo ""
