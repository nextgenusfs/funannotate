#!/bin/bash
#SBATCH --job-name=bench_trinity
#SBATCH --array=0-5
#SBATCH --time=0-14:00:00
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --account=stajichlab
#SBATCH --partition=intel
#SBATCH --constraint=broadwell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=96G

# Clean Trinity benchmark: Rust utilities vs TRINITY_NO_RUST=1, 3 replicates.
# Array index 0..5 -> rust,perl,rust,perl,rust,perl (replicate 1,1,2,2,3,3).
# Stage: funannotate train --stop_after_trinity (BENCH_STAGE=trinity).
#
#   sbatch benchmark_trinity_rust.sbatch
#
# The script tests/benchmark_trinity_rust.sh proves which engine ran and exits 4
# (INVALID) if the Rust-call count in the logs does not match the engine.

set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live/benchmark_env.sh
mkdir -p "$PROJECT_ROOT/slurm_logs"
bench_pick "rust perl"
bench_banner "Trinity benchmark"

OUT="$BENCH_ROOT/trinity"
cd "$PROJECT_ROOT"
# --frozen: use the environment as it is. Never `pixi install` in a benchmark job.
BENCH_STAGE=trinity pixi run --frozen -- \
    bash tests/benchmark_trinity_rust.sh "$DATASET_DIR" "$OUT" "$ENGINE" "$BENCH_CPUS" "$REP"
