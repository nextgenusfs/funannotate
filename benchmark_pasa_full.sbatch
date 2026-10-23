#!/bin/bash
#SBATCH --job-name=bench_pasa_full
#SBATCH --array=0-8
#SBATCH --time=0-04:00:00
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --account=stajichlab
#SBATCH --partition=intel
#SBATCH --constraint=broadwell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=96G

# PASA benchmark on the inputs of a real `funannotate train`: rust (as shipped), perl (all Rust PASA
# tools hidden) and rustasm (plus the Rust assembler), 3 replicates each. The first PASA benchmark ran
# train with --trinity, which skips hisat2 and StringTie and gives PASA a smaller job. Here a finished
# train output (benchmark_prep_train) is copied, PASA is made to run again, and train runs without
# --trinity. See the header of tests/benchmark_pasa_rust.sh.
# Array index 0..8 -> rust,perl,rustasm per replicate.
#
#   sbatch benchmark_pasa_full.sbatch

set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live/benchmark_env.sh
mkdir -p "$PROJECT_ROOT/slurm_logs"
bench_pick "rust perl rustasm"
bench_banner "PASA benchmark, full train inputs"

export PREP_OUTPUT="${PREP_OUTPUT:-$BENCH_ROOT/prep/work_rust_rep1/output}"
OUT="$BENCH_ROOT/pasa_full"
cd "$PROJECT_ROOT"
pixi run --frozen -- \
    bash tests/benchmark_pasa_rust.sh "$DATASET_DIR" "$OUT" "$ENGINE" "$BENCH_CPUS" "$REP"
