#!/bin/bash
#SBATCH --job-name=bench_train_full
#SBATCH --array=0-5
#SBATCH --time=0-06:00:00
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --account=stajichlab
#SBATCH --partition=intel
#SBATCH --constraint=broadwell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=96G

# The whole `funannotate train` step, Rust tools against the original tools, 3 replicates each
# (hisat2, genome-guided Trinity, StringTie, PASA, TransDecoder, kallisto, final model choice).
#   rust  everything as shipped
#   perl  TRINITY_NO_RUST=1 and all Rust PASA tools hidden
# Array index 0..5 -> rust,perl,rust,perl,rust,perl (replicate 1,1,2,2,3,3).
#
#   sbatch benchmark_train_full.sbatch

set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live/benchmark_env.sh
mkdir -p "$PROJECT_ROOT/slurm_logs"
bench_pick "rust perl"
bench_banner "Full train benchmark"

OUT="$BENCH_ROOT/train_full"
cd "$PROJECT_ROOT"
pixi run --frozen -- \
    bash tests/benchmark_train_full.sh "$DATASET_DIR" "$OUT" "$ENGINE" "$BENCH_CPUS" "$REP"
