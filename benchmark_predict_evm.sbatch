#!/bin/bash
#SBATCH --job-name=bench_evm
#SBATCH --array=0-5
#SBATCH --time=0-12:00:00
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_logs/%x_%A_%a.out
#SBATCH --account=stajichlab
#SBATCH --partition=intel
#SBATCH --constraint=broadwell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=96G

# Clean EVM benchmark: FUNANNOTATE_EVM_ENGINE=rust vs perl for `funannotate
# predict`, 3 replicates. Each task predicts on its own copy of the train output
# from benchmark_prep_train. The engine switch is the environment variable that
# predict.py reads, so it is reliable. Needs benchmark_prep_train to have finished.
#
#   sbatch --dependency=afterok:<prep_jobid> benchmark_predict_evm.sbatch

set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live/benchmark_env.sh
mkdir -p "$PROJECT_ROOT/slurm_logs"
bench_pick "rust perl"
bench_banner "EVM benchmark"

TRAIN_DIR="${TRAIN_DIR:-$BENCH_ROOT/prep/work_rust_rep1/output}"
OUT="$BENCH_ROOT/evm"
cd "$PROJECT_ROOT"
pixi run --frozen -- \
    bash tests/benchmark_predict_evm.sh "$TRAIN_DIR" "$OUT" "$ENGINE" "$PERL_EVM_HOME" "$BENCH_CPUS" "$REP"
