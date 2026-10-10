#!/bin/bash
#SBATCH --job-name=bench_pasa
#SBATCH --array=0-14
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
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jasonst@ucr.edu

# Clean PASA benchmark, 3 replicates of each engine:
#   rust       as shipped (slclust_rust, cdbyank_rust, faidx_rust; C++ assembler)
#   perl       all Rust PASA tools hidden
#   rustasm    rust plus the Rust assembler (PASA_ASSEMBLER=pasa_rust)
#   nocdb      rust with cdbyank_rust and faidx_rust hidden
#   noslclust  rust with slclust_rust hidden
# rust vs nocdb isolates cdbyank_rust; rust vs noslclust isolates slclust_rust.
# Array index 0..14 -> the five engines in the order above, per replicate. Needs benchmark_prep_train
# to have finished (its Trinity transcripts are the input).
#
#   sbatch --dependency=afterok:<prep_jobid> benchmark_pasa_rust.sbatch

set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live/benchmark_env.sh
mkdir -p "$PROJECT_ROOT/slurm_logs"
bench_pick "rust perl rustasm nocdb noslclust"
bench_banner "PASA benchmark"

export TRINITY_FASTA="${TRINITY_FASTA:-$BENCH_ROOT/prep/work_rust_rep1/output/training/funannotate_train.trinity-GG.fasta}"
OUT="$BENCH_ROOT/pasa"
cd "$PROJECT_ROOT"
pixi run --frozen -- \
    bash tests/benchmark_pasa_rust.sh "$DATASET_DIR" "$OUT" "$ENGINE" "$BENCH_CPUS" "$REP"
