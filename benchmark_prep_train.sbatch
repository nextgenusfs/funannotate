#!/bin/bash
#SBATCH --job-name=bench_prep_train
#SBATCH --time=1-00:00:00
#SBATCH --cpus-per-task=32
#SBATCH --output=slurm_logs/%x_%j.out
#SBATCH --account=stajichlab
#SBATCH --partition=intel
#SBATCH --constraint=broadwell
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --exclusive
#SBATCH --mem=0
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=jasonst@ucr.edu

# One full `funannotate train` with the shipped (Rust) configuration. It makes
# the inputs the PASA and EVM benchmarks reuse: the Trinity-GG transcripts and a
# complete train output directory. This run is NOT a timing leg.
#
#   sbatch benchmark_prep_train.sbatch

set -euo pipefail
source /bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live/benchmark_env.sh
mkdir -p "$PROJECT_ROOT/slurm_logs"
bench_banner "Prep: full train"

OUT="$BENCH_ROOT/prep"
mkdir -p "$OUT"
cd "$PROJECT_ROOT"
BENCH_STAGE=full pixi run --frozen -- \
    bash tests/benchmark_trinity_rust.sh "$DATASET_DIR" "$OUT" rust "$BENCH_CPUS" 1

TRAIN="$OUT/work_rust_rep1/output"
test -s "$TRAIN/training/funannotate_train.trinity-GG.fasta"
test -s "$TRAIN/training/funannotate_train.pasa.gff3"
echo "prep outputs OK: $TRAIN"
