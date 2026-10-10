# Shared settings for the clean Rust-vs-baseline benchmark (sourced by the
# benchmark_*.sbatch files by absolute path; no BASH_SOURCE, which is empty or
# wrong in SLURM work directories).
#
# Hardware: all jobs run on the `intel` partition with --constraint=broadwell,
# --exclusive, so every leg gets an identical whole node. Each job records its
# host and CPU model in results/<engine>_rep<N>.env, and
# tests/benchmark_summarize.py refuses to compare legs from different CPUs.

export PROJECT_ROOT="${PROJECT_ROOT:-/bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live}"
export DATASET_DIR="${DATASET_DIR:-/bigdata/stajichlab/jstajich/projects/funannotate/trinity_example/Cordyceps_militaris}"
export BENCH_ROOT="${BENCH_ROOT:-/bigdata/stajichlab/jstajich/projects/funannotate/BENCHMARK/clean_rust_2026-10}"
export PERL_EVM_HOME="${PERL_EVM_HOME:-/bigdata/stajichlab/jstajich/projects/funannotate/EVidenceModeler}"
export FUNANNOTATE_DB="${FUNANNOTATE_DB:-/srv/projects/db/funannotate/1.8}"
# Threads handed to funannotate. Fixed, so every leg uses the same count.
export BENCH_CPUS="${BENCH_CPUS:-32}"

# Node-local scratch; SLURM normally sets it.
export SCRATCH="${SCRATCH:-/scratch/${USER}/${SLURM_JOB_ID:-manual}}"
mkdir -p "$SCRATCH"

# Map a job-array index to (engine, replicate). Engines alternate within each
# replicate so that time-of-day and node effects are spread over all engines.
#   bench_pick "rust perl" 3   -> sets ENGINE and REP from $SLURM_ARRAY_TASK_ID
bench_pick() {
    local -a engines=($1)
    local idx="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not set}"
    local n="${#engines[@]}"
    ENGINE="${engines[$((idx % n))]}"
    REP=$((idx / n + 1))
    export ENGINE REP
}

bench_banner() {
    echo "=== $1 ==="
    echo "job:       ${SLURM_JOB_ID:-none} array task ${SLURM_ARRAY_TASK_ID:-none}"
    echo "host:      $(hostname)"
    echo "cpu:       $(lscpu | awk -F: '/Model name/ {gsub(/^ +/, "", $2); print $2; exit}')"
    echo "partition: ${SLURM_JOB_PARTITION:-none}"
    echo "date:      $(date -Is)"
    echo "repo:      $(git -C "$PROJECT_ROOT" rev-parse --short HEAD) $(git -C "$PROJECT_ROOT" status --porcelain --untracked-files=no | wc -l) modified tracked files"
    echo "engine:    ${ENGINE:-n/a} rep ${REP:-n/a}"
}
