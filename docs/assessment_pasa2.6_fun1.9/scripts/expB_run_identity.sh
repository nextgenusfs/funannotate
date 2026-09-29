#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 32G -t 4:00:00
#SBATCH -J expB_identity
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B/read_identity/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B/read_identity/logs/%x_%j.err
# Median read identity + map rate per experiment B genome (funannotate train gate code, 200k reads).
set -euo pipefail
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B
FL=/bigdata/stajichlab/jstajich/projects/funannotate/funannotate-live
export PATH=$FL/.pixi/envs/default/bin:$PATH
export PYTHONPATH=$FL
TMP=${SCRATCH:?}/identity; mkdir -p $TMP
python $E/read_identity/measure_identity.py --genomes $E/genomes.tsv --expdir $E \
    --out $E/read_identity/read_identity.tsv --n_reads 200000 --cpus 16 --tmpdir $TMP
