#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 32G -t 4:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/logs/%x_%j.err
# EVM refit: rerun EVM for every weight set on one experiment B run, then score the holdout.
# Usage: sbatch -J NAME refit_job.sh GENOME ARM   (ARM = pasa.B or busco_r1.B)
# 4 EVM runs in parallel, 4 CPUs each. Weight sets: weight_sets.tsv / weights/*.txt.
set -euo pipefail
module load apptainer
GENOME=$1; ARM=$2
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit
B=$E/../experiment_B
OUT=$B/$GENOME/$ARM/out
REF=$(awk -F'\t' -v g=$GENOME '$1==g{print $10}' $B/genomes.tsv)
R=$E/runs/$GENOME/$ARM; mkdir -p $R
ls $E/weights/*.txt | xargs -P 4 -I{} bash -c 'n=$(basename {} .txt); [ -s '$R'/$n.gff3 ] || '$E'/evm_rerun.sh '$OUT' {} '$R'/$n.gff3 4'
rm -f $R/scores.tsv
/usr/bin/python3.12 $E/evm_score.py --ref-gff3 $REF --holdout $B/$GENOME/split.holdout_chroms.txt \
    --genome $GENOME --out $R/scores.tsv $(for f in $R/*.gff3; do echo "$(basename $f .gff3)=$f"; done)
echo done
