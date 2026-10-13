#!/usr/bin/bash -l
#SBATCH -p epyc -c 16 --mem 48G -t 6:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/logs/%x_%j.err
# EVM refit validation: rerun EVM with the finalist weight sets on held-out genomes.
# Usage: sbatch -J NAME validate_job.sh GENOME_LIST SET_LIST
#   GENOME_LIST: file of experiment B genome names; arms from $ARMS (default "pasa.B busco_r1.B").
#   SET_LIST:    file of weight-set names (weights/<name>.txt); must include base.
set -uo pipefail
module load apptainer
GL=$1; SL=$2
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit
B=$E/../experiment_B
while read -r GENOME; do
  REF=$(awk -F'\t' -v g=$GENOME '$1==g{print $10}' $B/genomes.tsv)
  for ARM in ${ARMS:-pasa.B busco_r1.B}; do
    OUT=$B/$GENOME/$ARM/out
    [ -s $OUT/predict_misc/gene_predictions.gff3 ] || { echo "skip $GENOME $ARM (no EVM inputs)"; continue; }
    R=$E/runs/$GENOME/$ARM; mkdir -p $R
    sed 's#.*#'$E'/weights/&.txt#' $SL | xargs -P 4 -I{} bash -c 'n=$(basename {} .txt); [ -s '$R'/$n.gff3 ] || '$E'/evm_rerun.sh '$OUT' {} '$R'/$n.gff3 4 || echo "EVM failed $n"'
    rm -f $R/scores.tsv
    /usr/bin/python3.12 $E/evm_score.py --ref-gff3 $REF --holdout $B/$GENOME/split.holdout_chroms.txt \
        --genome $GENOME --out $R/scores.tsv $(for n in $(cat $SL); do [ -s $R/$n.gff3 ] && echo "$n=$R/$n.gff3"; done)
    # compress EVM outputs once scored
    zstd -q -f --rm -T4 $R/*.gff3 2>/dev/null || true
  done
done < $GL
echo done
