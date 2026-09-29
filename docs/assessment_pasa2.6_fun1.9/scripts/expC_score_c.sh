#!/usr/bin/bash -l
#SBATCH -p exfab -A exfab -c 2 --mem 16G -t 4:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_C/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_C/logs/%x_%j.err
# Experiment C scoring: per genome, mask = gene spans of the training models of ALL arms and repeats
# (PASA final_training_models.gff3, BUSCO busco.final.gff3); score every repeat on the same gene set.
set -uo pipefail
module load gffcompare
C=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_C
OUT=$C/scores.tsv; rm -f $OUT
tail -n +2 $C/genomes.tsv | while IFS=$'\t' read -r NAME ASMID ACC SPECIES STRAIN LOCUSTAG GROUP STRATUM TD REF HAS BUSCODB; do
  MASK=( $(ls $C/$NAME/pasa_r*/out/predict_misc/final_training_models.gff3 $C/$NAME/busco_r*/out/predict_misc/busco.final.gff3 2>/dev/null) )
  echo "$NAME mask files: ${#MASK[@]}"
  for W in $C/$NAME/pasa_r* $C/$NAME/busco_r*; do
    [ -d "$W" ] || continue
    grep -q "^exit_status	0" $W/run_status.tsv 2>/dev/null || { echo "skip $W (failed or unfinished)"; continue; }
    PRED=$(ls $W/out/predict_results/*.gff3 | grep -v '\.tbl' | head -1)
    /usr/bin/python3.12 $C/masked_score.py --ref-gff3 $REF --pred-gff3 $PRED --mask "${MASK[@]}" \
        --genome $NAME --variant $(basename $W) --out $OUT > /dev/null
  done
done
echo "wrote $OUT"
