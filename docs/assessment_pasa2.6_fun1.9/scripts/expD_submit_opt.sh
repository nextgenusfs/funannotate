#!/usr/bin/bash
# Submit experiment D: for each genome in expD_genomes.tsv, arms pasa and busco, as VARIANT=opt
# (--optimize_augustus) and VARIANT=ctl (same image, no optimize); step A, then step B after it.
# Genomes flagged `repeat` get a second opt run (_r2). Jobs go to the stajichlab partition.
# Prints the sbatch commands only; set SUBMIT=1 to submit. Run from the experiment_B directory
# with arm_opt.sh and expD_genomes.tsv next to this script's copy there:
#   cd .../pasa_train_performance_evaluate/experiment_B && SUBMIT=1 bash submit_opt.sh
set -euo pipefail
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B
cd $E
mkdir -p logs
: > jobs_opt.tsv
run() { if [ "${SUBMIT:-0}" = 1 ]; then sbatch --parsable "$@"; else echo "sbatch $*" >&2; echo DRY; fi; }
n=0
tail -n +2 expD_genomes.tsv | while IFS=$'\t' read -r G GROUP SN REPEAT; do
  for VARIANT in opt ctl; do
    reps=("")
    [ "$REPEAT" = 1 ] && [ $VARIANT = opt ] && reps+=(2)
    for ARM in pasa busco; do
      for REP in "${reps[@]}"; do
        JA=$(run -p stajichlab --export=ALL,VARIANT=$VARIANT -J "D_${G}_${ARM}_${VARIANT}${REP:+_r$REP}_A" arm_opt.sh $G $ARM A $REP)
        JB=$(run -p stajichlab --export=ALL,VARIANT=$VARIANT -J "D_${G}_${ARM}_${VARIANT}${REP:+_r$REP}_B" --dependency=afterok:$JA arm_opt.sh $G $ARM B $REP)
        echo -e "$G\t$VARIANT\t$ARM\t$REP\t$JA\t$JB" >> jobs_opt.tsv
      done
    done
  done
done
echo "$(wc -l < jobs_opt.tsv) arm pairs ($(( $(wc -l < jobs_opt.tsv) * 2 )) jobs)"
