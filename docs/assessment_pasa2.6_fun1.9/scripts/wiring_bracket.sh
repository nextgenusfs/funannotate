#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]
#SBATCH -p epyc -c 16 --mem 96G -t 12:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/gate_wiring_test/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/gate_wiring_test/logs/%x_%j.err
# Bracket runs for one genome: thresholds just below/above the measured values (decision-only).
# Identity: median read identity from train_new's report, ±0.1. PASA gate: the complete-model
# count recorded by pred_new's pasa_gate row, at the count (expect PASA) and count+1 (expect BUSCO).
set -uo pipefail
GENOME=$1
BD=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/gate_wiring_test
M=$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++)h[$i]=i} NR==2{print $h["median_identity_pct"]}' $BD/$GENOME/train_new/out/logfiles/train_rnaseq_gate.tsv)
N=$(awk -F'\t' '$1=="predict" && $2=="pasa_gate"{print $4}' $BD/$GENOME/pred_new/out/logfiles/training_decisions.tsv | head -1)
echo "median identity $M; complete models $N"
LO=$(awk -v m=$M 'BEGIN{printf "%.2f", m-0.1}'); HI=$(awk -v m=$M 'BEGIN{printf "%.2f", m+0.1}')
run() { local name=$1 expect=$2; shift 2
  bash $BD/wiring.sh decide new $GENOME $name "$@"
  printf "%s/%s\t%s\n" $GENOME $name $expect >> $BD/expect.tsv; }
run decide_identity_below pasa --min_rnaseq_identity $LO
run decide_identity_above busco --min_rnaseq_identity $HI
run decide_pasa_at pasa --min_pasa_complete_models $N
run decide_pasa_above busco --min_pasa_complete_models $((N+1))
