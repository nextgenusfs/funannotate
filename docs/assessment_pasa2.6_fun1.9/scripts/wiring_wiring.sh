#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]
#SBATCH -p epyc -c 16 --mem 96G -t 24:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/gate_wiring_test/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/gate_wiring_test/logs/%x_%j.err
# Gate wiring test (D118): rc.3 vs rc.3 + funannotate d18e67c (same code as a9b814f).
# Usage:
#   wiring.sh train   IMG GENOME                 -> $BD/GENOME/train_IMG
#   wiring.sh predict IMG GENOME RUN [ARGS...]   -> $BD/GENOME/RUN, trained from train_new,
#            in the production layout: training/ pruned with the BFD find rule, logfiles/ copied.
#   wiring.sh decide  IMG GENOME RUN [ARGS...]   -> as predict, but stops once the training
#            gate decisions are logged (pasa_gate or training_mode_switch row), to test triggers.
#   IMG: old = funannotate-1.9.0-rc.3.sif, new = funannotate-1.9.0-rc.3+d18e67c.sif
set -uo pipefail
module load apptainer
MODE=$1; IMG=$2; GENOME=$3; shift 3
BD=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/gate_wiring_test
B=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
READS=$R/do_annotation_rc3_pilot/rnaseq_reads
FDB=/bigdata/stajichlab/shared/lib/funannotate_db
case $IMG in
  old) SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.3.sif ;;
  new) SIF=$BD/image/funannotate-1.9.0-rc.3+d18e67c.sif ;;
  *) echo "unknown image $IMG"; exit 2 ;;
esac
IFS=$'\t' read -r NAME ASMID ACC SPECIES STRAIN LOCUSTAG GROUP STRATUM TD REF HAS BUSCODB \
    < <(awk -F'\t' -v g="$GENOME" '$1==g' $BD/genomes.tsv)
[ -n "${NAME:-}" ] || { echo "unknown genome $GENOME"; exit 2; }
G=$B/$NAME; TAG=$(echo "$SPECIES" | tr ' ' '_')
export TMPDIR=${SCRATCH:?}/wiring_${NAME}_${MODE}_$$; mkdir -p $TMPDIR
export APPTAINERENV_FUNANNOTATE_DB=$FDB
BINDS="$BD:$BD,$G:$G,$R/lib:$R/lib,$READS:$READS,$FDB:$FDB,$TMPDIR:$TMPDIR"

if [ $MODE = train ]; then
  W=$BD/$NAME/train_$IMG; rm -rf $W; mkdir -p $W; cd $W
  echo -e "mode\ttrain\nimage\t$SIF\ngenome\t$NAME\nreads\t$READS/${TAG}_norm_R1.fastq.gz" > $W/run_info.tsv
  start=$(date +%s)
  apptainer exec --bind $BINDS $SIF funannotate train -i $G/genome.fa -o $W/out \
      --left_norm $READS/${TAG}_norm_R1.fastq.gz --right_norm $READS/${TAG}_norm_R2.fastq.gz --aligners minimap2 \
      --species "$SPECIES" --strain "$STRAIN" --cpus 16 --memory 96G --header_length 24 \
      --jaccard_clip --no-progress --min_coverage 4 --max_intronlen 3000 --min_rnaseq_map_rate 10 \
      --pasa_db sqlite > $W/train.capture.log 2>&1
  status=$?
  echo -e "exit_status\t$status\nwall_seconds\t$(( $(date +%s) - start ))" > $W/run_status.tsv
  rm -rf $TMPDIR; exit 0
fi

RUN=$1; shift
W=$BD/$NAME/$RUN; rm -rf $W; mkdir -p $W/augustus
T=$BD/$NAME/train_new/out
grep -q "^exit_status	0" $BD/$NAME/train_new/run_status.tsv || { echo "train_new missing or failed for $NAME"; exit 3; }
# production layout: pruned training/ (FUNANNOTATE_TRAIN keep list) + copied logfiles/
cp -a $T/training $W/pred_training_full; mkdir -p $W/out; mv $W/pred_training_full $W/out/training
find "$W/out/training" -mindepth 1 -maxdepth 1 ! -name '.*' \
    ! -name funannotate_train.pasa.gff3 ! -name funannotate_train.coordSorted.bam \
    ! -name funannotate_train.transcripts.gff3 ! -name funannotate_train.trinity-GG.fasta \
    ! -name kallisto.tsv ! -name transcript.alignments.bam ! -name transcript.alignments.gff3 \
    ! -name trinity.alignments.bam ! -name trinity.alignments.gff3 ! -name trinity.fasta \
    -exec rm -rf {} +
cp -a $T/logfiles $W/out/logfiles
ls -la $W/out/training > $W/training_layout.txt
cp -r $R/lib/augustus/3.5/config $W/augustus/config
export APPTAINERENV_AUGUSTUS_CONFIG_PATH=$W/augustus/config
BINDS="$BINDS,$W:$W"
echo -e "mode\t$MODE\nimage\t$SIF\ngenome\t$NAME\nextra\t$*" > $W/run_info.tsv
CMD=(apptainer exec --bind $BINDS $SIF funannotate predict --name $LOCUSTAG -i $G/genome.fa --strain "$STRAIN"
    -o $W/out -s "$SPECIES" --cpu 16 --busco_db $BUSCODB
    --AUGUSTUS_CONFIG_PATH $W/augustus/config -w codingquarry:0 glimmerhmm:0 genemark:1
    --min_training_models 30 --tmpdir $TMPDIR --SeqCenter NCBI --keep_no_stops --header_length 24
    --protein_evidence $R/lib/swissprot_fungi.faa --max_intronlen 3000 --min_intronlen 10
    --tbl2asn "-l paired-ends" --table 1 --auto-skip-genemark --genemark_gtf $G/genemark.genome.gtf "$@")
start=$(date +%s)
if [ $MODE = predict ]; then
  "${CMD[@]}" > $W/predict.capture.log 2>&1; status=$?
else
  setsid "${CMD[@]}" > $W/predict.capture.log 2>&1 & pid=$!
  D=$W/out/logfiles/training_decisions.tsv; status=running
  while kill -0 $pid 2>/dev/null; do
    if [ -f $D ] && awk -F'\t' '$1=="predict" && ($2=="pasa_gate" || $2=="training_mode_switch")' $D | grep -q .; then
      sleep 20; kill -TERM -- -$pid 2>/dev/null; sleep 10; kill -KILL -- -$pid 2>/dev/null; status=stopped_after_decision; break
    fi
    sleep 15
  done
  [ $status = running ] && { wait $pid; status="exited_$?"; }
fi
echo -e "exit_status\t$status\nwall_seconds\t$(( $(date +%s) - start ))" > $W/run_status.tsv
rm -rf $TMPDIR
