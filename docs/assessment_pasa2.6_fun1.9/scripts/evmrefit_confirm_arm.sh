#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]
#SBATCH -p epyc -c 16 --mem 64G -t 16:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/confirm/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/confirm/logs/%x_%j.err
# EVM weight confirmation (D127): full funannotate predict (post-EVM filters included), experiment B
# design, PASA training arm only.
# Usage: sbatch -J NAME confirm_arm.sh GENOME A            -> GENOME_DIR/pasa.A (train chromosomes)
#        sbatch -J NAME confirm_arm.sh GENOME B WTS REP    -> confirm/GENOME/WTS_rREP.B (whole genome,
#                                                           -p GENOME_DIR/pasa.A parameters)
#   WTS base = current BFD weights (-w codingquarry:0 glimmerhmm:0 genemark:1)
#       new  = refit weights (D126): augustus:1 hiq:3 genemark:2 snap:1 pasa:4
set -euo pipefail
module load apptainer
GENOME=$1; STEP=$2; WTS=${3:-base}; REP=${4:-1}
C=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/confirm
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
X=$R/pasa_train_performance_evaluate
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.3.sif
FDB=/bigdata/stajichlab/shared/lib/funannotate_db
IFS=$'\t' read -r NAME ASMID ACC SPECIES STRAIN LOCUSTAG GROUP STRATUM TD REF HAS BUSCODB G \
    < <(awk -F'\t' -v g="$GENOME" '$1==g' $C/genomes.tsv)
[ -n "${NAME:-}" ] || { echo "unknown genome $GENOME"; exit 2; }
case $WTS in
  base) W_ARGS=(codingquarry:0 glimmerhmm:0 genemark:1) ;;
  new)  W_ARGS=(codingquarry:0 glimmerhmm:0 augustus:1 hiq:3 genemark:2 snap:1 pasa:4) ;;
  *) echo "unknown weights $WTS"; exit 2 ;;
esac
if [ $STEP = A ]; then
  W=$G/pasa.A
  GEN=$G/genome_train.fa; BAM=$G/trinity.genome_train.bam; TRX=$G/transcripts.genome_train.gff3
  GMK=$G/genemark.genome_train.gtf; PASA=$G/pasa.genome_train.gff3
  EXTRA=(--min_pasa_complete_models 0)
else
  W=$C/$NAME/${WTS}_r$REP.B
  GEN=$G/genome.fa; BAM=$G/trinity.genome.bam; TRX=$G/transcripts.genome.gff3
  GMK=$G/genemark.genome.gtf; PASA=$G/pasa.genome.gff3
  EXTRA=(-p "$(ls $G/pasa.A/out/predict_results/*.parameters.json)")
fi
rm -rf $W; mkdir -p $W/augustus; cd $W
cp -r $R/lib/augustus/3.5/config $W/augustus/config
export TMPDIR=${SCRATCH:?}
echo -e "genome\t$NAME\nstep\t$STEP\nweights\t${W_ARGS[*]}\nrep\t$REP\nextra\t${EXTRA[*]}" > $W/arm_info.tsv
export APPTAINERENV_FUNANNOTATE_DB=$FDB
export APPTAINERENV_AUGUSTUS_CONFIG_PATH=$W/augustus/config
BINDS="$W:$W,$G:$G,$R/lib:$R/lib,$FDB:$FDB,$TMPDIR:$TMPDIR"
start=$(date +%s)
set +e
apptainer exec --bind $BINDS $SIF funannotate predict --name $LOCUSTAG -i $GEN --strain "$STRAIN" \
    -o $W/out -s "$SPECIES" --cpu 16 --busco_db $BUSCODB \
    --AUGUSTUS_CONFIG_PATH $W/augustus/config -w "${W_ARGS[@]}" \
    --min_training_models 30 --tmpdir $TMPDIR --SeqCenter NCBI --keep_no_stops --header_length 24 \
    --protein_evidence $R/lib/swissprot_fungi.faa --max_intronlen 3000 --min_intronlen 10 \
    --tbl2asn "-l paired-ends" --table 1 --auto-skip-genemark --genemark_gtf $GMK \
    --rna_bam $BAM --pasa_gff $PASA --transcript_alignments $TRX "${EXTRA[@]}" > $W/predict.capture.log 2>&1
status=$?
set -e
echo -e "exit_status\t$status\nwall_seconds\t$(( $(date +%s) - start ))" > $W/run_status.tsv
[ $status -eq 0 ] || { echo "predict failed ($status)"; exit $status; }
if [ $STEP = B ]; then
  PRED=$(ls $W/out/predict_results/*.gff3 | grep -v '\.tbl' | head -1)
  module load gffcompare
  /usr/bin/python3.12 $X/refseq_benchmark/predict_scorer.py score --ref-gff3 $REF --pred-gff3 $PRED \
      --holdout $G/split.holdout_chroms.txt --genome $NAME --variant ${WTS}_r$REP --out $W/score.tsv
  /usr/bin/python3.12 $X/evm_refit/evm_score.py --ref-gff3 $REF --holdout $G/split.holdout_chroms.txt \
      --genome $NAME --out $W/gene_score.tsv final=$PRED evm_round1=$W/out/predict_misc/evm.round1.gff3
  rm -rf $W/out/predict_misc/EVM $W/out/predict_misc/busco $W/out/predict_misc/tbl2asn
  grep "EVM Weights" $W/out/logfiles/funannotate-predict.log | tail -1 >> $W/arm_info.tsv
fi
echo done
