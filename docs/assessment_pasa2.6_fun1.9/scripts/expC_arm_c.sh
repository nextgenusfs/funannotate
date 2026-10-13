#!/usr/bin/bash -l
#SBATCH -p exfab -A exfab -c 16 --mem 64G -t 12:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_C/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_C/logs/%x_%j.err
# Experiment C (D118): train Augustus/SNAP on the WHOLE genome and predict it, 3 repeats per job.
# Usage: sbatch -J NAME arm_c.sh GENOME ARM
#   ARM pasa  : train from the whole-genome rc.3 PASA models (--min_pasa_complete_models 0)
#       busco : BUSCO-forced training (--min_pasa_complete_models 1e9), same evidence
# Inputs are the experiment B whole-genome files (genome.fa, pasa.genome.gff3, trinity.genome.bam,
# transcripts.genome.gff3, genemark.genome.gtf). Flags match experiment B (arm_b.sh) and the BFD
# predict call. Scoring is a separate step (score_c.sh) because the mask is the union of all arms.
set -euo pipefail
module load apptainer
GENOME=$1; ARM=$2
C=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_C
B=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.3.sif
FDB=/bigdata/stajichlab/shared/lib/funannotate_db
IFS=$'\t' read -r NAME ASMID ACC SPECIES STRAIN LOCUSTAG GROUP STRATUM TD REF HAS BUSCODB \
    < <(awk -F'\t' -v g="$GENOME" '$1==g' $C/genomes.tsv)
[ -n "${NAME:-}" ] || { echo "unknown genome $GENOME"; exit 2; }
G=$B/$NAME
case $ARM in
  pasa)  GATE=(--min_pasa_complete_models 0) ;;
  busco) GATE=(--min_pasa_complete_models 1000000000) ;;
  *) echo "unknown arm $ARM"; exit 2 ;;
esac
export APPTAINERENV_FUNANNOTATE_DB=$FDB
for REP in 1 2 3; do
  W=$C/$NAME/${ARM}_r$REP
  if [ -s $W/run_status.tsv ] && grep -q "^exit_status	0" $W/run_status.tsv; then echo "skip $W (done)"; continue; fi
  rm -rf $W; mkdir -p $W/augustus; cd $W
  cp -r $R/lib/augustus/3.5/config $W/augustus/config
  export TMPDIR=${SCRATCH:?}/expC_${NAME}_${ARM}_$REP; mkdir -p $TMPDIR
  export APPTAINERENV_AUGUSTUS_CONFIG_PATH=$W/augustus/config
  echo -e "genome\t$NAME\narm\t$ARM\nrep\t$REP\nsif\t$SIF\ngate\t${GATE[*]}" > $W/arm_info.tsv
  BINDS="$W:$W,$G:$G,$R/lib:$R/lib,$FDB:$FDB,$TMPDIR:$TMPDIR"
  start=$(date +%s)
  set +e
  apptainer exec --bind $BINDS $SIF funannotate predict --name $LOCUSTAG -i $G/genome.fa --strain "$STRAIN" \
      -o $W/out -s "$SPECIES" --cpu 16 --busco_db $BUSCODB \
      --AUGUSTUS_CONFIG_PATH $W/augustus/config -w codingquarry:0 glimmerhmm:0 genemark:1 \
      --min_training_models 30 --tmpdir $TMPDIR --SeqCenter NCBI --keep_no_stops --header_length 24 \
      --protein_evidence $R/lib/swissprot_fungi.faa --max_intronlen 3000 --min_intronlen 10 \
      --tbl2asn "-l paired-ends" --table 1 --auto-skip-genemark --genemark_gtf $G/genemark.genome.gtf \
      --rna_bam $G/trinity.genome.bam --pasa_gff $G/pasa.genome.gff3 --transcript_alignments $G/transcripts.genome.gff3 \
      "${GATE[@]}" > $W/predict.capture.log 2>&1
  status=$?
  set -e
  echo -e "exit_status\t$status\nwall_seconds\t$(( $(date +%s) - start ))" > $W/run_status.tsv
  rm -rf $TMPDIR
  [ $status -eq 0 ] || echo "predict failed ($status) for $W"
done
echo done
