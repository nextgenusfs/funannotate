#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]
#SBATCH -p stajichlab -c 16 --mem 64G -t 2-00:00:00
#SBATCH -o logs/%x_%j.out -e logs/%x_%j.err
# One arm of experiment D: arm_b.sh with `predict --optimize_augustus` added in step A (VARIANT=opt), or
# the same run without it (VARIANT=ctl, an image-matched control). The rc.3 image used for experiment B is
# gone, so both use $SIF (default 1.9.0-rc.5); experiment B's plain runs (rc.3) are a second baseline.
# Everything else (image, evidence, flags, scoring) is identical to arm_b.sh, so the matching plain
# arm of experiment B (pasa.B, busco_r1.B) is a baseline. Results go to <arm>_<VARIANT>[_rN].{A,B}.
# optimize_augustus.pl runs for hours, hence the 2-day limit. Only arms pasa and busco are used.
# Original header of arm_b.sh follows.
# Usage: sbatch -J NAME arm_b.sh GENOME ARM STEP [REP]
#   ARM   pasa  : train from the rc.3 PASA models (--min_pasa_complete_models 0)
#         busco : BUSCO-forced training (--min_pasa_complete_models 1e9), same evidence
#         norna : BUSCO-forced, no RNA-seq evidence at all (sidecar; predict_misc kept)
#   STEP  A     : train on the train chromosomes (genome_train.fa)
#         B     : whole genome, -p parameters.json from step A; pasa/busco use the same
#                 fixed evidence (pilot PASA + transcripts + Trinity BAM); then score the
#                 holdout chromosomes against RefSeq and run proteome BUSCO
#   REP   repeat number for busco/norna arms (directory <arm>_r<REP>)
# Image: funannotate 1.9.0-rc.3 (PASA 2.6.1-rc.2). Flags match the BFD pipeline predict call.
set -euo pipefail
module load apptainer
GENOME=$1; ARM=$2; STEP=$3; REP=${4:-}
E=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/experiment_B
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
X=$R/pasa_train_performance_evaluate
SIF=${SIF:-/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.5.sif}
[ -s "$SIF" ] || { echo "image not found: $SIF"; exit 2; }
VARIANT=${VARIANT:-opt}
case $VARIANT in opt|ctl) ;; *) echo "VARIANT must be opt or ctl"; exit 2 ;; esac
FDB=/bigdata/stajichlab/shared/lib/funannotate_db
IFS=$'\t' read -r NAME ASMID ACC SPECIES STRAIN LOCUSTAG GROUP STRATUM TD REF HAS BUSCODB \
    < <(awk -F'\t' -v g="$GENOME" '$1==g' $E/genomes.tsv)
[ -n "${NAME:-}" ] || { echo "unknown genome $GENOME"; exit 2; }
G=$E/$NAME
AD=${ARM}_${VARIANT}${REP:+_r$REP}
W=$G/$AD.$STEP
mkdir -p $W; cd $W
export TMPDIR=${SCRATCH:?}
rm -rf $W/augustus $W/out
mkdir -p $W/augustus; cp -r $R/lib/augustus/3.5/config $W/augustus/config

EXTRA=()
if [ $STEP = A ]; then
  GEN=$G/genome_train.fa; BAM=$G/trinity.genome_train.bam; TRX=$G/transcripts.genome_train.gff3
  GMK=$G/genemark.genome_train.gtf; PASA=$G/pasa.genome_train.gff3
  case $ARM in
    pasa)  EXTRA+=(--min_pasa_complete_models 0) ;;
    busco|norna) EXTRA+=(--min_pasa_complete_models 1000000000) ;;
    *) echo "unknown arm $ARM"; exit 2 ;;
  esac
  [ $VARIANT = opt ] && EXTRA+=(--optimize_augustus)
else
  GEN=$G/genome.fa; BAM=$G/trinity.genome.bam; TRX=$G/transcripts.genome.gff3
  GMK=$G/genemark.genome.gtf; PASA=$G/pasa.genome.gff3
  EXTRA+=(-p "$(ls $G/$AD.A/out/predict_results/*.parameters.json)")
fi
if [ $ARM = norna ]; then
  EVID=()
else
  EVID=(--rna_bam $BAM --pasa_gff $PASA --transcript_alignments $TRX)
fi
echo -e "genome\t$NAME\narm\t$ARM\nrep\t${REP:-}\nstep\t$STEP\ngenome_fa\t$GEN\nevidence\t${EVID[*]:-none}\nextra\t${EXTRA[*]:-}\nvariant\t$VARIANT\nsif\t$SIF" > $W/arm_info.tsv

export APPTAINERENV_FUNANNOTATE_DB=$FDB
export APPTAINERENV_AUGUSTUS_CONFIG_PATH=$W/augustus/config
BINDS="$W:$W,$G:$G,$R/lib:$R/lib,$FDB:$FDB,$TMPDIR:$TMPDIR"
start=$(date +%s)
set +e
apptainer exec --bind $BINDS $SIF funannotate predict --name $LOCUSTAG -i $GEN --strain "$STRAIN" \
    -o $W/out -s "$SPECIES" --cpu 16 --busco_db $BUSCODB \
    --AUGUSTUS_CONFIG_PATH $W/augustus/config -w codingquarry:0 glimmerhmm:0 genemark:1 \
    --min_training_models 30 --tmpdir $TMPDIR --SeqCenter NCBI --keep_no_stops --header_length 24 \
    --protein_evidence $R/lib/swissprot_fungi.faa --max_intronlen 3000 --min_intronlen 10 \
    --tbl2asn "-l paired-ends" --table 1 --auto-skip-genemark --genemark_gtf $GMK \
    "${EVID[@]}" "${EXTRA[@]}" > $W/predict.capture.log 2>&1
status=$?
set -e
echo -e "exit_status\t$status\nwall_seconds\t$(( $(date +%s) - start ))" > $W/run_status.tsv
[ $status -eq 0 ] || { echo "predict failed ($status)"; exit $status; }

if [ $STEP = B ]; then
  PRED=$(ls $W/out/predict_results/*.gff3 | grep -v '\.tbl' | head -1)
  module load gffcompare
  /usr/bin/python3.12 $X/refseq_benchmark/predict_scorer.py score \
      --ref-gff3 $REF --pred-gff3 $PRED --holdout $G/split.holdout_chroms.txt \
      --genome $NAME --variant "$AD" --out $W/score.tsv
  module load busco/5.8.0
  busco -i $(ls $W/out/predict_results/*.proteins.fa) -m proteins -l /srv/projects/db/BUSCO/v10/lineages/fungi_odb10 \
      -o busco -c 16 --offline -f --out_path $TMPDIR > $W/busco.log 2>&1 || true
  cp $TMPDIR/busco/short_summary*.txt $W/busco_short_summary.txt 2>/dev/null || true
fi
echo done
