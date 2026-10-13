#!/usr/bin/bash -l
#SBATCH --exclude=c[01-30]
#SBATCH -p epyc -c 16 --mem 96G -t 24:00:00
#SBATCH -o /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/confirm/logs/%x_%j.out
#SBATCH -e /bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/confirm/logs/%x_%j.err
# EVM weight confirmation (D127): production-style train (shared species Trinity + normalized reads)
# on the whole masked genome, as the rc.3 pilot ran it (pilot used --pasa_db mysql; sqlite here).
# Usage: sbatch -J NAME confirm_train.sh GENOME
set -euo pipefail
module load apptainer
C=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs/pasa_train_performance_evaluate/evm_refit/confirm
R=/bigdata/stajichlab/shared/projects/BFD/Fungi_BFD_runs
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.3.sif
IFS=$'\t' read -r NAME ASMID ACC SPECIES STRAIN LOCUSTAG GROUP STRATUM TD REF HAS BUSCODB GDIR \
    < <(awk -F'\t' -v g="$1" '$1==g' $C/genomes.tsv)
TAG=$(echo "$SPECIES" | tr ' ' '_')
G=$GDIR; W=$G/train; mkdir -p $W; cd $G
export TMPDIR=${SCRATCH:?}
[ -s genome.fa ] || zcat $R/input_clean_genomes/$ASMID.masked.fasta.gz > genome.fa
cp $R/rnaseq_data/$TAG.trinity-GG.fasta $W/$TAG.trinity-GG.fasta
BINDS="$G:$G,$R/rnaseq_reads:$R/rnaseq_reads,$TMPDIR:$TMPDIR"
cd $W
apptainer exec --bind $BINDS $SIF funannotate train -i $G/genome.fa -o $W/out \
    --trinity $W/$TAG.trinity-GG.fasta \
    --left_norm $R/rnaseq_reads/${TAG}_norm_R1.fastq.gz --right_norm $R/rnaseq_reads/${TAG}_norm_R2.fastq.gz \
    --species "$SPECIES" --strain "$STRAIN" --cpus 16 --memory 96G --header_length 24 \
    --jaccard_clip --no-progress --max_intronlen 3000 --pasa_min_avg_per_id 85 \
    --pasa_min_pct_aligned 70 --pasa_num_bp_splice 1 --pasa_db sqlite > $W/train.capture.log 2>&1
rm -f $W/$TAG.trinity-GG.fasta
echo done
