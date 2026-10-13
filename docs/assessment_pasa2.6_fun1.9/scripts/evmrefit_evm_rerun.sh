#!/usr/bin/bash
# Rerun EVM only, with a new weights file, on the saved inputs of one predict run.
# Usage: evm_rerun.sh RUN_OUT_DIR WEIGHTS_FILE RESULT_GFF3 [CPUS]
#   RUN_OUT_DIR: a predict -o folder (uses its predict_misc/gene_predictions.gff3, protein and
#   transcript alignments, genome.softmasked.fa). Same runEVM flags as predict (-m 10 -i 1500).
# Work goes to $SCRATCH; only the evm.round1-equivalent GFF3 is copied to RESULT_GFF3.
set -euo pipefail
OUT=$1; WTS=$2; RES=$3; CPUS=${4:-4}
SIF=/bigdata/stajichlab/shared/singularity_cache/funannotate-1.9.0-rc.3.sif
M=$OUT/predict_misc
T=$(mktemp -d ${SCRATCH:?}/evmre.XXXXXX)
TX=(); [ -s $M/transcript_alignments.gff3 ] && TX=(-t $M/transcript_alignments.gff3)
if [ ${#TX[@]} -eq 0 ]; then grep -v "^TRANSCRIPT" $WTS > $T/weights.evm.txt; else cp $WTS $T/weights.evm.txt; fi
apptainer exec --bind $M:$M,$T:$T $SIF python3 \
  /pixi/.pixi/envs/base/lib/python3.8/site-packages/funannotate/aux_scripts/funannotate-runEVM.py \
  -w $T/weights.evm.txt -c $CPUS -g $M/gene_predictions.gff3 -d $T/EVM -f $M/genome.softmasked.fa \
  -l $T/evm.log -m 10 -i 1500 -o $T/evm.round1.gff3 \
  -p $M/protein_alignments.gff3 "${TX[@]}" > $T/stdout.log 2>&1
mkdir -p $(dirname $RES); cp $T/evm.round1.gff3 $RES
rm -rf $T
