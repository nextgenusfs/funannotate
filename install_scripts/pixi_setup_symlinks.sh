#!/bin/bash
# bioconda's fasta3 package installs fasta36 only. Add a "fasta" link for tools
# that expect that name. funannotate and PASA also accept fasta36 and fasta3
# directly, so this link is a convenience. Link only when fasta36 exists, and
# replace a dangling link left by an earlier version of this script.
BIN=$CONDA_PREFIX/bin
if [ -x "$BIN/fasta36" ] && [ ! -e "$BIN/fasta" ]; then
    ln -sfn fasta36 "$BIN/fasta"
fi
