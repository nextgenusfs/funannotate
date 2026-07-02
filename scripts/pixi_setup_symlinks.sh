#!/bin/bash
BIN=$CONDA_PREFIX/bin
if [ ! -L "$BIN/fasta" ]; then
    ln -s "$BIN/fasta36" "$BIN/fasta"
fi
