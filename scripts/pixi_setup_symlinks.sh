#!/bin/bash
BIN=$CONDA_PREFIX/bin
if [ ! -L "$BIN/fasta" ]; then
    ln -s "$BIN/fasta36" "$BIN/fasta"
fi

# funannotate expects PASA binaries (seqclean, etc.) under $PASAHOME/bin, but
# PASA_rust's install.sh installs them into pasa-rust-3.0/bin, a sibling of
# pasa-rust-3.0/src (which is what PASAHOME points at). Symlink so both layouts work.
PASA_SRC_DIR="${CONDA_PREFIX}/opt/pasa-rust-3.0/src"
if [ -d "$PASA_SRC_DIR" ] && [ ! -e "$PASA_SRC_DIR/bin" ]; then
    ln -s ../bin "$PASA_SRC_DIR/bin"
fi
