#!/usr/bin/env bash
#SBATCH -J build_sif
#SBATCH -p epyc
#SBATCH -c 8
#SBATCH --mem=32G
#SBATCH -t 02:00:00
# Build a funannotate beta image as a layer on top of a released image, without
# rebuilding the conda env, Trinity or PASA:
#   - EVidenceModeler is built from a git ref (tag) and its binaries and EvmUtils
#     replace the ones in the base image;
#   - funannotate is pip-installed from a git ref of this repository.
# The EVM binaries are built on the host. Host glibc (2.28) is older than the
# image's (2.41), so they run in the image.
#
# For the real release use CI (.github/workflows/container.yml); this is for
# testing a beta before a tag exists.
#
# Usage:
#   sbatch scripts/build_overlay_sif.sh [FUN_REPO] [FUN_REF] [EVM_REPO] [EVM_REF] [BASE_SIF] [OUT_SIF]
set -euo pipefail

FUN_REPO=${1:-$PWD}
FUN_REF=${2:-HEAD}
EVM_REPO=${3:-$HOME/projects/funannotate/EVidenceModeler}
EVM_REF=${4:-v2.5.1_rust-beta1}
CACHE=/bigdata/stajichlab/shared/singularity_cache
BASE_SIF=${5:-$CACHE/funannotate-1.9.0.sif}
OUT_SIF=${6:-$CACHE/funannotate-1.9.1-beta.1.sif}
W=${SCRATCH:?SCRATCH is not set; run inside a SLURM job}/sifbuild
APPTAINER_BIN=/opt/linux/rocky/8.x/x86_64/pkgs/apptainer/1.5.3/bin
export PATH="$APPTAINER_BIN:$PATH"

[ -f "$BASE_SIF" ] || { echo "missing base image $BASE_SIF" >&2; exit 1; }
FUN_COMMIT=$(git -C "$FUN_REPO" rev-parse "$FUN_REF^{commit}")
EVM_COMMIT=$(git -C "$EVM_REPO" rev-parse "$EVM_REF^{commit}")
echo "funannotate $FUN_REF = $FUN_COMMIT"
echo "EVidenceModeler $EVM_REF = $EVM_COMMIT"

rm -rf "$W"; mkdir -p "$W/stage/bin"
# Sources come from git archive, so uncommitted files never reach the image.
git -C "$EVM_REPO" archive "$EVM_COMMIT" --prefix=evm_src/ | tar -x -C "$W"
git -C "$FUN_REPO" archive "$FUN_COMMIT" --prefix=fun_src/ | tar -x -C "$W"

echo "== building EVM"
( cd "$W/evm_src/evm" && CARGO_TARGET_DIR="$W/evm_target" cargo build --release --locked )
for b in evidence_modeler EVidenceModeler convert_EVM_outputs_to_GFF3 partition_evm_inputs \
         recombine_evm_outputs gff3_file_to_proteins augustus_to_evm_gff3 \
         gff3_gene_prediction_file_validator; do
    if [ -x "$W/evm_target/release/$b" ]; then cp "$W/evm_target/release/$b" "$W/stage/bin/"
    else echo "note: $b not built" >&2; fi
done
cp -r "$W/evm_src/EvmUtils" "$W/stage/EvmUtils"
[ -d "$W/evm_src/PerlLib" ] && cp -r "$W/evm_src/PerlLib" "$W/stage/PerlLib"

cat > "$W/overlay.def" <<DEF
Bootstrap: localimage
From: $BASE_SIF

%files
    $W/stage /opt/beta_stage
    $W/fun_src /opt/fun_src

%post
    PIXI=\${PIXI_ENV_PATH:-/pixi/.pixi/envs/base}
    cp -f /opt/beta_stage/bin/* "\$PIXI/bin/"
    mkdir -p "\$PIXI/opt/evm"
    cp -rf /opt/beta_stage/EvmUtils "\$PIXI/opt/evm/"
    [ -d /opt/beta_stage/PerlLib ] && cp -rf /opt/beta_stage/PerlLib "\$PIXI/opt/evm/"
    "\$PIXI/bin/pip" install --no-deps --no-cache-dir --force-reinstall /opt/fun_src
    rm -rf /opt/beta_stage /opt/fun_src

%labels
    overlay.base $(basename "$BASE_SIF")
    overlay.funannotate_commit $FUN_COMMIT
    overlay.evm_ref $EVM_REF
    overlay.evm_commit $EVM_COMMIT

DEF

echo "== apptainer build"
apptainer build --fakeroot --force "$OUT_SIF.tmp" "$W/overlay.def"
mv "$OUT_SIF.tmp" "$OUT_SIF"
ls -la "$OUT_SIF"

# Checks run after the build, so a failed check names itself and leaves the image to inspect.
echo "== checks"
apptainer exec "$OUT_SIF" funannotate --version
fail=0
conv_help=$(apptainer exec "$OUT_SIF" convert_EVM_outputs_to_GFF3 --help 2>&1 || true)
for o in --skip_empty --concat_out --threads; do
    if grep -q -- "$o" <<< "$conv_help"; then echo "ok   convert_EVM_outputs_to_GFF3 has $o"
    else echo "FAIL convert_EVM_outputs_to_GFF3 lacks $o"; fail=1; fi
done
rec_help=$(apptainer exec "$OUT_SIF" recombine_evm_outputs --help 2>&1 || true)
if grep -q -- --threads <<< "$rec_help"; then echo "ok   recombine_evm_outputs has --threads"
else echo "FAIL recombine_evm_outputs lacks --threads"; fail=1; fi
pred_help=$(apptainer exec "$OUT_SIF" funannotate predict --help 2>&1 || true)
if grep -q -- --evm_tmpdir <<< "$pred_help"; then echo "ok   funannotate predict has --evm_tmpdir"
else echo "FAIL funannotate predict lacks --evm_tmpdir"; fail=1; fi
apptainer exec "$OUT_SIF" bash -c 'command -v convert_EVM_outputs_to_GFF3 evidence_modeler'
apptainer inspect --labels "$OUT_SIF"
exit $fail
