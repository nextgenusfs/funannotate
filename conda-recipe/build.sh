#!/usr/bin/env bash
set -euo pipefail

# Install funannotate using setuptools
$PYTHON setup.py install --single-version-externally-managed --record=record.txt

# Verify installation
python -c "import funannotate; print(f'funannotate {funannotate.__version__} installed successfully')"

# Build the source-built Rust components (trinity/PASA/EVM rust_optimize
# forks) using the same install scripts pixi and Docker use, so all three
# packaging paths stay in sync. These are NOT on bioconda/conda-forge --
# there is no package to depend on -- so they must be built from source
# here. NOTE: this clones from GitHub (pinned to specific commits inside
# each script), which means this recipe is only buildable locally
# (`conda build .` / `mamba build .`) with network access; it is NOT
# submittable to bioconda as-is, since bioconda's CI sandboxes build.sh
# without network access and expects source fetches to happen via the
# recipe's own `source:` section. Actual bioconda submission would require
# splitting pasa-rust/evidencemodeler-rust/trinity-rust into standalone
# recipes (each with its own `source: git_url/git_rev`) published first,
# matching how PASApipeline's own upstream conda-recipe is structured.
export CONDA_PREFIX="$PREFIX"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(dirname "${SCRIPT_DIR}")"
source "${REPO_ROOT}/scripts/pixi_install_trinity_rust.sh"
source "${REPO_ROOT}/scripts/pixi_install_evm_rust.sh"
source "${REPO_ROOT}/scripts/pixi_install_pasa_rust.sh"

echo "✓ funannotate conda package built successfully"
echo "  Python: $PYTHON"
echo "  Prefix: $PREFIX"
