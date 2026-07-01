#!/usr/bin/env bash
set -euo pipefail

# Install funannotate using setuptools
$PYTHON setup.py install --single-version-externally-managed --record=record.txt

# Verify installation
python -c "import funannotate; print(f'funannotate {funannotate.__version__} installed successfully')"

echo "✓ funannotate conda package built successfully"
echo "  Python: $PYTHON"
echo "  Prefix: $PREFIX"
