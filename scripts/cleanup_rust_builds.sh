#!/usr/bin/env bash
# Cleanup failed/incomplete Rust tool installations
#
# Use this script if you encounter errors like:
#   "fatal: destination path '...' already exists and is not an empty directory"
#
# This removes incomplete installations so they can be rebuilt from scratch.
#
# Usage:
#   ./scripts/cleanup_rust_builds.sh                 # Show what would be removed
#   ./scripts/cleanup_rust_builds.sh --force         # Actually remove

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "$SCRIPT_DIR")"
FORCE="${1:-}"

# Check if in pixi environment
if [ -z "${CONDA_PREFIX:-}" ]; then
    echo "ERROR: Not in a pixi/conda environment"
    echo "Please run: pixi shell"
    exit 1
fi

echo "Cleanup Script for Rust Tool Installations"
echo "=========================================="
echo ""
echo "CONDA_PREFIX: $CONDA_PREFIX"
echo ""

# Define cleanup targets
declare -a CLEANUP_DIRS=(
    "${CONDA_PREFIX}/opt/evm-rust/src"
    "${CONDA_PREFIX}/opt/pasa-rust-3.0/src"
)

declare -a CARGO_DIRS=(
    "${PROJECT_ROOT}/rust_bio_utils/target"
)

echo "Items to clean:"
echo ""

total_size=0
for dir in "${CLEANUP_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        size=$(du -sh "$dir" 2>/dev/null | cut -f1)
        echo "  ✓ $dir ($size)"
        total_size=$((total_size + $(du -s "$dir" 2>/dev/null | cut -f1)))
    fi
done

echo ""
echo "Build artifacts to clean:"
echo ""
for dir in "${CARGO_DIRS[@]}"; do
    if [ -d "$dir" ]; then
        size=$(du -sh "$dir" 2>/dev/null | cut -f1)
        echo "  ✓ $dir ($size)"
    fi
done

echo ""
echo "Total disk space to free: ~${total_size}KB"
echo ""

if [ "$FORCE" == "--force" ] || [ "$FORCE" == "-f" ]; then
    echo "Removing directories..."
    echo ""

    removed_count=0
    for dir in "${CLEANUP_DIRS[@]}"; do
        if [ -d "$dir" ]; then
            echo "Removing: $dir"
            rm -rf "$dir"
            ((removed_count++))
        fi
    done

    echo ""
    echo "Cleaning build artifacts..."
    for dir in "${CARGO_DIRS[@]}"; do
        if [ -d "$dir" ]; then
            echo "Removing: $dir"
            rm -rf "$dir"
            ((removed_count++))
        fi
    done

    echo ""
    echo "✓ Cleanup complete ($removed_count items removed)"
    echo ""
    echo "Next steps:"
    echo "  1. Rebuild: pixi install"
    echo "  2. Verify: which sam_to_read_coords"
    echo "  3. Retry:  sbatch benchmark_trinity_rust.sbatch"

else
    echo "To actually remove these files, run:"
    echo "  ./scripts/cleanup_rust_builds.sh --force"
    echo ""
    echo "This will free up disk space and allow fresh installation."
fi
