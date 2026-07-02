#!/usr/bin/env bash
# Build funannotate Docker images
# Usage: ./build.sh [DOCKERFILE] [TAG] [EXTRA_ARGS]
#
# Examples:
#   ./build.sh                                      # Build Dockerfile -> funannotate:latest
#   ./build.sh Dockerfile.dev                       # Build dev image
#   ./build.sh Dockerfile2 funannotate:with-db      # Build slim with databases

set -euo pipefail

DOCKERFILE="${1:-Dockerfile}"
TAG="${2:-funannotate:latest}"
EXTRA_ARGS="${3:---no-cache}"

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_ROOT="$(dirname "${SCRIPT_DIR}")"

# Validate Dockerfile exists
if [ ! -f "${PROJECT_ROOT}/${DOCKERFILE}" ]; then
    echo "Error: ${DOCKERFILE} not found in ${PROJECT_ROOT}" >&2
    exit 1
fi

echo "=== Building funannotate Docker image ==="
echo "Dockerfile:   ${DOCKERFILE}"
echo "Tag:          ${TAG}"
echo "Extra args:   ${EXTRA_ARGS}"
echo ""

# Check if we need to copy local checkouts for Dockerfile.dev. Dockerfile.dev
# COPYs trinityrnaseq/PASApipeline/EVidenceModeler (rust_optimize branch)
# directly from the build context root, so stage them there under those
# exact names (matching scripts/pixi_install_*.sh's local-checkout search).
if [ "${DOCKERFILE}" = "Dockerfile.dev" ]; then
    for repo in trinityrnaseq PASApipeline EVidenceModeler; do
        if [ ! -d "${PROJECT_ROOT}/../${repo}" ]; then
            echo "Error: Dockerfile.dev requires ../${repo} (rust_optimize branch) to be available" >&2
            exit 1
        fi
    done
    echo "Using local checkouts from ../trinityrnaseq, ../PASApipeline, ../EVidenceModeler"
    CLEANUP_DIRS=()
    for repo in trinityrnaseq PASApipeline EVidenceModeler; do
        if [ ! -e "${PROJECT_ROOT}/${repo}" ]; then
            cp -r "${PROJECT_ROOT}/../${repo}" "${PROJECT_ROOT}/${repo}"
            CLEANUP_DIRS+=("${PROJECT_ROOT}/${repo}")
        fi
    done
    BUILD_CONTEXT="${PROJECT_ROOT}"
else
    BUILD_CONTEXT="${PROJECT_ROOT}"
    CLEANUP_DIRS=()
fi

# Build the image
docker build ${EXTRA_ARGS} \
    -f "${PROJECT_ROOT}/${DOCKERFILE}" \
    -t "${TAG}" \
    "${BUILD_CONTEXT}"

# Cleanup temporary directories
for dir in "${CLEANUP_DIRS[@]}"; do
    rm -rf "$dir" 2>/dev/null || true
done

echo ""
echo "✓ Build complete!"
echo "  Image: ${TAG}"
echo ""
echo "Run with:"
echo "  docker run -it ${TAG} funannotate --help"
echo "  docker run -it -v \$(pwd)/data:/work ${TAG} funannotate predict -i genome.fasta ..."
