# App image: adds the funannotate Python package on top of the pre-built
# funannotate-base image (conda env + Rust-optimized Trinity/PASA/EVM --
# see Dockerfile.base). This is the fast-changing layer: ordinary code
# releases only rebuild this file, not the expensive base. BASE_TAG lets
# container.yml pin a specific known-good base at release time
# (--build-arg BASE_TAG=vX); it defaults to latest for local builds.
ARG BASE_TAG=latest
FROM ghcr.io/nextgenusfs/funannotate-base:${BASE_TAG} AS runtime

# FUNANNOTATE_REF is the git tag/branch/sha to pip-install. container.yml
# sets it to the pushed v* tag itself (e.g. v1.9.0) on a release build, so
# the image always contains exactly the tagged/released code, not a
# hardcoded branch. Default here is only a fallback for ad-hoc local
# `docker build .` runs with no --build-arg.
ARG FUNANNOTATE_REF=master

# funannotate-base ends as USER funannotate with PIXI_ENV_PATH owned by
# root (PIXI_ENV_PATH itself is inherited as an ENV var from the base
# image), so installing into it needs root back briefly.
USER root
RUN "${PIXI_ENV_PATH}/bin/pip" install --no-cache-dir "git+https://github.com/nextgenusfs/funannotate.git@${FUNANNOTATE_REF}"

# Verify installation as the runtime user. Base image already checked
# EVM/PASA/Trinity; this only needs to confirm funannotate itself resolves.
USER funannotate
WORKDIR /work
SHELL ["/bin/bash", "-c"]
RUN set -euo pipefail && funannotate --version

CMD ["funannotate", "--help"]
