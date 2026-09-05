#!/usr/bin/env bash
# Build the three pipeline images.
#
# Why three: the stacks are mutually exclusive. econml/dowhy pin numpy 1.x and
# scikit-learn 1.4; h2o's environment here is numpy 2.2 and scikit-learn 1.7.
# One image cannot hold both, which is the same split that forced separate
# interpreters on the laptop (params.python / python_causal / python_h2o).
#
# Where: any host with a docker daemon and push access to $MTB_REGISTRY. An
# earlier arrangement used a host-local NAMING
# CONVENTION for images built on the node these apps run on, not a registry;
# convention with no registry behind it, which could not be pushed to. Nomad's
# docker driver finds the tag in the local image cache.
# (brainstorms/mtb-resistotyper-ml/2026-07-03-crypticdb-llm-abc-cluster-deploy-findings.md §2)
#
#   ./build.sh              build all three
#   ./build.sh causal       build one
#   TAG=v0.2.0 ./build.sh   override the tag

set -euo pipefail

# A real registry. The previous default was a local-daemon naming convention with no
# registry behind it, and Nomad's docker driver garbage-collects host-only images
# three minutes after the last task using them ends, leaving nothing to pull.
REGISTRY="${REGISTRY:-${MTB_REGISTRY:-ghcr.io/abhi18av-phd-projects/mtb-resistotyper-ml-pipeline}}"
TAG="${TAG:-v0.1.0}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/.." && pwd)"

IMAGES=("${@:-fe causal h2o}")
read -r -a IMAGES <<< "${IMAGES[*]}"

if ! docker info >/dev/null 2>&1; then
    echo "no docker daemon reachable." >&2
    echo "Start Docker, or build on a host with a daemon and access to $REGISTRY." >&2
    exit 1
fi

for name in "${IMAGES[@]}"; do
    dockerfile="$HERE/$name/Dockerfile"
    [ -f "$dockerfile" ] || { echo "unknown image '$name'" >&2; exit 2; }
    tag="${REGISTRY}/mtb-${name}:${TAG}"
    echo
    echo "=== building ${tag} ==="
    # Context is the repo root so analysis/ can be COPYed in; the image carries
    # the analysis code, which is what keeps the pipeline repo pure Nextflow.
    docker build --build-arg GIT_REVISION="$(git -C "$ROOT" rev-parse --short HEAD 2>/dev/null || echo unknown)" \
        -f "$dockerfile" -t "$tag" "$ROOT"
    echo "--- ${tag} built"
done

echo
echo "images:"
docker images --filter "reference=${REGISTRY}/mtb-*" \
    --format '  {{.Repository}}:{{.Tag}}  {{.Size}}'
echo
if [[ "${PUSH:-0}" == "1" ]]; then
    for tag in "${BUILT[@]:-}"; do [[ -n "$tag" ]] && docker push "$tag"; done
else
    echo "Built but not pushed. Set PUSH=1 to publish to $REGISTRY."
fi
