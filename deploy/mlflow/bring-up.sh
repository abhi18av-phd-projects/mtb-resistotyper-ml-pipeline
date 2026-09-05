#!/usr/bin/env bash
# Build the MLflow tracking server on $HOST and deploy it as an abc-app.
#
# Built ON $HOST because "$HOST.local/..." names an image already present in
# that host's docker daemon; it is not a registry, and `docker push` to it always
# fails. Build elsewhere and the deploy places a job that can never pull.
#
# Per-command context via ABC_CLI_CONTEXT — the global active context is never
# modified (CLAUDE.md rule: other sessions share it).
set -euo pipefail

HERE_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
if [ -f "$HERE_ROOT/site.env" ]; then set -a; . "$HERE_ROOT/site.env"; set +a; fi
: "${ABC_CLI_CONTEXT:?set ABC_CLI_CONTEXT in site.env}"
CTX="$ABC_CLI_CONTEXT"
HOST="${MTB_BUILD_HOST:?set MTB_BUILD_HOST in site.env — the host with a docker daemon}"
TAG="${MLFLOW_TAG:-v0.2.0}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

echo "==> shipping the build context to $HOST"
ssh "$HOST" 'rm -rf ~/mlflow-build && mkdir -p ~/mlflow-build'
scp -q "$HERE/Dockerfile" "$HERE/entrypoint.sh" "$HOST:~/mlflow-build/"

echo "==> building $HOST.local/mlflow-tracking:$TAG"
ssh "$HOST" "cd ~/mlflow-build && docker build -t $HOST.local/mlflow-tracking:$TAG ."

echo "==> deploying"
cd "$HERE"
# Render the descriptor from the template using site.env.
envsubst < "$HERE/abc-app.yaml.template" > "$HERE/abc-app.yaml"
ABC_CLI_CONTEXT="$CTX" abc app validate
# --node-pool is not optional: this context supplies no head_pool, so without it
# the job carries no node_pool, lands in the empty "default" pool, and is
# accepted but never placed ("No nodes were eligible for evaluation"). $HOST is
# the platform pool, which is also where the image lives.
ABC_CLI_CONTEXT="$CTX" abc app deploy --node-pool "${MTB_HEAD_POOL:-platform}" --health-timeout 5m
ABC_CLI_CONTEXT="$CTX" abc app show mlflow
