#!/usr/bin/env bash
# Restore the backfilled mlflow.db from MinIO, then serve the MLflow UI.
#
# The FE->metric experiments are produced locally by
# analysis/scripts/feature_mart/tracking.py (backfill from the committed
# result manifests) and pushed to MinIO by bring-up.sh. This entrypoint
# pulls that sqlite db on start. Artifacts (the per-candidate result
# parquets logged via mlflow.log_artifact) live directly in MinIO and are
# served from there — no local copy needed.
#
# Platform notes (from the crypticdb-llm deploy, verified 2026-07):
#   - No persistent volumes -> restore-on-start from MinIO.
#   - The MinIO endpoint var is ABC_MINIO_ENDPOINT (NOT AWS_ENDPOINT_URL);
#     MLflow's S3 artifact client reads MLFLOW_S3_ENDPOINT_URL, so we map it.
#   - Served at ROOT (no --static-prefix) to match the shared-plane Tailscale
#     exposure, which is a direct host:port at root
#     (http://<host>:<port>/), NOT the /apps/<name>/
#     reverse-proxy path. MLflow's SPA works at root, so no static-prefix and
#     no public plane are needed. (MLFLOW_STATIC_PREFIX is left unset; the
#     flag is still wired below in case a proxy-path deploy is ever wanted.)
set -euo pipefail

: "${MLFLOW_DB_BUCKET:?MLFLOW_DB_BUCKET must be set (MinIO bucket holding mlflow.db)}"
: "${MLFLOW_DB_OBJECT:=mlflow/mlflow.db}"
: "${MLFLOW_DB_PATH:=/data/mlflow.db}"
: "${MLFLOW_ARTIFACT_URI:=s3://${MLFLOW_DB_BUCKET}/mlflow/artifacts}"
: "${MLFLOW_STATIC_PREFIX:=}"        # e.g. /apps/mtb-resistotyper-ml-mlflow-tracking
: "${AWS_ACCESS_KEY_ID:?Vault-minted AWS_ACCESS_KEY_ID missing — data: block wrong}"
: "${AWS_SECRET_ACCESS_KEY:?Vault-minted AWS_SECRET_ACCESS_KEY missing}"
: "${ABC_MINIO_ENDPOINT:?Vault-minted ABC_MINIO_ENDPOINT missing}"

# MLflow's boto3 S3 client (for artifacts) reads MLFLOW_S3_ENDPOINT_URL.
export MLFLOW_S3_ENDPOINT_URL="${ABC_MINIO_ENDPOINT}"

mkdir -p "$(dirname "$MLFLOW_DB_PATH")"

echo "[entrypoint] restoring s3://${MLFLOW_DB_BUCKET}/${MLFLOW_DB_OBJECT} -> ${MLFLOW_DB_PATH}"
python - <<PY || echo "[entrypoint] no existing mlflow.db in MinIO — starting empty (backfill + re-push to populate)"
import os, boto3, botocore
s3 = boto3.client("s3", endpoint_url=os.environ["ABC_MINIO_ENDPOINT"])
try:
    s3.download_file(os.environ["MLFLOW_DB_BUCKET"], os.environ["MLFLOW_DB_OBJECT"], os.environ["MLFLOW_DB_PATH"])
    print(f"[entrypoint] restored {os.path.getsize(os.environ['MLFLOW_DB_PATH'])} bytes")
except botocore.exceptions.ClientError as e:
    raise SystemExit(1)
PY

STATIC_ARGS=()
if [[ -n "$MLFLOW_STATIC_PREFIX" ]]; then
  STATIC_ARGS+=(--static-prefix "$MLFLOW_STATIC_PREFIX")
fi

# One worker, deliberately. MLflow 3.x defaults to four uvicorn workers, each
# loading the whole app: all four reached "Application startup complete" and the
# task was gone seconds later, before a single health probe arrived, which is
# what an OOM kill looks like from the outside. A tracking UI serving one
# research group has no use for four workers, and one worker plus headroom is a
# better trade than four workers and a restart loop.
: "${MLFLOW_WORKERS:=1}"

# The snapshot writer runs alongside the server, and the server is NOT exec'd:
# exec would replace this shell, and with it the only process able to take a
# final snapshot when Nomad sends SIGTERM. Losing the last minute of runs on
# every planned restart is avoidable, so it is avoided.
python /snapshot.py &
SNAPSHOT_PID=$!

shutdown() {
  echo "[entrypoint] SIGTERM: flushing the database before exit"
  kill -TERM "$SNAPSHOT_PID" 2>/dev/null || true
  wait "$SNAPSHOT_PID" 2>/dev/null || true
  kill -TERM "$MLFLOW_PID" 2>/dev/null || true
  wait "$MLFLOW_PID" 2>/dev/null || true
  exit 0
}
trap shutdown TERM INT

echo "[entrypoint] starting mlflow server on 0.0.0.0:${PORT:-5000} " \
     "(workers=${MLFLOW_WORKERS}, static-prefix='${MLFLOW_STATIC_PREFIX}')"
mlflow server \
  --host 0.0.0.0 \
  --port "${PORT:-5000}" \
  --workers "${MLFLOW_WORKERS}" \
  --backend-store-uri "sqlite:///${MLFLOW_DB_PATH}" \
  --artifacts-destination "$MLFLOW_ARTIFACT_URI" \
  "${STATIC_ARGS[@]}" &
MLFLOW_PID=$!
wait "$MLFLOW_PID"
