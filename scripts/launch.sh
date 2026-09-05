#!/usr/bin/env bash
# Launch one experiment arm on the abc cluster, with the four corrections a bare
# `abc pipeline run` needs on this deployment. Each exists because omitting it
# fails in a way that looks like something else.
#
#   --plugin nf-nomad@…      abc-cli otherwise resolves 0.4.4 regardless of what
#                            the nomad profile pins, so the profile's pin is
#                            silently ignored and the run gets an older executor.
#
#   --head-pool / --worker-pool
#                            from site.env. The WORKER pool is the one that
#                            matters. The context sets worker_pool=compute, and
#                            this cluster's only client is in platform, so every
#                            worker job is accepted and then never placed:
#                            "No nodes were eligible for evaluation".
#
#   NF_MINIO_ENDPOINT        without it the head is handed the PUBLIC edge for
#                            object storage, so a database sitting on the same
#                            physical host is fetched Stellenbosch -> Johannesburg
#                            -> Stellenbosch. Measured: 0.4 MB/s against 85 MB/s,
#                            a ~200x difference, on a 10 GB file.
#
# Lives in scripts/, NOT bin/. A bin/ directory at the project root is staged by
# Nextflow for every task, and on Nextflow 26.04.3 the s5cmd work-dir provider
# cannot enumerate it (session.getBinDirs() has no signature), falls back to
# uploading it whole, and the allocation fails. nf-nomad then reads the absent
# .exitcode as 0 and reports the failure as a task that succeeded and produced
# nothing. The project therefore has no bin/ at all.
#
# Usage: scripts/launch.sh experiments/<arm>.yml [extra abc flags...]
set -euo pipefail

ARM="${1:?usage: bin/launch.sh experiments/<arm>.yml [extra flags]}"; shift || true
[ -f "$ARM" ] || { echo "no such params file: $ARM" >&2; exit 1; }

# Site settings live in site.env (gitignored; see site.env.example). Sourced
# here rather than hardcoded, so this script names no cluster.
HERE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
if [ -f "$HERE_DIR/site.env" ]; then set -a; . "$HERE_DIR/site.env"; set +a; fi

REPO="${MTB_PIPELINE_REPO:-https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-pipeline}"
NAME="$(basename "$ARM" .yml | tr '.' '-')"
: "${ABC_CLI_CONTEXT:?set ABC_CLI_CONTEXT in site.env}"; export ABC_CLI_CONTEXT
: "${MTB_S3_ENDPOINT:?set MTB_S3_ENDPOINT in site.env — the INTERNAL object-store address}"
NF_MINIO_ENDPOINT="$MTB_S3_ENDPOINT"; export NF_MINIO_ENDPOINT
: "${NF_NOMAD_PLUGIN:=nf-nomad@0.5.0-edge7}"
: "${MTB_GROUP_BUCKET:?set MTB_GROUP_BUCKET in site.env}"
: "${WORK_DIR:=s3://${MTB_GROUP_BUCKET}/pipelines/mtb-resistotyper-ml/workdir/}"

# A private repo needs a clone credential. Read from the environment rather than
# taken as an argument: --git-token would put it on the command line, and it
# lands in the Nomad job spec either way, which is worth knowing about.
if [ -z "${GITHUB_TOKEN:-}" ] && command -v gh >/dev/null 2>&1; then
  GITHUB_TOKEN="$(gh auth token 2>/dev/null || true)"; export GITHUB_TOKEN
fi
[ -n "${GITHUB_TOKEN:-}" ] || echo "warning: no GITHUB_TOKEN; a private repo will fail to clone" >&2

REV="${REVISION:-$(git rev-parse HEAD)}"
echo "arm      $ARM"
echo "revision $REV"
echo "plugin   $NF_NOMAD_PLUGIN"
echo "s3       $NF_MINIO_ENDPOINT"

exec abc pipeline run "$REPO" \
  --params-file "$ARM" \
  --profile nomad,containers \
  --revision "$REV" \
  `# nextflow.config resolves these with System.getenv INSIDE THE HEAD JOB, so
   # sourcing site.env here is not enough: without forwarding them the config
   # falls back to its placeholder defaults and the run reads
   # s3://su-example-group/, which does not exist. That is the intended failure
   # mode for an unset variable, but it has to be forwarded to be set.` \
  --env GITHUB_TOKEN --env NF_MINIO_ENDPOINT \
  --env MTB_GROUP_BUCKET --env MTB_NOMAD_NAMESPACE --env MTB_REGISTRY \
  --env MTB_MLFLOW_URI \
  --head-pool "${MTB_HEAD_POOL:-platform}" --worker-pool "${MTB_WORKER_POOL:-platform}" \
  --plugin "$NF_NOMAD_PLUGIN" --plugin "${NF_S5CMD_PLUGIN:=nf-nomad-s5cmd@0.1.8}" \
  --work-dir "$WORK_DIR" \
  --name "$NAME" "$@"
