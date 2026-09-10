#!/usr/bin/env bash
# Train and package a deployable bundle set for ONE CRyPTIC release.
#
# The webapp's release selector can only offer what is installed, and a
# comparison needs two. This produces the second and third sets -- the full
# v3.4.0 build and v2.1.2 -- against the same code that produced the shipped
# v3.4.0-slim one, so a difference between them is a difference in the data
# rather than in the pipeline.
#
# Runs the pipeline's own containers rather than a host environment: the host
# has no numpy, and a bundle trained in an environment nobody can reconstruct
# is not reproducible whatever its AUC says. The repository's analysis/ is
# bind-mounted OVER the image's copy, because the images bake analysis/ in and
# the release-labelling fix is newer than every published tag.
#
# Usage:  train-release-bundles.sh <db-path> <workdir>
set -euo pipefail

DB="${1:?usage: train-release-bundles.sh <db> <workdir>}"
WORK="${2:?usage: train-release-bundles.sh <db> <workdir>}"
REPO="${REPO:-$HOME/mtb-pipeline-src}"
REG="ghcr.io/abhi18av-phd-projects/mtb-resistotyper-ml-pipeline"
FE="$REG/mtb-fe:v0.4.2"
CAUSAL="$REG/mtb-causal:v0.3.0"
DRUGS="${DRUGS:-RIF INH EMB MXF LEV ETH KAN AMI}"
POOL="${POOL:-300}"

mkdir -p "$WORK/marts" "$WORK/causal" "$WORK/models" "$WORK/release" "$WORK/logs"
DBDIR="$(cd "$(dirname "$DB")" && pwd)"; DBFILE="$(basename "$DB")"

run() {  # run <image> <command...>
  local img="$1"; shift
  docker run --rm \
    -v "$REPO:/opt/mtb" -v "$DBDIR:/db:ro" -v "$WORK:/work" \
    -w /opt/mtb "$img" "$@"
}

echo "==> [1/4] feature marts"
for d in $DRUGS; do
  if ls "$WORK"/marts/feature_mart_"$d"_*_vFULL.parquet >/dev/null 2>&1; then
    echo "    $d mart present, skipping"; continue
  fi
  echo "    $d"
  # --version vFULL is not cosmetic: train_logistic globs for "*vFULL.parquet",
  # so a mart carrying build_mart's default version label is invisible to it and
  # the drug is silently skipped with no error anywhere.
  run "$FE" python -m analysis.scripts.feature_mart.build_mart \
      --drug "$d" --db "/db/$DBFILE" --out /work/marts --version vFULL \
      > "$WORK/logs/mart_$d.log" 2>&1
done

echo "==> [2/4] causal concordance (deployment arm, refit on all)"
for d in $DRUGS; do
  out="$WORK/causal/$d"
  if [ -f "$out/causal_concordance_${d}_mutation_level.parquet" ]; then
    echo "    $d concordance present, skipping"; continue
  fi
  mart="$(ls "$WORK"/marts/feature_mart_"$d"_*_vFULL.parquet 2>/dev/null | head -1)"
  [ -n "$mart" ] || { echo "    $d has no mart, skipping"; continue; }
  mkdir -p "$out"
  echo "    $d"
  run "$CAUSAL" python -m analysis.scripts.feature_mart.causal \
      --mart "/work/marts/$(basename "$mart")" --drug "$d" --level mutation \
      --db "/db/$DBFILE" --pool-size "$POOL" --out "/work/causal/$d" \
      --refit-on-all \
      > "$WORK/logs/causal_$d.log" 2>&1
done

echo "==> [3/4] tuned L2-logistic"
run "$CAUSAL" python -m analysis.scripts.feature_mart.train_logistic \
    --marts /work/marts --causal /work/causal --out /work/models \
    --baseline /work/none.json 2>&1 | tee "$WORK/logs/train_logistic.log"

echo "==> [4/4] package a depositable release"
# One self-contained directory per release -- models, cards, schemas, training
# and selection provenance, licence, citation, Zenodo metadata, checksums -- that
# can be zipped and deposited as one record, and that the webapp installs
# unchanged. The pipeline revision is recorded from $REVISION when the source
# tree carries no .git (an rsync'd copy does not).
REVISION="${REVISION:-$(git -C "$REPO" rev-parse --short HEAD 2>/dev/null || echo unrecorded)}"
run "$CAUSAL" python -m analysis.scripts.feature_mart.package_bundles \
    --models /work/models --marts /work/marts --causal /work/causal \
    --out /work/release --licence-file /opt/mtb/LICENSE \
    --pipeline-revision "$REVISION" --images "$FE, $CAUSAL" \
    2>&1 | tee "$WORK/logs/package.log"

echo "==> done: $(ls -d "$WORK"/release/mtb-resistotyper-ml-models-* 2>/dev/null | head -1)"
