#!/usr/bin/env bash
# Sequential orchestrator — runs every numbered stage in order.
#
# Usage:
#   ./run_all.sh <src_dir> <db_path>
#   MTB_CRYPTIC_SRC=... MTB_CRYPTIC_DB=... ./run_all.sh
#
# For Nextflow integration, prefer running stages as separate processes so
# each one is independently cacheable and resumable. This wrapper is for
# interactive / one-shot builds.
#
# Total expected wall time: ~14 minutes on Apple-silicon laptop.
# Final DB size: ~19 GB.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

SRC_DIR="${1:-${MTB_CRYPTIC_SRC:-}}"
DB_PATH="${2:-${MTB_CRYPTIC_DB:-}}"

if [[ -z "$SRC_DIR" || -z "$DB_PATH" ]]; then
    cat <<EOF >&2
USAGE: $0 <src_dir> <db_path>

Where <src_dir> is the cryptic-tables-vX.Y.Z directory containing the
parquet/CSV files, and <db_path> is the output .duckdb file path.

You can also export MTB_CRYPTIC_SRC and MTB_CRYPTIC_DB instead of
passing positional args.
EOF
    exit 2
fi

T0=$(date +%s)
echo "============================================================"
echo "CRyPTIC DuckDB build  —  $(date -u +%Y-%m-%dT%H:%M:%SZ)"
echo "  src: $SRC_DIR"
echo "  db:  $DB_PATH"
echo "============================================================"

for stage in \
    00_init_db.sh \
    01_reference.sh \
    02_sample_anchors.sh \
    03_one_to_one.sh \
    04_effects_predictions.sh \
    05_phenotypes.sh \
    06_bashthebug.sh \
    07_ukmyc_growth.sh \
    08_mutations.sh \
    09_variants.sh \
    10_views.sh \
    11_ri_checks.sh \
    12_metadata.sh
do
    bash "$SCRIPT_DIR/$stage" "$SRC_DIR" "$DB_PATH"
done

T1=$(date +%s)
echo
echo "============================================================"
echo "BUILD COMPLETE  —  $((T1 - T0))s total"
echo "  $(ls -lh "$DB_PATH" | awk '{print $5"  "$9}')"
echo "============================================================"
