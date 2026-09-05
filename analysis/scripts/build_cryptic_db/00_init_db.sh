#!/usr/bin/env bash
# Stage 00 — initialise an empty DuckDB file at $DB_PATH.
#
# Removes any existing file so the build is fully reproducible. Skip this
# stage if you want incremental top-ups; later stages will fail loudly if a
# table already exists.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "00_init_db"

command -v duckdb >/dev/null || { echo "ERROR: duckdb CLI not on PATH" >&2; exit 1; }

if [[ -f "$DB_PATH" ]]; then
    echo "Removing existing $DB_PATH"
    rm -f "$DB_PATH"
fi

run_sql "SELECT version() AS duckdb_version;"
echo "Initialised empty DB at $DB_PATH"

stage_end "00_init_db"
