#!/usr/bin/env bash
# Shared helpers for the staged CRyPTIC DuckDB build.
# Sourced by every numbered stage script. Not executable on its own.

set -euo pipefail

# --- argument resolution -----------------------------------------------------
# Each stage script accepts (SRC_DIR, DB_PATH) as positional args, falling back
# to MTB_CRYPTIC_SRC and MTB_CRYPTIC_DB env vars, then the project defaults.
resolve_args() {
    SRC_DIR="${1:-${MTB_CRYPTIC_SRC:-}}"
    DB_PATH="${2:-${MTB_CRYPTIC_DB:-}}"

    if [[ -z "$SRC_DIR" || -z "$DB_PATH" ]]; then
        echo "USAGE: $0 <src_dir> <db_path>" >&2
        echo "   or set MTB_CRYPTIC_SRC and MTB_CRYPTIC_DB env vars" >&2
        exit 2
    fi

    [[ -d "$SRC_DIR" ]] || { echo "ERROR: SRC_DIR not a directory: $SRC_DIR" >&2; exit 2; }

    mkdir -p "$(dirname "$DB_PATH")"
    export SRC_DIR DB_PATH
}

# --- pretty banners + timing -------------------------------------------------
stage_start() {
    local name="$1"
    echo
    echo "================================================================"
    printf "STAGE %s   started at %s\n" "$name" "$(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "================================================================"
    STAGE_T0=$(date +%s)
}

stage_end() {
    local name="$1"
    local t1=$(date +%s)
    local elapsed=$((t1 - STAGE_T0))
    local size
    size=$(du -h "$DB_PATH" 2>/dev/null | awk '{print $1}')
    printf "\nSTAGE %s   done in %ds   db_size=%s\n" "$name" "$elapsed" "$size"
}

# --- duckdb wrapper ----------------------------------------------------------
# Run a heredoc of SQL against the build DB with timer enabled. Stderr is kept
# so SQL errors propagate; stdout shows row counts / timer lines.
run_sql() {
    duckdb "$DB_PATH" "$@"
}

# --- idempotency helper ------------------------------------------------------
# Returns 0 if the named table already exists in the DB.
table_exists() {
    local tbl="$1"
    duckdb "$DB_PATH" -csv -noheader -c "
      SELECT COUNT(*) FROM duckdb_tables()
      WHERE schema_name='main' AND table_name='$tbl'
    " | grep -q '^1$'
}

# --- prerequisite check ------------------------------------------------------
# Confirms the duckdb CLI is on PATH and the source file exists.
require_source_file() {
    local rel="$1"  # e.g. "EFFECTS.parquet"
    [[ -f "$SRC_DIR/$rel" ]] || {
        echo "ERROR: source file missing: $SRC_DIR/$rel" >&2
        exit 3
    }
}
