#!/usr/bin/env bash
# Stage 12 — build metadata table.
#
# Records build provenance inside the database itself so the file is self-
# describing: source release version, build timestamp, DuckDB version, and a
# free-form RI status string. Pipeline runs that read this DB can audit which
# release version they're working against without external metadata files.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "12_metadata"

# Allow the caller to inject a source-release version. Defaults to whatever
# the source directory is named, which conventionally encodes the version
# (e.g. "cryptic-tables-v3.4.0").
SOURCE_VERSION="${MTB_CRYPTIC_VERSION:-$(basename "$SRC_DIR")}"

run_sql <<SQL
CREATE OR REPLACE TABLE _database_metadata (
    key   VARCHAR PRIMARY KEY,
    value VARCHAR
);
INSERT INTO _database_metadata VALUES
    ('source_version',      '$SOURCE_VERSION'),
    ('source_dir',          '$SRC_DIR'),
    ('build_timestamp_utc', strftime('%Y-%m-%dT%H:%M:%SZ', current_timestamp::TIMESTAMPTZ AT TIME ZONE 'UTC')),
    ('duckdb_version',      version()),
    ('build_method',        'analysis/scripts/build_cryptic_db/run_all.sh');

SELECT * FROM _database_metadata ORDER BY key;
SQL

stage_end "12_metadata"
