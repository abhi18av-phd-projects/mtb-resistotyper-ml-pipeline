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
# The release this database was built from. Falling back to the source dir's
# BASENAME produced "src" for every build, because the convention is
# <build-dir>/src -- so the one field that distinguishes two releases was
# identical in both, and a downstream mart could not tell them apart.
#
# Order: an explicit MTB_CRYPTIC_VERSION, then a cryptic-tables-vX.Y.Z component
# anywhere in the path, then a parent directory that looks like a version. If
# none of those yields something usable the build FAILS: a database that cannot
# name its release will produce marts that cannot be compared, and finding that
# out at mart-build time is far more expensive than here.
_derive_release() {
    [[ -n "${MTB_CRYPTIC_VERSION:-}" ]] && { echo "$MTB_CRYPTIC_VERSION"; return; }
    local p; p="$(cd "$SRC_DIR" && pwd)"
    local m; m="$(printf '%s' "$p" | grep -oE 'cryptic-tables-v[0-9]+(\.[0-9]+)*' | tail -1)"
    [[ -n "$m" ]] && { echo "${m#cryptic-tables-}"; return; }
    local base; base="$(basename "$p")"
    if [[ "$base" =~ ^(src|source|data|\.)$ ]]; then
        base="$(basename "$(dirname "$p")")"
    fi
    m="$(printf '%s' "$base" | grep -oE 'v?[0-9]+\.[0-9]+(\.[0-9]+)?' | tail -1)"
    [[ -n "$m" ]] && { echo "$m"; return; }
    echo ""
}
SOURCE_VERSION="$(_derive_release)"
if [[ -z "$SOURCE_VERSION" ]]; then
    echo "ERROR: cannot determine the CRyPTIC release for $SRC_DIR." >&2
    echo "       Set MTB_CRYPTIC_VERSION, or name the source dir after the release." >&2
    exit 4
fi
echo "  release: $SOURCE_VERSION"

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
