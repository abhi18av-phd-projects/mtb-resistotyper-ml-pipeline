#!/usr/bin/env bash
# Stage 02 — sample anchor tables.
#
# Loads: wgs_samples (54,565), dst_samples (65,842)
#
# These are the two parent tables every fact table joins to. They are kept
# separate because the populations are different:
#
#   - wgs_samples = isolates that were *sequenced* through Pathogena/Clockwork.
#     Every genetic table (mutations, variants, effects, predictions, genomes)
#     has UNIQUEID values that are a subset of wgs_samples.UNIQUEID.
#
#   - dst_samples = isolates that have *clinical drug-susceptibility test*
#     measurements (broth microdilution, MGIT, LJ, etc.). dst_measurements
#     is a child of this table.
#
# The two populations overlap heavily but neither contains the other.
# 54,354 UNIQUEIDs appear in both. ML pipelines that need genotype + phenotype
# implicitly take this intersection via INNER JOIN.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "02_sample_anchors"

require_source_file "WGS_SAMPLES.parquet"
require_source_file "DST_SAMPLES.parquet"

run_sql <<SQL
.timer on

CREATE TABLE wgs_samples AS SELECT * FROM read_parquet('$SRC_DIR/WGS_SAMPLES.parquet');
ALTER TABLE wgs_samples ADD PRIMARY KEY (UNIQUEID);
-- run_accession is also unique (each row of WGS_SAMPLES is one ENA run)
CREATE UNIQUE INDEX ux_wgs_run_accession ON wgs_samples(run_accession);

CREATE TABLE dst_samples AS SELECT * FROM read_parquet('$SRC_DIR/DST_SAMPLES.parquet');
ALTER TABLE dst_samples ADD PRIMARY KEY (UNIQUEID);

SELECT 'wgs_samples' AS tbl, COUNT(*) AS n FROM wgs_samples
UNION ALL SELECT 'dst_samples', COUNT(*) FROM dst_samples;
SQL

stage_end "02_sample_anchors"
