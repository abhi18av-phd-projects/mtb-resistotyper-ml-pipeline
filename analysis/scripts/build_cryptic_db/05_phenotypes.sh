#!/usr/bin/env bash
# Stage 05 — phenotype tables (clinical DST + UKMYC plate-derived).
#
# Loads: dst_measurements (660,961), ukmyc_phenotypes (288,904)
#
# Two parallel measurement streams:
#
#   - dst_measurements — clinical DST results: MIC + binary R/S derived from
#                        broth microdilution, MGIT960, LJ, MYCOTB, and other
#                        platforms. Multiple measurements per (UNIQUEID, DRUG)
#                        are common (different methods, different labs).
#                        FK: UNIQUEID → dst_samples.
#                        PK: surrogate row_id; (UNIQUEID, DRUG, SOURCE) is
#                        a near-key but not strictly unique once methods
#                        + quality flags are considered.
#
#   - ukmyc_phenotypes — One row per (UNIQUEID, DRUG) for samples that were
#                        tested on the CRyPTIC UKMYC5/6 microtitre plate,
#                        with both quantitative (MIC, log2MIC) and binary
#                        (R/S) interpretations plus quality flags.
#                        FK: UNIQUEID → ukmyc_plates.
#                        PK: (UNIQUEID, DRUG) is strictly unique.
#
# Use ukmyc_phenotypes when you want the high-quality CRyPTIC MICs;
# use dst_measurements when you need broader clinical coverage including
# samples that don't have a UKMYC plate.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "05_phenotypes"

require_source_file "DST_MEASUREMENTS.parquet"
require_source_file "UKMYC_PHENOTYPES.parquet"

run_sql <<SQL
.timer on

CREATE TABLE dst_measurements AS
SELECT row_number() OVER () AS row_id, *
FROM read_parquet('$SRC_DIR/DST_MEASUREMENTS.parquet');
ALTER TABLE dst_measurements ADD PRIMARY KEY (row_id);
CREATE INDEX ix_dstm_uid       ON dst_measurements(UNIQUEID);
CREATE INDEX ix_dstm_uid_drug  ON dst_measurements(UNIQUEID, DRUG);
CREATE INDEX ix_dstm_drug      ON dst_measurements(DRUG);

CREATE TABLE ukmyc_phenotypes AS SELECT * FROM read_parquet('$SRC_DIR/UKMYC_PHENOTYPES.parquet');
ALTER TABLE ukmyc_phenotypes ADD PRIMARY KEY (UNIQUEID, DRUG);
CREATE INDEX ix_ukmycph_drug   ON ukmyc_phenotypes(DRUG);
CREATE INDEX ix_ukmycph_design ON ukmyc_phenotypes(PLATEDESIGN);

SELECT 'dst_measurements' AS tbl, COUNT(*) AS n FROM dst_measurements
UNION ALL SELECT 'ukmyc_phenotypes', COUNT(*) FROM ukmyc_phenotypes;
SQL

stage_end "05_phenotypes"
