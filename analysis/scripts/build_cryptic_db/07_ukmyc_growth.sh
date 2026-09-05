#!/usr/bin/env bash
# Stage 07 — UKMYC well-level growth measurements.
#
# Loads: ukmyc_growth (3,193,728)
#
# One row per (UNIQUEID, DRUG, READINGDAY, DILUTION, SITEID, PLATEDESIGN, ...)
# describing the growth measurement at a specific dilution well on a specific
# day for a specific sample. ~33k near-duplicates exist (plates re-imaged on
# different days), so we use a surrogate row_id PK.
#
# Joins:
#   - UNIQUEID                                   → ukmyc_plates
#   - (PLATEDESIGN, DRUG, DILUTION)              → plate_layout (gives CONC,
#                                                  ROW, COL, BINARY_PHENOTYPE)
#
# Use this table when you need raw growth values rather than the aggregated
# MIC in ukmyc_phenotypes — e.g., for image-feature ML or growth-curve
# modelling.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "07_ukmyc_growth"

require_source_file "UKMYC_GROWTH.parquet"

run_sql <<SQL
.timer on

CREATE TABLE ukmyc_growth AS
SELECT row_number() OVER () AS row_id, *
FROM read_parquet('$SRC_DIR/UKMYC_GROWTH.parquet');
ALTER TABLE ukmyc_growth ADD PRIMARY KEY (row_id);
CREATE INDEX ix_ukmycg_uid_drug ON ukmyc_growth(UNIQUEID, DRUG);
CREATE INDEX ix_ukmycg_layout   ON ukmyc_growth(PLATEDESIGN, DRUG, DILUTION);

SELECT COUNT(*) AS n FROM ukmyc_growth;
SQL

stage_end "07_ukmyc_growth"
