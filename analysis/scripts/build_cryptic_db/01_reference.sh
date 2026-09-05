#!/usr/bin/env bash
# Stage 01 — reference / lookup tables.
#
# Loads: drugs (39), countries (243), sites (31), plate_layout (327)
# These are tiny, deterministic, and define the value spaces every fact table
# joins back to. Always loaded first so PK/FK relationships make sense.
#
# Notes on schema decisions:
#   - drugs.drug         — three-letter code, the join key for every drug
#                          column elsewhere (effects.DRUG, predictions.DRUG, etc.)
#   - countries.code_3   — ISO 3-letter, the cleanest country FK target
#   - sites.site_id      — string ID encoded inside UNIQUEID (e.g. "01"..."13")
#   - plate_layout.well_id (surrogate) — the natural composite
#       (PLATEDESIGN, DRUG, ROW, COL) is a true PK only when ROW/COL are not
#       NULL; positive-control wells (DRUG='POS') and the >/< boundary rows
#       break that, so we use a surrogate and add indexes for the natural key.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "01_reference"

require_source_file "DRUG_CODES.csv"
require_source_file "COUNTRIES_LOOKUP.csv"
require_source_file "SITES.csv"
require_source_file "PLATE_LAYOUT.parquet"

run_sql <<SQL
.timer on

CREATE TABLE drugs AS
SELECT DRUG_3_LETTER_CODE AS drug, DRUG_NAME AS drug_name
FROM read_csv_auto('$SRC_DIR/DRUG_CODES.csv');
ALTER TABLE drugs ADD PRIMARY KEY (drug);

CREATE TABLE countries AS
SELECT COUNTRY_NAME           AS country_name,
       COUNTRY_CODE_2_LETTER  AS code_2,
       COUNTRY_CODE_3_LETTER  AS code_3,
       COUNTRY_CODE_NUMERIC   AS code_numeric,
       LAT                    AS latitude,
       LONG                   AS longitude
FROM read_csv_auto('$SRC_DIR/COUNTRIES_LOOKUP.csv');
ALTER TABLE countries ADD PRIMARY KEY (code_3);

CREATE TABLE sites AS
SELECT SITEID                AS site_id,
       DESCRIPTION            AS description,
       ABBREVIATION           AS abbreviation,
       CITY                   AS city,
       COUNTRY                AS country_name,
       COUNTRY_3_LETTER_CODE  AS country_code_3,
       LONG                   AS longitude,
       LAT                    AS latitude
FROM read_csv_auto('$SRC_DIR/SITES.csv');
ALTER TABLE sites ADD PRIMARY KEY (site_id);

CREATE TABLE plate_layout AS
SELECT row_number() OVER () AS well_id,
       PLATEDESIGN, DRUG, DILUTION, CONC, ROW, COL, BINARY_PHENOTYPE
FROM read_parquet('$SRC_DIR/PLATE_LAYOUT.parquet');
ALTER TABLE plate_layout ADD PRIMARY KEY (well_id);
CREATE INDEX ix_plate_layout_design_drug_dil ON plate_layout(PLATEDESIGN, DRUG, DILUTION);
CREATE INDEX ix_plate_layout_design_drug_rc  ON plate_layout(PLATEDESIGN, DRUG, ROW, COL);

SELECT 'drugs'        AS tbl, COUNT(*) AS n FROM drugs
UNION ALL SELECT 'countries',    COUNT(*) FROM countries
UNION ALL SELECT 'sites',        COUNT(*) FROM sites
UNION ALL SELECT 'plate_layout', COUNT(*) FROM plate_layout;
SQL

stage_end "01_reference"
