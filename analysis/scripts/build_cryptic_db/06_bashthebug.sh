#!/usr/bin/env bash
# Stage 06 — BashTheBug citizen-science classifications.
#
# Loads: bashthebug (307,540), bashthebug_classifications (4,911,886)
#
#   - bashthebug                  — Aggregated reading per (UNIQUEID, DRUG,
#                                    READINGDAY, PLATEIMAGE): consensus MIC,
#                                    confidence, number of valid classifications.
#                                    PK: the four-column composite is unique.
#
#   - bashthebug_classifications  — Raw per-volunteer classifications (4.9M
#                                    rows). PK: classification_id. Linked to
#                                    bashthebug via filename / IMAGEFILENAME.
#                                    No UNIQUEID column — must go through
#                                    bashthebug to reach a sample.
#
# These are mostly used for QC / auditing; they're not on the standard
# resistotyping ML hot path, but downstream consumers need them for image-
# level analyses or for comparing AMyGDA / TMAS / BashTheBug agreement.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "06_bashthebug"

require_source_file "BASHTHEBUG.parquet"
require_source_file "BASHTHEBUG_CLASSIFICATIONS.parquet"

run_sql <<SQL
.timer on

CREATE TABLE bashthebug AS SELECT * FROM read_parquet('$SRC_DIR/BASHTHEBUG.parquet');
ALTER TABLE bashthebug ADD PRIMARY KEY (UNIQUEID, DRUG, READINGDAY, PLATEIMAGE);
CREATE INDEX ix_btb_uid       ON bashthebug(UNIQUEID);
CREATE INDEX ix_btb_filename  ON bashthebug(IMAGEFILENAME);

CREATE TABLE bashthebug_classifications AS
SELECT * FROM read_parquet('$SRC_DIR/BASHTHEBUG_CLASSIFICATIONS.parquet');
ALTER TABLE bashthebug_classifications ADD PRIMARY KEY (classification_id);
CREATE INDEX ix_btbc_filename ON bashthebug_classifications(filename);
CREATE INDEX ix_btbc_drug     ON bashthebug_classifications(drug);

SELECT 'bashthebug' AS tbl, COUNT(*) AS n FROM bashthebug
UNION ALL SELECT 'bashthebug_classifications', COUNT(*) FROM bashthebug_classifications;
SQL

stage_end "06_bashthebug"
