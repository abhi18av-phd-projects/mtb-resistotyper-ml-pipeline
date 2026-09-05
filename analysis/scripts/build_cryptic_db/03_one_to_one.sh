#!/usr/bin/env bash
# Stage 03 — 1:1 sample-attribute tables.
#
# Loads: genomes (54,057), ukmyc_plates (26,768)
#
# Both have UNIQUEID as a clean primary key. Each sequenced sample has at most
# one row in each:
#
#   - genomes        — species/lineage/coverage stats from Mykrobe/Clockwork.
#                      54,057 rows for 54,565 sequenced samples ⇒ a few hundred
#                      WGS samples lack a genome record (typically the
#                      "cannot speciate" or "cannot assemble" subset called out
#                      in the v3.x release notes).
#
#   - ukmyc_plates   — provenance of an isolate's UKMYC5/6 microtitre plate
#                      image. Only ~26,768 of the sequenced samples have one,
#                      because not every sample was tested on a CRyPTIC plate
#                      (the rest only have clinical DST in dst_measurements).
#
# Both link to wgs_samples by UNIQUEID. ukmyc_plates is itself the parent of
# ukmyc_phenotypes, ukmyc_growth, and bashthebug.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "03_one_to_one"

require_source_file "GENOMES.parquet"
require_source_file "UKMYC_PLATES.parquet"

run_sql <<SQL
.timer on

CREATE TABLE genomes AS SELECT * FROM read_parquet('$SRC_DIR/GENOMES.parquet');
ALTER TABLE genomes ADD PRIMARY KEY (UNIQUEID);

CREATE TABLE ukmyc_plates AS SELECT * FROM read_parquet('$SRC_DIR/UKMYC_PLATES.parquet');
ALTER TABLE ukmyc_plates ADD PRIMARY KEY (UNIQUEID);
CREATE INDEX ix_ukmycp_image ON ukmyc_plates(IMAGEFILENAME);

SELECT 'genomes' AS tbl, COUNT(*) AS n FROM genomes
UNION ALL SELECT 'ukmyc_plates', COUNT(*) FROM ukmyc_plates;
SQL

stage_end "03_one_to_one"
