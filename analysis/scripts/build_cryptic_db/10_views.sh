#!/usr/bin/env bash
# Stage 10 — convenience views.
#
# Creates three views that downstream notebooks rely on:
#
#   - v_catalogues       — distinct (CATALOGUE_NAME, CATALOGUE_VERSION) tuples
#                          across effects + predictions, with row counts.
#                          Use this to discover which catalogues are present
#                          in the current release without grepping the docs.
#
#   - v_sample_coverage  — one row per UNIQUEID with boolean flags for each
#                          downstream table (has_wgs, has_dst, has_genome,
#                          has_mutations, ...). Lets you compute the
#                          intersection of any combination of data types in
#                          one query.
#
#   - v_sample_fanout    — long-format (table_name, UNIQUEID, n_rows) showing
#                          how many child rows each sample has in each fact
#                          table. Useful for QC and for picking sample
#                          subsets with adequate signal.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "10_views"

# `variants` is release-dependent: v3.4.0 ships no VARIANTS.parquet, so stage 09
# does not run and the table is absent. Omit the two references rather than
# fabricate an empty table -- an empty one would report has_variants = FALSE for
# every isolate and silently drop the fanout row, which is a wrong answer
# wearing the shape of a right one.
if table_exists variants; then
    HAS_VARIANTS_EXPR="EXISTS (SELECT 1 FROM variants    v WHERE v.UNIQUEID = coalesce(w.UNIQUEID, d.UNIQUEID)) AS has_variants,"
    VARIANTS_FANOUT_BRANCH="UNION ALL SELECT 'variants',         UNIQUEID, COUNT(*) FROM variants         GROUP BY UNIQUEID
"
else
    echo "NOTE: table 'variants' absent (no VARIANTS.parquet in this release) —"
    echo "      v_sample_coverage omits has_variants; v_sample_fanout omits its row."
    HAS_VARIANTS_EXPR=""
    VARIANTS_FANOUT_BRANCH=""
fi

run_sql <<SQL
.timer on

CREATE OR REPLACE VIEW v_catalogues AS
SELECT CATALOGUE_NAME, CATALOGUE_VERSION,
       COUNT(*) FILTER (WHERE source = 'effects')     AS n_effects_rows,
       COUNT(*) FILTER (WHERE source = 'predictions') AS n_predictions_rows
FROM (SELECT 'effects'     AS source, CATALOGUE_NAME, CATALOGUE_VERSION FROM effects
      UNION ALL
      SELECT 'predictions' AS source, CATALOGUE_NAME, CATALOGUE_VERSION FROM predictions)
GROUP BY 1, 2;

CREATE OR REPLACE VIEW v_sample_coverage AS
SELECT coalesce(w.UNIQUEID, d.UNIQUEID) AS UNIQUEID,
       w.UNIQUEID  IS NOT NULL  AS has_wgs,
       d.UNIQUEID  IS NOT NULL  AS has_dst,
       g.UNIQUEID  IS NOT NULL  AS has_genome,
       up.UNIQUEID IS NOT NULL  AS has_ukmyc_plate,
       EXISTS (SELECT 1 FROM mutations   m WHERE m.UNIQUEID = coalesce(w.UNIQUEID, d.UNIQUEID)) AS has_mutations,
       ${HAS_VARIANTS_EXPR}
       EXISTS (SELECT 1 FROM effects     e WHERE e.UNIQUEID = coalesce(w.UNIQUEID, d.UNIQUEID)) AS has_effects,
       EXISTS (SELECT 1 FROM predictions p WHERE p.UNIQUEID = coalesce(w.UNIQUEID, d.UNIQUEID)) AS has_predictions
FROM wgs_samples w
FULL OUTER JOIN dst_samples d USING (UNIQUEID)
LEFT  JOIN genomes      g  USING (UNIQUEID)
LEFT  JOIN ukmyc_plates up USING (UNIQUEID);

CREATE OR REPLACE VIEW v_sample_fanout AS
SELECT 'mutations'        AS table_name, UNIQUEID, COUNT(*) AS n FROM mutations         GROUP BY UNIQUEID
${VARIANTS_FANOUT_BRANCH}UNION ALL SELECT 'effects',          UNIQUEID, COUNT(*) FROM effects          GROUP BY UNIQUEID
UNION ALL SELECT 'predictions',      UNIQUEID, COUNT(*) FROM predictions      GROUP BY UNIQUEID
UNION ALL SELECT 'dst_measurements', UNIQUEID, COUNT(*) FROM dst_measurements GROUP BY UNIQUEID
UNION ALL SELECT 'ukmyc_phenotypes', UNIQUEID, COUNT(*) FROM ukmyc_phenotypes GROUP BY UNIQUEID
UNION ALL SELECT 'ukmyc_growth',     UNIQUEID, COUNT(*) FROM ukmyc_growth     GROUP BY UNIQUEID
UNION ALL SELECT 'bashthebug',       UNIQUEID, COUNT(*) FROM bashthebug       GROUP BY UNIQUEID;

SELECT * FROM v_catalogues ORDER BY 1, 2;
SQL

stage_end "10_views"
