#!/usr/bin/env bash
# Stage 11 — referential-integrity check views + summary.
#
# DuckDB v1.x parses but does not reliably enforce FOREIGN KEY constraints
# across bulk-loaded fact tables in this size range. Instead of declaring
# (and then ignoring) FK constraints, we materialise verifiable check
# views that surface orphan rows. A clean DB has every v_ri_* view return
# zero rows.
#
# Known cryptic-tables-v3.4.0 RI quirks (not bugs in this build):
#   - ukmyc_growth_to_plates  ⇒ 43,872 rows / 416 distinct UIDs
#   - bashthebug_to_plates    ⇒ 1,190 rows / 214 distinct UIDs
#   These reflect plate-reading rows whose UNIQUEID was never registered in
#   UKMYC_PLATES at the source. Downstream code that joins through
#   ukmyc_plates should LEFT JOIN and filter, not INNER JOIN.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "11_ri_checks"

# `variants` is release-dependent (v3.4.0 ships no VARIANTS.parquet). Omit its
# integrity check rather than report it as passing: a check that silently does
# not run is worse than one that is visibly absent, because the summary would
# read clean while one relation went unexamined.
if table_exists variants; then
    RI_VARIANTS_VIEW="CREATE OR REPLACE VIEW v_ri_variants_to_wgs        AS SELECT row_id, UNIQUEID FROM variants   WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM wgs_samples);"
    RI_VARIANTS_ROW="UNION ALL SELECT 'variants_to_wgs',        COUNT(*) FROM v_ri_variants_to_wgs
"
else
    echo "NOTE: table 'variants' absent — the variants_to_wgs integrity check is SKIPPED, not passed."
    RI_VARIANTS_VIEW=""
    RI_VARIANTS_ROW=""
fi

run_sql <<SQL
.timer on

CREATE OR REPLACE VIEW v_ri_genomes_to_wgs         AS SELECT * FROM genomes      WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM wgs_samples);
CREATE OR REPLACE VIEW v_ri_mutations_to_wgs       AS SELECT row_id, UNIQUEID FROM mutations  WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM wgs_samples);
${RI_VARIANTS_VIEW}
CREATE OR REPLACE VIEW v_ri_effects_to_wgs         AS SELECT row_id, UNIQUEID FROM effects    WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM wgs_samples);
CREATE OR REPLACE VIEW v_ri_predictions_to_wgs     AS SELECT * FROM predictions  WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM wgs_samples);
CREATE OR REPLACE VIEW v_ri_dstmeas_to_dst         AS SELECT row_id, UNIQUEID FROM dst_measurements WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM dst_samples);
CREATE OR REPLACE VIEW v_ri_ukmyc_phen_to_plates   AS SELECT * FROM ukmyc_phenotypes WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM ukmyc_plates);
CREATE OR REPLACE VIEW v_ri_ukmyc_growth_to_plates AS SELECT row_id, UNIQUEID FROM ukmyc_growth WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM ukmyc_plates);
CREATE OR REPLACE VIEW v_ri_bashthebug_to_plates   AS SELECT * FROM bashthebug   WHERE UNIQUEID NOT IN (SELECT UNIQUEID FROM ukmyc_plates);
CREATE OR REPLACE VIEW v_ri_effects_to_drugs       AS SELECT DISTINCT DRUG FROM effects     WHERE DRUG NOT IN (SELECT drug FROM drugs);
CREATE OR REPLACE VIEW v_ri_predictions_to_drugs   AS SELECT DISTINCT DRUG FROM predictions WHERE DRUG NOT IN (SELECT drug FROM drugs);
CREATE OR REPLACE VIEW v_ri_dstmeas_to_drugs       AS SELECT DISTINCT DRUG FROM dst_measurements WHERE DRUG NOT IN (SELECT drug FROM drugs);
CREATE OR REPLACE VIEW v_ri_ukmyc_phen_to_drugs    AS SELECT DISTINCT DRUG FROM ukmyc_phenotypes WHERE DRUG NOT IN (SELECT drug FROM drugs);
CREATE OR REPLACE VIEW v_ri_ukmyc_growth_to_layout AS SELECT DISTINCT PLATEDESIGN, DRUG, DILUTION FROM ukmyc_growth WHERE (PLATEDESIGN, DRUG, DILUTION) NOT IN (SELECT PLATEDESIGN, DRUG, DILUTION FROM plate_layout);

SELECT 'genomes_to_wgs'         AS check_name, COUNT(*) AS orphans FROM v_ri_genomes_to_wgs
UNION ALL SELECT 'mutations_to_wgs',       COUNT(*) FROM v_ri_mutations_to_wgs
${RI_VARIANTS_ROW}UNION ALL SELECT 'effects_to_wgs',         COUNT(*) FROM v_ri_effects_to_wgs
UNION ALL SELECT 'predictions_to_wgs',     COUNT(*) FROM v_ri_predictions_to_wgs
UNION ALL SELECT 'dstmeas_to_dst',         COUNT(*) FROM v_ri_dstmeas_to_dst
UNION ALL SELECT 'ukmyc_phen_to_plates',   COUNT(*) FROM v_ri_ukmyc_phen_to_plates
UNION ALL SELECT 'ukmyc_growth_to_plates', COUNT(*) FROM v_ri_ukmyc_growth_to_plates
UNION ALL SELECT 'bashthebug_to_plates',   COUNT(*) FROM v_ri_bashthebug_to_plates
UNION ALL SELECT 'effects_to_drugs',       COUNT(*) FROM v_ri_effects_to_drugs
UNION ALL SELECT 'predictions_to_drugs',   COUNT(*) FROM v_ri_predictions_to_drugs
UNION ALL SELECT 'dstmeas_to_drugs',       COUNT(*) FROM v_ri_dstmeas_to_drugs
UNION ALL SELECT 'ukmyc_phen_to_drugs',    COUNT(*) FROM v_ri_ukmyc_phen_to_drugs
UNION ALL SELECT 'ukmyc_growth_to_layout', COUNT(*) FROM v_ri_ukmyc_growth_to_layout
ORDER BY orphans DESC, 1;
SQL

stage_end "11_ri_checks"
