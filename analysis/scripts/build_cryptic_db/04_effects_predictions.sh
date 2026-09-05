#!/usr/bin/env bash
# Stage 04 — catalogue-based prediction tables.
#
# Loads: effects (1,154,127), predictions (810,615)
#
# These two are the *output* of running each sample's variant set through a
# WHO resistance catalogue (`gnomonicus` → `piezo`). They are the typical
# "labels" or "ground-truth" used in resistotyping ML.
#
#   - predictions   — ONE row per (UNIQUEID, DRUG, CATALOGUE_NAME, CATALOGUE_VERSION).
#                     PREDICTION ∈ {S, R, U, F}: the drug-level call. This is
#                     the primary classification target.
#
#   - effects       — Why each prediction is what it is. ONE row per
#                     (UNIQUEID, DRUG, GENE, MUTATION, catalogue): each gene/
#                     mutation that contributed evidence to the drug call.
#                     Fan-out is ~21 rows per sample.
#
# Catalogue versioning is critical:
#   - effects/predictions can carry results from MULTIPLE catalogue versions
#     simultaneously (WHO1 + WHO2 historically). v3.4.0 has only WHO2 right
#     now, but every analysis SHOULD filter explicitly on
#       (CATALOGUE_NAME='WHO-UCN-GTB-PCI-2023.5' AND CATALOGUE_VERSION=2.0)
#     so it survives a future release that re-adds WHO1.
#
# Schema decisions:
#   - predictions PK = (UNIQUEID, DRUG, CATALOGUE_NAME, CATALOGUE_VERSION) is
#     a true unique composite, used as the PK directly.
#   - effects gets a surrogate row_id PK because GENE can be NULL for
#     multi-gene mutations (e.g. "fabG1@G241G&fabG1@c723t"), which DuckDB
#     forbids in a PK column. The natural near-key is the same six columns,
#     and (UNIQUEID, DRUG) is indexed for join performance.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "04_effects_predictions"

require_source_file "EFFECTS.parquet"
require_source_file "PREDICTIONS.parquet"

run_sql <<SQL
.timer on

CREATE TABLE effects AS
SELECT row_number() OVER () AS row_id, *
FROM read_parquet('$SRC_DIR/EFFECTS.parquet');
ALTER TABLE effects ADD PRIMARY KEY (row_id);
CREATE INDEX ix_eff_uid_drug ON effects(UNIQUEID, DRUG);
CREATE INDEX ix_eff_drug_gm  ON effects(DRUG, GENE, MUTATION);
CREATE INDEX ix_eff_uid_cat  ON effects(UNIQUEID, DRUG, CATALOGUE_NAME);

CREATE TABLE predictions AS SELECT * FROM read_parquet('$SRC_DIR/PREDICTIONS.parquet');
ALTER TABLE predictions ADD PRIMARY KEY (UNIQUEID, DRUG, CATALOGUE_NAME, CATALOGUE_VERSION);
CREATE INDEX ix_pred_uid_drug ON predictions(UNIQUEID, DRUG);

SELECT 'effects' AS tbl, COUNT(*) AS n FROM effects
UNION ALL SELECT 'predictions', COUNT(*) FROM predictions;
SQL

stage_end "04_effects_predictions"
