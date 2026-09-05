#!/usr/bin/env bash
# Stage 09 — VCF-level variant calls (LARGE).
#
# Loads: variants (83,581,726 rows ⇒ ~9 GB on-disk after indexes)
#
# Lower-level than mutations: one row per VCF record (per-base SNP / indel
# / null call). The mutations table is gnomonicus's interpretation of these
# variants in terms of named gene mutations. Variants is what you want for
# raw genotype features, mutations is what you want for catalogue-aware
# features.
#
# Same surrogate-row_id pattern as mutations (~14k near-dups on the
# (UNIQUEID, GENE, VARIANT, MINOR_VARIANT) composite).
#
# Indexes:
#   - (UNIQUEID)        — for "all variants in this sample"
#   - (GENE, VARIANT)   — for "all samples carrying variant X"
#
# Expected wall time: ~2 min load + ~2 min indexes (~227 s total).

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "09_variants"

require_source_file "VARIANTS.parquet"

run_sql <<SQL
.timer on

CREATE TABLE variants AS
SELECT row_number() OVER () AS row_id, *
FROM read_parquet('$SRC_DIR/VARIANTS.parquet');
ALTER TABLE variants ADD PRIMARY KEY (row_id);
CREATE INDEX ix_var_uid       ON variants(UNIQUEID);
CREATE INDEX ix_var_gene_var  ON variants(GENE, VARIANT);

SELECT COUNT(*) AS n FROM variants;
SQL

stage_end "09_variants"
