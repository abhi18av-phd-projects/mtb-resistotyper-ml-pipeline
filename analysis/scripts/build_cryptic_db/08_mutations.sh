#!/usr/bin/env bash
# Stage 08 — gene-level mutation calls (LARGE).
#
# Loads: mutations (80,069,097 rows ⇒ ~2.5 GB on-disk after indexes)
#
# This is the table the resistotyping pipeline encodes as features. One row
# per (UNIQUEID, GENE, MUTATION) detected by gnomonicus, with quality fields
# (FRS, COVERAGE, IS_MINOR, etc.).
#
# Why a surrogate row_id PK:
#   The natural near-key (UNIQUEID, GENE, MUTATION, IS_MINOR) has ~469
#   duplicates out of 80M (each carrying distinct quality fields). Strict PK
#   would fail; surrogate PK keeps every row addressable, indexes restore
#   join performance.
#
# Indexes (load 4× longer than table itself):
#   - (UNIQUEID)        — for "all mutations in this sample" queries
#   - (GENE, MUTATION)  — for "all samples carrying mutation X" queries
#   - (UNIQUEID, GENE)  — for gene-level aggregation
#
# Expected wall time: ~4 min load + ~3 min indexes (~395 s total on Apple M-series).
#
# This stage benefits from running with as much RAM as possible; the parquet
# file is 1 GB and the in-memory representation during pivot operations can be
# several GB.

set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "$SCRIPT_DIR/_lib.sh"
resolve_args "$@"
stage_start "08_mutations"

require_source_file "MUTATIONS.parquet"

run_sql <<SQL
.timer on

CREATE TABLE mutations AS
SELECT row_number() OVER () AS row_id, *
FROM read_parquet('$SRC_DIR/MUTATIONS.parquet');
ALTER TABLE mutations ADD PRIMARY KEY (row_id);
CREATE INDEX ix_mut_uid       ON mutations(UNIQUEID);
CREATE INDEX ix_mut_gene_mut  ON mutations(GENE, MUTATION);
CREATE INDEX ix_mut_uid_gene  ON mutations(UNIQUEID, GENE);

SELECT COUNT(*) AS n FROM mutations;
SQL

stage_end "08_mutations"
