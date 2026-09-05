/*
 * Mirror every table to parquet.
 *
 * The training pipeline reads DuckDB, but a parquet mirror is what makes the
 * release readable without DuckDB at all -- by another language, another tool,
 * or a reviewer who should not have to install anything to look at the data.
 */
process EXPORT_PARQUET {
    tag "parquet"
    label 'load_large'

    input:
    path db

    output:
    path 'parquet/*.parquet', emit: parquet

    script:
    """
    mkdir -p parquet
    ${params.python} - <<'PY'
import duckdb
con = duckdb.connect("${db}", read_only=True)
tables = [r[0] for r in con.execute(
    "SELECT table_name FROM information_schema.tables "
    "WHERE table_schema='main' AND table_type='BASE TABLE'").fetchall()]
for t in tables:
    con.execute(f"COPY (SELECT * FROM \\"{t}\\") TO 'parquet/{t}.parquet' (FORMAT PARQUET)")
print(f"exported {len(tables)} table(s)")
PY
    """
}
