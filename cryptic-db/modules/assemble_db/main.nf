/*
 * Attach every shard and copy its tables into one database.
 *
 * DuckDB cannot have several writers on one file, which is why the loads
 * sharded in the first place. ATTACH is read-only here and the copy is a single
 * writer, so the concurrency stays where it is safe.
 */
process ASSEMBLE_DB {
    tag "assemble"
    label 'load_large'

    input:
    path shards, stageAs: 'shards/*'

    output:
    path 'cryptic.duckdb', emit: db

    script:
    """
    ${params.python} ${params.pipeline_dir}/bin/assemble_db.py \\
        --shards shards/*.duckdb \\
        --out cryptic.duckdb
    """
}
