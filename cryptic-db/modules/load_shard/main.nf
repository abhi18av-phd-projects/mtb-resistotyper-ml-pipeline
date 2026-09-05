/*
 * One shell stage, one shard database.
 *
 * Stages 01-09 each read only the source parquet they load and the tables they
 * create, so they are independent and run concurrently. Each writes its OWN
 * DuckDB file; ASSEMBLE_DB attaches them afterwards. The shell scripts already
 * take (src, db) positionally, so they are reused unmodified -- the pipeline
 * changes where they write, not what they do.
 */
process LOAD_SHARD {
    tag "${stage}"
    label { stage in ['08_mutations', '09_variants'] ? 'load_large' : 'load' }

    input:
    tuple val(stage), path(src)

    output:
    tuple val(stage), path("${stage}.duckdb"), path("${stage}.log"), emit: shard

    script:
    """
    bash ${params.scripts_dir}/${stage}.sh ${src} ${stage}.duckdb > ${stage}.log 2>&1
    """
}
