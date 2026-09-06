// Assemble the publishable DuckDB: everything after feature engineering and
// before training, in one citable file.
//
// Off by default and read by nothing downstream. The pipeline does not need it;
// a reader who wants the data without the cluster does. Enable with
// `--create_duckdb <name>` (or `create_duckdb:` in an experiment file).

process CREATE_DUCKDB {
    tag "${params.experiment}"
    label 'mart'
    publishDir "${params.outdir}/publish", mode: 'copy'

    input:
    path marts
    path metas
    path selections

    output:
    path "${params.create_duckdb}", emit: db

    when:
    params.create_duckdb

    script:
    def sel = selections instanceof List ? selections.join(' ') : "${selections}"
    def mart = marts instanceof List ? marts.join(' ') : "${marts}"
    """
    cd ${params.project_root} && \\
    ${params.python} -m analysis.scripts.feature_mart.create_duckdb \\
        --marts ${mart.split(' ').collect { "\$OLDPWD/${it}" }.join(' ')} \\
        --selections ${sel ? sel.split(' ').collect { "\$OLDPWD/${it}" }.join(' ') : ''} \\
        --arm ${params.experiment} \\
        --out \$OLDPWD/${params.create_duckdb}
    """
}
