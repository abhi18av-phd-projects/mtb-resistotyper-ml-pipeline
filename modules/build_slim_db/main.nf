process BUILD_SLIM_DB {
    tag "cryptic-slim"
    label 'io_heavy'
    publishDir "${params.outdir}/database", mode: 'copy'

    input:
    path cryptic_src

    output:
    tuple path('cryptic-slim.duckdb'), path('build.log'), emit: db

    script:
    """
    bash ${params.project_root}/analysis/scripts/build_cryptic_db/run_all.sh \\
        ${cryptic_src} cryptic-slim.duckdb > build.log 2>&1
    """
}
