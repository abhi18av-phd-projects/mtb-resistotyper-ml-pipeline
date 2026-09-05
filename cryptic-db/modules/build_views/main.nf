/* Stage 10 — convenience views. Depends on the assembled whole. */
process BUILD_VIEWS {
    tag "views"
    label 'load'

    input:
    path db

    output:
    path db, emit: db

    script:
    """
    bash ${params.scripts_dir}/10_views.sh "${params.cryptic_src}" ${db} > views.log 2>&1
    """
}
