// Assemble the publishable DuckDB: everything after feature engineering and
// before training, in one citable file.
//
// Off by default and read by nothing downstream. The pipeline does not need it;
// a reader who wants the data without the cluster does. Enable with
// `--create_duckdb <name>` (or `create_duckdb:` in an experiment file).

process CREATE_DUCKDB {
    tag "${params.experiment}"
    label 'mart'
    // A deposit stage must not be able to kill a campaign. Nothing downstream
    // reads this file, so its failure costs a deposit and nothing else -- but
    // an un-ignored failure propagates to the head, the head exits non-zero,
    // Nomad restarts it, and with no cross-run cache the ENTIRE pipeline re-runs.
    // That happened: six full re-runs across two allocations, ~13 hours of a
    // shared single-node cluster, because an optional artefact could not be
    // written.
    errorStrategy 'ignore'
    publishDir "${params.outdir}/publish", mode: 'copy'

    input:
    path marts
    path metas
    path selections, stageAs: 'conc*/*'
    path manifests,  stageAs: 'man*/*'

    output:
    path "${params.create_duckdb}", emit: db

    when:
    params.create_duckdb

    script:
    def sel = selections instanceof List ? selections.join(' ') : "${selections}"
    def man = manifests instanceof List ? manifests.join(' ') : "${manifests}"
    def mart = marts instanceof List ? marts.join(' ') : "${marts}"
    """
    cd ${params.project_root} && \\
    ${params.python} -m analysis.scripts.feature_mart.create_duckdb \\
        --marts ${mart.split(' ').collect { "\$OLDPWD/${it}" }.join(' ')} \\
        --selections ${sel ? sel.split(' ').collect { "\$OLDPWD/${it}" }.join(' ') : ''} \\
        --manifests ${man ? man.split(' ').collect { "\$OLDPWD/${it}" }.join(' ') : ''} \\
        --arm ${params.experiment} \\
        --out \$OLDPWD/${params.create_duckdb}
    """
}
