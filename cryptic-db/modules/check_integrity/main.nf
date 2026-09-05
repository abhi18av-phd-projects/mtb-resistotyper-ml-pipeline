/*
 * Stage 11 as a GATE, scoped to what it protects.
 *
 * The shell build created the referential-integrity check views but never
 * failed on them, so a database with orphaned rows could go on to produce marts
 * that look perfectly well formed.
 *
 * Scope matters as much as the check. On the v3.4.0 build, blocking on ANY
 * orphan means failing on 43,872 in ukmyc_growth_to_plates and 1,190 in
 * bashthebug_to_plates -- plate readings and citizen-science classifications
 * that no feature mart reads. A gate that fails on findings outside what it
 * protects is a gate people switch off. So orphans in the analysis path fail
 * the run and orphans elsewhere are reported; --strict collapses the two.
 */
process CHECK_INTEGRITY {
    tag "ri-checks"
    label 'load'

    input:
    path db

    output:
    path db,                emit: db
    path 'ri_summary.tsv',  emit: summary
    path 'ri_summary.json', emit: report

    script:
    def strict = params.strict_integrity ? '--strict' : ''
    """
    bash ${params.scripts_dir}/11_ri_checks.sh "${params.cryptic_src}" ${db} > ri.log 2>&1

    ${params.python} ${params.pipeline_dir}/bin/check_integrity.py \\
        --db ${db} \\
        --out ri_summary.tsv \\
        --json ri_summary.json \\
        ${strict}
    """
}
