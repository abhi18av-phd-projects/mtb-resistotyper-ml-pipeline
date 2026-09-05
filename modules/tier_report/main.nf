/*
 * Collect the fold-internal evaluations into the per-drug operating range a
 * model card must carry.
 *
 * Only the lineage a fold's selection never saw contributes that fold's number.
 * The report also emits what the same folds would claim if every lineage were
 * counted, so the cost of selection leakage is reported rather than removed
 * silently.
 */
process TIER_REPORT {
    tag "${drug}/${fe_name}"
    label 'light'
    publishDir "${params.outdir}/tiers", mode: 'copy'

    input:
    tuple val(drug), val(fe_name), path(folds, stageAs: 'folds/*')

    output:
    tuple val(drug), path("tiers_${drug}_${fe_name}.json"), emit: tiers

    script:
    """
    ${params.python} ${params.project_root}/analysis/scripts/feature_mart/tier_report.py \\
        --drug ${drug} \\
        --feature-set ${params.h2o_feature_set} \\
        --folds ${folds} \\
        --out tiers_${drug}_${fe_name}.json
    """
}
