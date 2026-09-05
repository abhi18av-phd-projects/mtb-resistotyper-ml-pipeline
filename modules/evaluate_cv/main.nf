/*
 * Random, lineage-leave-one-out and country-held-out on identical features and
 * folds. The random regime is computed only so the inflation it introduces can
 * be reported as a quantity rather than asserted as a concern.
 *
 * The output file is named for the lineage this fold's selection withheld, so
 * the held-out identity travels with the artefact and TIER_REPORT can take the
 * one honest number from each fold rather than trusting argument order.
 */
process EVALUATE_CV {
    tag "${drug}/held-out-${held_out}"
    label 'evaluate'

    input:
    tuple val(drug), val(fe), val(held_out), path(mart), path(meta),
          path(concordance), path(selection)

    output:
    tuple val(drug), val(fe), val(held_out), path("fold_${held_out}.json"), emit: metrics

    script:
    """
    cd ${params.project_root} && \\
    ${params.python_causal} -m analysis.scripts.feature_mart.evaluate_cv \\
        --mart \$OLDPWD/${mart} \\
        --concordance \$OLDPWD/${concordance} \\
        --drug ${drug} \\
        --out \$OLDPWD/cv/

    cd \$OLDPWD && cp cv/manifest.json fold_${held_out}.json
    """
}
