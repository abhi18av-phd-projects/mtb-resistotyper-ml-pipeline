/*
 * Apply a chain of feature-engineering steps to a mart.
 *
 * `phase` is passed to the registry, which refuses a label-using step in the
 * pre-mart phase. That check is why this process exists: a misplacement becomes
 * a configuration error rather than a mart that looks fine and was built while
 * looking at the lineage it will be tested on.
 *
 * Instantiated twice -- once before the label boundary for the label-free layers,
 * once inside the fold for the label-using ones -- so both chains are the same
 * code and cannot drift apart.
 */
process FE_STEP {
    tag "${drug}/${fe_name}/${phase}${held_out == 'none' ? '' : '/held-out-' + held_out}"
    label 'mart'

    input:
    tuple val(drug), val(fe), val(fe_name), path(mart), path(meta), val(held_out)
    val phase
    val steps
    path annotations

    output:
    tuple val(drug), val(fe), val(fe_name), path('fe_mart.parquet'), path(meta),
          val(held_out), path('fe_mart.steps.json'), emit: mart

    script:
    def arm = held_out == 'none' ? '' : "--held-out-lineage ${held_out}"
    """
    cd ${params.project_root} && \\
    ${params.python} -m analysis.scripts.fe_steps.run_step \\
        --mart \$OLDPWD/${mart} \\
        --out \$OLDPWD/fe_mart.parquet \\
        --notes \$OLDPWD/fe_mart.steps.json \\
        --steps '${steps}' \\
        --phase ${phase} \\
        --annotations \$OLDPWD/${annotations} \\
        ${arm}
    """
}
