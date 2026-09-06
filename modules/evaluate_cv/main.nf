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
    # Without these the tracking in evaluate_cv.py degrades to a silent no-op:
    # it keys off MLFLOW_TRACKING_URI, so an unset variable means the honest
    # evaluation is computed, written to the manifest, and never reaches the
    # dashboard. Tracking stays non-fatal; the manifest remains the truth.
    export MLFLOW_TRACKING_URI='${params.mlflow_uri ?: ""}'
    export MLFLOW_EXPERIMENT='${params.mlflow_experiment}'
    # container_tag is correct here: this task runs mtb-causal.
    export MTB_CONTAINER_TAG='${params.container_tag}'
    export MTB_GIT_REVISION='${workflow.commitId ?: workflow.scriptId}'
    export NXF_UUID='${workflow.sessionId}'

    cd ${params.project_root} && \\
    ${params.python_causal} -m analysis.scripts.feature_mart.evaluate_cv \\
        --mart \$OLDPWD/${mart} \\
        --concordance \$OLDPWD/${concordance} \\
        --drug ${drug} \\
        --held-out ${held_out} \\
        --out \$OLDPWD/cv/

    cd \$OLDPWD && cp cv/manifest.json fold_${held_out}.json
    """
}
