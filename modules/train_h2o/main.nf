/*
 * H2O AutoML -> base models -> stacked ensembles.
 *
 * The JVM is booted inside the task and torn down with it. A shared long-lived
 * H2O server would be mutable state across tasks: one run's memory pressure
 * becomes another run's failure, a crash takes out unrelated work, and -resume
 * stops meaning anything. Per-task keeps the heap bounded by the allocation
 * that was actually reserved.
 *
 * h2o_mem is derived from the task's own memory rather than set independently,
 * so the JVM cannot overcommit what the scheduler promised it.
 */
process TRAIN_H2O {
    tag "${drug}/${fe.name}"
    label 'jvm'

    input:
    tuple val(drug), val(fe), val(held_out), path(mart), path(meta),
          path(concordance), path(selection)

    output:
    tuple val(drug), val(fe), path('h2o/**'), path('h2o_manifest.json'), emit: model

    script:
    def heap = "${(task.memory.toGiga() * 0.8) as int}G"
    """
    export JAVA_HOME=${params.java_home}
    export PATH=\$JAVA_HOME/bin:\$PATH

    # Tracking is opt-in and non-fatal. Unset mlflow_uri and the training script
    # takes a no-op handle: the run still writes manifest.json, which stays the
    # source of truth. A tracking server that is down must not cost an hour of
    # cluster time.
    export MLFLOW_TRACKING_URI='${params.mlflow_uri ?: ""}'
    export MLFLOW_EXPERIMENT='${params.mlflow_experiment}'
    # Provenance the metric cannot be reproduced without, carried as tags so a
    # number in the UI can be traced back to the run and image that made it.
    # h2o_container_tag, not container_tag: this task runs mtb-h2o, and
    # params.container_tag is the mtb-causal image. Stamping the wrong tag
    # made every logged run name an image it had not run in.
    export MTB_CONTAINER_TAG='${params.h2o_container_tag}'
    export MTB_GIT_REVISION='${workflow.commitId ?: workflow.scriptId}'
    export NXF_UUID='${workflow.sessionId}'

    cd ${params.project_root} && \\
    ${params.python_h2o} -m analysis.scripts.feature_mart.${params.h2o_entrypoint} \\
        --mart \$OLDPWD/${mart} \\
        --concordance \$OLDPWD/${concordance} \\
        --feature-set ${params.h2o_feature_set} \\
        --max-runtime-secs ${params.h2o_runtime_secs} \\
        --drug ${drug} \\
        --sort-metric ${params.h2o_sort_metric} \\
        ${params.h2o_balance_classes ? '--balance-classes' : ''} \\
        --h2o-mem ${heap} \\
        --out \$OLDPWD/h2o/

    cd \$OLDPWD
    find h2o -name 'manifest.json' -print -quit | xargs -I{} cp {} h2o_manifest.json
    [ -f h2o_manifest.json ] || echo '{}' > h2o_manifest.json
    """
}
