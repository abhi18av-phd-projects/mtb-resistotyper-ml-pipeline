/*
 * The feature mart: the reproducibility unit.
 *
 * Every FE choice arrives as a parameter and is recorded in the provenance
 * sidecar, so a mart states how it was made and two marts are diffable. The
 * fe_hash in the output path is derived from the config, so a swept sweep
 * lands in distinct directories rather than overwriting itself.
 */
process BUILD_MART {
    tag "${drug}/${fe.name}"
    label 'mart'

    input:
    tuple val(drug), val(fe), path(db), path(cohort)

    output:
    tuple val(drug), val(fe), val("${fe.name}"),
          path("feature_mart_${drug}.parquet"),
          path("feature_mart_${drug}.metadata.json"), emit: mart

    script:
    """
    cd ${params.project_root} && \\
    ${params.python} -m analysis.scripts.feature_mart.build_mart \\
        --drug ${drug} \\
        --db \$OLDPWD/${db} \\
        --out \$OLDPWD/ \\
        --mutation-selector ${fe.mutation_selector} \\
        --top-n ${fe.top_n} \\
        --max-carrier-frac ${fe.max_carrier_frac} \\
        --min-carrier-count ${fe.min_carrier_count} \\
        --fold-strategy ${fe.fold_strategy} \\
        --n-folds ${fe.n_folds} \\
        ${fe.gene_level_flags ? '--gene-level-flags' : '--no-gene-level-flags'} \\
        --fe-name ${fe.name}

    cd \$OLDPWD
    mv feature_mart_${drug}_*.parquet      feature_mart_${drug}.parquet
    mv feature_mart_${drug}_*.metadata.json feature_mart_${drug}.metadata.json
    """
}
