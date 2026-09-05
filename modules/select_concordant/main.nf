/*
 * Concordance selection: EBM importance, causal-forest CATE with refutation,
 * and an FWL lineage-residual partial correlation. Two of three votes admit a
 * mutation.
 *
 * held_out names the lineage withheld from selection. Anything other than
 * 'none' is the EVALUATION arm and the selection must not see that lineage;
 * 'none' is the DEPLOYMENT refit, where using all data is correct. The two are
 * kept in one process so they cannot drift apart, and distinguished by a value
 * that ends up in the output path.
 */
process SELECT_CONCORDANT {
    tag "${drug}/${fe.name}/held-out-${held_out}"
    label 'select'

    input:
    tuple val(drug), val(fe), val(fe_name), path(mart), path(meta), val(held_out)
    path db

    output:
    tuple val(drug), val(fe), val(held_out), path(mart), path(meta),
          path('concordance.parquet'), path('selection.json'), emit: selected

    script:
    def arm = held_out == 'none' ? '--refit-on-all' : "--held-out-lineage ${held_out}"
    // Null means causal.py's own default (300). The test profile shrinks it so a
    // validation run costs minutes rather than the half hour a full pool takes:
    // each candidate carries a CATE fit and a refutation pass, so wall time is
    // close to linear in the pool.
    // causal.py's --db defaults to a hardcoded slim-database path that exists on
    // exactly one laptop. Its only use is the catalogue benchmark, which compares
    // the concordant set against the WHO catalogue -- so pass the database the run
    // is actually using rather than letting the default decide.
    def pool = params.candidate_pool_size ? "--pool-size ${params.candidate_pool_size}" : ''
    """
    cd ${params.project_root} && \\
    ${params.python_causal} -m analysis.scripts.feature_mart.causal \\
        --mart \$OLDPWD/${mart} \\
        --drug ${drug} \\
        --level ${params.candidate_level} \
        --db \$OLDPWD/${db} \
        ${pool} \\
        --out \$OLDPWD/ \\
        ${arm}

    cd \$OLDPWD
    # Name the primary output rather than globbing for it. Passing --db turned on
    # the catalogue benchmark, which writes a SECOND causal_concordance_*.parquet
    # (…_vs_catalogue). The glob then matched two files and mv read the last
    # argument as a destination directory: "target 'concordance.parquet': No such
    # file or directory". The benchmark stays in the work directory; it is a
    # comparison against the WHO catalogue, not the selection this process emits.
    mv causal_concordance_${drug}_${params.candidate_level}_level.parquet concordance.parquet
    mv *manifest*.json selection.json 2>/dev/null || echo '{}' > selection.json
    """
}
