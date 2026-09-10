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

    // Named by drug AND fold. Five instances of this process each emitted
    // 'concordance.parquet' and 'selection.json', so collecting them into one task
    // collided, and the per-file staging directories added to dodge that
    // collision then separated each concordance from its manifest -- which the
    // deposit pairs BY DIRECTORY. Naming by fold fixed one drug and broke two:
    // CREATE_DUCKDB collects every drug's selections, so RIF and INH both emitted
    // `concordance_lineage1.parquet` and the task refused the collision. The
    // single-drug campaign never reached it; a two-drug -stub-run does in
    // seconds. Drug and fold together make the files unique at the source, keep
    // a pair together, and let a reader recover both from the filename alone.
    output:
    tuple val(drug), val(fe), val(held_out), path(mart), path(meta),
          path("concordance_${drug}_${held_out}.parquet"),
          path("selection_${drug}_${held_out}.json"),
          emit: selected

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
    mv causal_concordance_${drug}_${params.candidate_level}_level.parquet concordance_${drug}_${held_out}.parquet
    mv *manifest*.json selection_${drug}_${held_out}.json 2>/dev/null || echo '{}' > selection_${drug}_${held_out}.json
    """

    stub:
    """
    touch concordance_${drug}_${held_out}.parquet
    echo '{"drug":"${drug}","held_out":"${held_out}","stub":true}' > selection_${drug}_${held_out}.json
    """
}
