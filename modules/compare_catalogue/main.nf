/*
 * The acceptance test: built catalogue against a published reference.
 *
 * Reported, not adjudicated. The reference ships at every parameter setting
 * (`frs/<drug>/bg_0.05_p_0.9_FRS_{0.1..0.9}.csv`), so a comparison is available
 * per drug AND per read-support threshold rather than only at the end.
 */
process COMPARE_CATALOGUE {
    tag "${drug}/frs${frs}"
    label 'catalogue'

    input:
    tuple val(drug), val(dataset_tag), val(frs), path(built), path(cat_json), path(build_json)
    path reference

    output:
    tuple val(drug), val(dataset_tag), val(frs),
          path("comparison_${drug}.json"), emit: comparison

    stub:
    """
    echo '{"drug":"${drug}","stub":true}' > comparison_${drug}.json
    """

    script:
    """
    cd ${params.project_root} && \\
    ${params.python_catalogue} -m analysis.scripts.catalogue.compare_catalogue \\
        --built \$OLDPWD/${built} \\
        --reference \$OLDPWD/${reference} \\
        --drug ${drug} \\
        --out \$OLDPWD/

    cd \$OLDPWD
    """
}
