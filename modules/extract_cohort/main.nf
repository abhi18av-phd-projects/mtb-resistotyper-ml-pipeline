/*
 * catomatic's two input tables, for one drug, from the compendium.
 *
 * The cohort is a release TAG (`wgs_samples.dataset`), not a sample list. That
 * is what makes the reproduction checkable rather than asserted, and it is why
 * this process fails loudly on a database that does not carry the column: a
 * silent fallback to the whole compendium yields a catalogue that looks like
 * the paper's and answers a different question.
 */
process EXTRACT_COHORT {
    tag "${drug}/${dataset_tag}/frs${frs}"
    label 'catalogue'

    input:
    tuple val(drug), val(genes), val(dataset_tag), val(frs), path(db)

    output:
    tuple val(drug), val(dataset_tag), val(frs),
          path("samples_${drug}.csv"),
          path("mutations_${drug}.csv"),
          path("cohort_${drug}.json"),
          path("wildcards_${drug}.json"), emit: cohort

    script:
    def frs_arg = frs == null || frs == 'none' ? '' : "--frs ${frs}"
    def gene_arg = genes ? "--genes ${genes}" : ''
    """
    cd ${params.project_root} && \\
    ${params.python_catalogue} -m analysis.scripts.catalogue.extract_cohort \\
        --db \$OLDPWD/${db} \\
        --drug ${drug} \\
        --dataset-tag ${dataset_tag} \\
        --phenotype-source ${params.catalogue_phenotype_source} \\
        ${frs_arg} \\
        ${gene_arg} \\
        --out \$OLDPWD/

    cd \$OLDPWD
    """

    stub:
    """
    printf 'UNIQUEID,PHENOTYPE\nstub.1,R\nstub.2,S\n' > samples_${drug}.csv
    printf 'UNIQUEID,MUTATION\nstub.1,rpoB@S450L\n'    > mutations_${drug}.csv
    echo '{"drug":"${drug}","dataset_tag":"${dataset_tag}","stub":true}' > cohort_${drug}.json
    echo '{"rpoB@*=":{"pred":"S"}}' > wildcards_${drug}.json
    """
}
