/*
 * The Binary Builder, one task per drug.
 *
 * catomatic is consumed as published: it is on PyPI, it exposes a command-line
 * interface, and its dependencies already include piezo. No wrapper is written,
 * which is the point -- the claim this stage supports is that a published,
 * notebook-bound analysis survives being decomposed into declared cluster jobs,
 * and rewriting it would forfeit that claim.
 *
 * Read support is filtered upstream in EXTRACT_COHORT. catomatic has no FRS
 * flag; the published sweep varies FRS at fixed background and p, so the
 * threshold is an input to the extraction, not to the builder.
 */
process RUN_CATOMATIC {
    tag "${drug}/frs${frs}"
    label 'catalogue'

    input:
    tuple val(drug), val(dataset_tag), val(frs), path(samples), path(mutations),
          path(cohort), path(wildcards)

    output:
    tuple val(drug), val(dataset_tag), val(frs),
          path("catalogue_${drug}.csv"),
          path("catalogue_${drug}.json"),
          path("build_${drug}.json"), emit: catalogue

    stub:
    """
    printf 'GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER\n' > catalogue_${drug}.csv
    printf 'NC_000962.3,stub,1,GARC1,RUS,${drug},rpoB@S450L,R,{},{},{}\n' >> catalogue_${drug}.csv
    echo '{"stub":true}' > catalogue_${drug}.json
    echo '{"drug":"${drug}","stub":true}' > build_${drug}.json
    """

    script:
    def name = "${params.catalogue_name}-${dataset_tag}"
    def version = "${params.catalogue_version}"
    """
    catomatic binary \\
        --samples ${samples} \\
        --mutations ${mutations} \\
        --to_piezo \\
        --outfile catalogue_${drug}.csv \\
        --genbank_ref '${params.catalogue_genbank_ref}' \\
        --catalogue_name '${name}' \\
        --version '${version}' \\
        --drug ${drug} \\
        --wildcards ${wildcards} \\
        --test ${params.catalogue_test} \\
        --background ${params.catalogue_background} \\
        --p ${params.catalogue_p} \\
        --tails ${params.catalogue_tails}

    catomatic binary \\
        --samples ${samples} \\
        --mutations ${mutations} \\
        --to_json \\
        --outfile catalogue_${drug}.json \\
        --test ${params.catalogue_test} \\
        --background ${params.catalogue_background} \\
        --p ${params.catalogue_p} \\
        --tails ${params.catalogue_tails}

    # The build manifest is what makes the catalogue interpretable later: the
    # compendium release and cohort come from the extraction's own report, and
    # the grading parameters, image tag and revision from this run. A catalogue
    # that cannot name these is not a release, it is a file.
    ${params.python_catalogue} - <<'PY'
import json, os
from pathlib import Path
cohort = json.loads(Path("${cohort}").read_text())
Path("build_${drug}.json").write_text(json.dumps({
    "drug": "${drug}",
    "catalogue_name": "${name}",
    "catalogue_version": "${version}",
    "genbank_ref": "${params.catalogue_genbank_ref}",
    "grading": {
        "test": "${params.catalogue_test}",
        "background": ${params.catalogue_background},
        "p": ${params.catalogue_p},
        "tails": "${params.catalogue_tails}",
        "frs_threshold": cohort.get("frs_threshold"),
    },
    "cohort": cohort,
    "provenance": {
        "container_tag": os.environ.get("MTB_CONTAINER_TAG", "unknown"),
        "git_revision": os.environ.get("MTB_GIT_REVISION", "unknown"),
        "nextflow_run": os.environ.get("NXF_UUID", "unknown"),
    },
}, indent=2) + "\\n")
PY
    """
}
