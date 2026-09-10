/*
 * The Binary Builder, one task per drug.
 *
 * catomatic is consumed as published: it is on PyPI, it exposes a command-line
 * interface, and its dependencies already include piezo. No wrapper is written,
 * which is the point -- the claim this stage supports is that a published,
 * notebook-bound analysis survives being decomposed into declared cluster jobs,
 * and rewriting it would forfeit that claim.
 *
 * Read support is filtered upstream in EXTRACT_COHORT rather than by catomatic's
 * own --frs. Both would filter; only the extraction can report what fraction of
 * the rows in scope carried an FRS value at all, and on this compendium that is
 * 10.19% -- the difference between a threshold that selects and one that
 * silently discards nearly everything.
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

    script:
    def name = "${params.catalogue_name}-${dataset_tag}"
    def version = "${params.catalogue_version}"
    // String param, and "false" is truthy in Groovy.
    def strict = params.catalogue_strict_unlock.toString().toBoolean() ? '--strict_unlock' : ''
    """
    # Exported here, not read from the ambient environment: a task runs in a
    # container that inherits nothing, so a manifest that reads these from
    # os.environ records "unknown" for all three unless the process puts them
    # there. The first run of this stage did exactly that.
    export MTB_CONTAINER_TAG='${params.catalogue_container_tag}'
    export MTB_GIT_REVISION='${workflow.commitId ?: workflow.scriptId}'
    export NXF_UUID='${workflow.sessionId}'

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
        --grammar ${params.catalogue_grammar} \\
        --values ${params.catalogue_values} \\
        --test ${params.catalogue_test} \\
        --background ${params.catalogue_background} \\
        --p ${params.catalogue_p} \\
        ${strict}

    catomatic binary \\
        --samples ${samples} \\
        --mutations ${mutations} \\
        --to_json \\
        --outfile catalogue_${drug}.json \\
        --test ${params.catalogue_test} \\
        --background ${params.catalogue_background} \\
        --p ${params.catalogue_p} \\
        ${strict}

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
        "strict_unlock": ${params.catalogue_strict_unlock},
        "grammar": "${params.catalogue_grammar}",
        "values": "${params.catalogue_values}",
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

    stub:
    """
    printf 'GENBANK_REFERENCE,CATALOGUE_NAME,CATALOGUE_VERSION,CATALOGUE_GRAMMAR,PREDICTION_VALUES,DRUG,MUTATION,PREDICTION,SOURCE,EVIDENCE,OTHER\n' > catalogue_${drug}.csv
    printf 'NC_000962.3,stub,1,GARC1,RUS,${drug},rpoB@S450L,R,{},{},{}\n' >> catalogue_${drug}.csv
    echo '{"stub":true}' > catalogue_${drug}.json
    echo '{"drug":"${drug}","stub":true}' > build_${drug}.json
    """
}
