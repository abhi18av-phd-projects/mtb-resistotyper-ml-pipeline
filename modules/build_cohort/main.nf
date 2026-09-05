process BUILD_COHORT {
    tag "cohort"
    label 'light'
    publishDir "${params.outdir}/cohort", mode: 'copy'

    input:
    path db

    output:
    // The DIRECTORY, not its contents. cohort.py writes four files, and
    // `cohort/*` expands to four separate path elements, so the tuple reaching
    // BUILD_MART carried seven elements against a four-element declaration.
    // BUILD_MART re-derives the cohort from the database anyway -- this edge
    // exists to order the two, not to hand over data.
    path 'cohort', emit: cohort

    script:
    """
    cd ${params.project_root} && \\
    ${params.python} -m analysis.scripts.feature_mart.cohort \\
        --db \$OLDPWD/${db} \\
        --out \$OLDPWD/cohort/
    """
}
