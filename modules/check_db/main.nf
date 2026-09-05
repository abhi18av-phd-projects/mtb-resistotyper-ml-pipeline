/*
 * Referential-integrity gate. Stage 11 of the shell build created the check
 * VIEWS but never failed on them, so a database could carry orphaned rows into
 * every mart downstream. Here a violation stops the run.
 *
 * This used to shell out to the `duckdb` CLI, which is not in the image: the
 * environment installs the duckdb PYTHON package, and pip does not ship the
 * binary. It also queried a table called `ri_summary` that does not exist --
 * stage 11 creates `v_ri_*` VIEWS -- so the gate would have failed every run
 * for the wrong reason even had the binary been there.
 *
 * The database is NOT re-emitted as an output. It used to be, so that the
 * checked file could be passed downstream, and on an S3 work dir that fails:
 * nf-nomad-s5cmd sweeps staged inputs after the task ("rm s3://.../inputs/*"),
 * so the file is gone before Nextflow collects the outputs and the process dies
 * with "Missing output file(s) cryptic.duckdb" after having exited 0. The
 * ordering it existed to enforce is now done in main.nf by joining on the
 * summary, which is a real output.
 *
 * The checker is invoked from the image at ${params.project_root}, NOT from the
 * pipeline's bin/ directory. Putting it in bin/ makes Nextflow stage that
 * directory for the task, and on Nextflow 26.04.3 the s5cmd work-dir provider
 * fails to enumerate it -- "session.getBinDirs() failed: No signature of method:
 * nextflow.Session.getBinDirs()" -- falls back to uploading the directory whole,
 * and the allocation then fails. nf-nomad reports that failure as success,
 * because it reads the ABSENT .exitcode as 0 ("reported Nomad alloc-state
 * failure but local .exitcode = 0; trusting the worker exit code"), so the head
 * sees a task that exited 0 and produced no outputs. None of that reaches
 * stdout. The analysis code is already baked into the image, so calling it by
 * path avoids the staging path entirely.
 *
 * Both duckdb bugs are fixed by check_integrity.py, which discovers the checks
 * from the catalogue and scopes blocking to the tables the marts actually read.
 * Orphans in the analysis path fail; orphans elsewhere are reported. On the
 * v3.4.0 build the difference is 45,062 plate-level and citizen-science orphans
 * that no mart reads: a gate that blocks on findings outside what it protects
 * is a gate people switch off, which is worse than a narrow one.
 */
process CHECK_DB {
    tag "ri-checks"
    label 'light'
    publishDir "${params.outdir}/database", mode: 'copy', pattern: '*.{tsv,json}'

    input:
    path db

    output:
    path 'ri_summary.tsv',  emit: summary
    path 'ri_summary.json', emit: json

    script:
    """
    bash ${params.project_root}/analysis/scripts/build_cryptic_db/11_ri_checks.sh \\
        "" ${db} > ri_run.log 2>&1 || true

    # NOT --strict. That flag fails on advisory orphans too, which is precisely
    # the behaviour the scoped gate exists to avoid: the 45,062 peripheral
    # orphans would block the run again, and the scoping would be decoration.
    python ${params.project_root}/analysis/scripts/build_cryptic_db/check_integrity.py \\
        --db ${db} --out ri_summary.tsv --json ri_summary.json
    """
}
