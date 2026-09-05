/*
 * OPTIONAL — a reduced copy of the database, for iterating quickly.
 *
 * Slim was only ever a speed hack, and it is not free: it drops EVIDENCE and
 * PREDICTION_VALUES from `effects`, and the plate-level tables entirely. Three
 * catalogued feature-engineering techniques depend on exactly those columns --
 * RFUS/SOLO grading, the EVIDENCE JSON, and real WHO tier grading, which
 * build_mart currently substitutes with per-mutation R/S/U/F counts because the
 * slim table cannot support it.
 *
 * So: build the full database, do the science on it, and reach for slim only
 * when the wait is what you are optimising. The reduced copy records what it
 * dropped so a mart built from it cannot silently claim capabilities it lacks.
 */
process MAKE_SLIM {
    tag "slim"
    label 'load_large'
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path db

    output:
    path 'cryptic-slim.duckdb', emit: db
    path 'slim_manifest.json',  emit: manifest

    script:
    """
    ${params.python} ${params.pipeline_dir}/bin/make_slim.py \\
        --full ${db} --out cryptic-slim.duckdb --manifest slim_manifest.json
    """
}
