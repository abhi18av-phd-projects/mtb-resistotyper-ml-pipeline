/*
 * Identify the inputs by content, not by path.
 *
 * A release directory edited in place is otherwise indistinguishable from one
 * that has not been, and every downstream claim of reproducibility rests on
 * knowing exactly which bytes were read.
 */
process VERIFY_SOURCE {
    tag "${src.name}"
    label 'light'

    input:
    path src

    output:
    path src,                        emit: src
    path 'source_manifest.json',     emit: manifest

    script:
    """
    ${params.python} ${params.pipeline_dir}/bin/source_manifest.py \\
        --src ${src} \\
        --release '${params.release ?: ""}' \\
        --out source_manifest.json
    """
}
