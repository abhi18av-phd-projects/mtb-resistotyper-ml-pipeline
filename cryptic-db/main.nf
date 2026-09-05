#!/usr/bin/env nextflow

/*
 * cryptic-db — CRyPTIC release tables into a reproducible DuckDB.
 *
 * Separate from the training pipeline on purpose. The database changes when
 * CRyPTIC cuts a release; the feature engineering changes every afternoon.
 * Rebuilding a 3.5 GB database to try a different top-N is the wrong unit of
 * work, so this emits a versioned, checksummed artefact that the training
 * pipeline consumes with --db.
 *
 * Stages 01-09 of the shell build each read only the source parquet they load
 * and the tables they themselves create -- verified, not assumed -- so they run
 * as independent shards here and are attached together afterwards. Stages 10-12
 * genuinely depend on the assembled whole and stay sequential.
 */

include { VERIFY_SOURCE   } from './modules/verify_source'
include { LOAD_SHARD      } from './modules/load_shard'
include { ASSEMBLE_DB     } from './modules/assemble_db'
include { BUILD_VIEWS     } from './modules/build_views'
include { CHECK_INTEGRITY } from './modules/check_integrity'
include { WRITE_METADATA  } from './modules/write_metadata'
include { EXPORT_PARQUET  } from './modules/export_parquet'
include { MAKE_SLIM       } from './modules/make_slim'

/*
 * The independent load stages. Order is irrelevant -- each reads only the source
 * parquet it loads and the tables it creates, so they run concurrently.
 */
def shards() {
    ['01_reference', '02_sample_anchors', '03_one_to_one', '04_effects_predictions',
     '05_phenotypes', '06_bashthebug', '07_ukmyc_growth', '08_mutations']
}

/*
 * Stages whose source table is not in every release. v3.4.0 ships no
 * VARIANTS.parquet, so running 09_variants against it fails on a missing file
 * rather than on anything interesting. Included only when the input is there.
 */
def optionalShards(Path src) {
    def present = []
    if( file("${src}/VARIANTS.parquet").exists() )
        present << '09_variants'
    else
        log.warn "VARIANTS.parquet absent from ${src} — skipping stage 09_variants"
    return present
}

workflow {
    main:

    log.info """\
      cryptic-db
      ==========
      source  : ${params.cryptic_src}
      release : ${params.release ?: '(read from the source directory name)'}
      outdir  : ${params.outdir}
      shards  : ${shards().size()} required, + any optional
      """.stripIndent()

    if( !params.cryptic_src )
        error "Provide --cryptic_src <cryptic-tables-vX.Y.Z/>"

    ch_src = channel.fromPath(params.cryptic_src, type: 'dir', checkIfExists: true)

    /* ---- provenance first -----------------------------------------------
     * Checksum every input before reading any of it. A database is only
     * reproducible if the inputs it was built from are identified by content
     * rather than by path, and a release directory that has been edited in
     * place is otherwise indistinguishable from one that has not.
     */
    ch_verified = VERIFY_SOURCE(ch_src)

    /* ---- independent loads ---------------------------------------------- */
    ch_shards = LOAD_SHARD(
        ch_verified.src.flatMap { src ->
            (shards() + optionalShards(src)).collect { stage -> tuple(stage, src) } })

    /* ---- assemble, finish, gate ----------------------------------------- */
    ch_db     = ASSEMBLE_DB(ch_shards.shard.collect())
    ch_views  = BUILD_VIEWS(ch_db)
    ch_meta   = WRITE_METADATA(ch_views, ch_verified.manifest)

    // Referential integrity is a GATE. A database that fails its own checks is
    // worse than no database, because every mart built on it looks fine.
    ch_final  = CHECK_INTEGRITY(ch_meta)

    ch_parquet = params.export_parquet ? EXPORT_PARQUET(ch_final.db) : channel.empty()

    /* ---- optional reduced copy -------------------------------------------
     * Off by default. Slim is a speed hack, and not a free one: it drops
     * EVIDENCE and PREDICTION_VALUES from `effects`, which is what RFUS/SOLO
     * grading, the EVIDENCE JSON and real WHO tier grading all need. Feature
     * engineering belongs on the full database; reach for slim when the wait is
     * what you are optimising.
     */
    ch_slim = params.make_slim ? MAKE_SLIM(ch_final.db).db : channel.empty()

    publish:
    database = ch_final.db
    checks   = ch_final.summary
    manifest = ch_verified.manifest
    parquet  = ch_parquet
    slim     = ch_slim
}

output {
    database { path '.' }
    checks   { path 'checks' }
    manifest { path 'provenance' }
    parquet  { path 'parquet' }
    slim     { path 'slim' }
}
