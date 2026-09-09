#!/usr/bin/env nextflow

/*
 * mtb-resistotyper-ml — feature engineering and model training
 *
 * Replaces the sequential shell orchestration (analysis/scripts/build_cryptic_db/
 * run_all.sh) and the by-hand invocation of the feature-mart and training scripts.
 *
 * The pipeline exists so that FEATURE ENGINEERING IS A PARAMETER. The mart, not
 * the model, is the reproducibility unit: change an FE choice, get a new mart with
 * a new provenance sidecar, and everything downstream re-runs while everything
 * untouched stays cached. Two FE configurations can then be compared on identical
 * folds rather than on recollected memory of how each was built.
 */

include { BUILD_SLIM_DB     } from './modules/build_slim_db'
include { CHECK_DB          } from './modules/check_db'
include { BUILD_COHORT      } from './modules/build_cohort'
include { BUILD_MART        } from './modules/build_mart'
include { FE_STEP           } from './modules/fe_step'
include { FE_STEP as FE_RANK } from './modules/fe_step'
include { FE_STEP as FE_PRE  } from './modules/fe_step'
include { FE_STEP as FE_FOLD } from './modules/fe_step'
include { SELECT_CONCORDANT } from './modules/select_concordant'
include { SELECT_CONCORDANT as SELECT_FULL } from './modules/select_concordant'
include { TRAIN_H2O         } from './modules/train_h2o'
include { TRAIN_H2O as TRAIN_H2O_FOLD } from './modules/train_h2o'
include { EVALUATE_CV       } from './modules/evaluate_cv'
include { TIER_REPORT       } from './modules/tier_report'
include { CREATE_DUCKDB     } from './modules/create_duckdb'
include { UNPACK_DUCKDB     } from './modules/unpack_duckdb'
include { EXTRACT_COHORT    } from './modules/extract_cohort'
include { RUN_CATOMATIC     } from './modules/run_catomatic'
include { COMPARE_CATALOGUE } from './modules/compare_catalogue'


/*
 * The FE configurations to build marts for.
 *
 * One configuration by default, taken from the flat fe_* params so every knob
 * is an ordinary command-line flag. --fe_configs points at a CSV or JSON of
 * several, each becoming its own mart against identical downstream folds; keys
 * absent from a row fall back to the fe_* default, so a sweep file only has to
 * state what it varies.
 */
def feDefaults() {
    [
        name              : params.fe_name,
        mutation_selector : params.fe_mutation_selector,
        top_n             : params.fe_top_n,
        max_carrier_frac  : params.fe_max_carrier_frac,
        min_carrier_count : params.fe_min_carrier_count,
        gene_level_flags  : params.fe_gene_level_flags,
        fold_strategy     : params.fe_fold_strategy,
        n_folds           : params.fe_n_folds,
    ]
}

def feConfigs() {
    if( !params.fe_configs )
        return [feDefaults()]

    def f = file(params.fe_configs, checkIfExists: true)
    def rows = f.name.endsWith('.json')
        ? new groovy.json.JsonSlurper().parse(f)
        : f.splitCsv(header: true)
    if( !rows )
        error "--fe_configs ${params.fe_configs} contained no configurations"

    return rows.collect { row ->
        def cfg = feDefaults() + row.collectEntries { k, v -> [k, coerce(k, v)] }
        if( !cfg.name ) error "every fe_configs entry needs a unique 'name'"
        cfg
    }
}

/* CSV gives every field back as a String; the process needs real types. */
def coerce(String key, value) {
    if( !(value instanceof String) )
        return value
    if( key in ['top_n', 'min_carrier_count', 'n_folds'] )
        return value as Integer
    if( key == 'max_carrier_frac' )
        return value as Double
    if( key == 'gene_level_flags' )
        return value.toBoolean()
    return value
}


workflow {
    main:

    log.info """\
      mtb-resistotyper-ml  ·  FE + training
      =====================================
      cryptic source : ${params.cryptic_src ?: '(none — using --db)'}
      database       : ${params.db ?: '(built by the pipeline)'}
      drugs          : ${params.drugs}
      held-out folds : ${params.lineages}
      experiment     : ${params.experiment}
      tracking       : ${params.mlflow_uri ?: '(off — results still written to results/)'}
      FE config      : ${feConfigs().collect{ it.name }.join(', ')}
      pre-mart FE    : ${params.fe_pre_steps  ?: '(none)'}
      per-fold FE    : ${params.fe_fold_steps ?: '(none)'}
      """.stripIndent()

    /* ---- A0. start from a published artefact instead -----------------------
     * `--from_duckdb <file>` skips the compendium, the thirteen-stage build,
     * feature engineering and causal selection: their results are unpacked from
     * the deposit and fed straight to training. This is how somebody continues
     * the work from the published artefact, without the cluster or the twelve
     * hours that produced it.
     *
     * A branch here rather than a named workflow: Nextflow's strict parser
     * refuses `-entry`, and says so — "use a param to run a named workflow from
     * the entry workflow".
     */
    /* ---- A00. build a catalogue instead ------------------------------------
     * `--build_catalogue true` builds a mutations catalogue with catomatic and
     * stops. It runs no feature engineering, no training and no evaluation:
     * those stages answer a different question and a catalogue build should not
     * pay for them.
     *
     * Third top-level branch rather than a stage inside the pipeline, for the
     * same reason `from_duckdb` is a branch: the strict parser refuses `-entry`.
     * Params arrive as STRINGS, and "false" is truthy in Groovy, so the flag is
     * coerced rather than tested.
     */
    if( params.build_catalogue.toString().toBoolean() ) {

        if( params.db ) {
            ch_cat_db = channel.fromPath(params.db, checkIfExists: true)
        }
        else if( params.cryptic_src ) {
            ch_cat_db = BUILD_SLIM_DB(
                channel.fromPath(params.cryptic_src, checkIfExists: true)).db
        }
        else {
            error "Provide --db <cryptic.duckdb> or --cryptic_src <cryptic-tables-vX.Y.Z/>"
        }

        // The RI gate applies here for the same reason it applies to a mart: a
        // catalogue built on a database that fails its own referential
        // integrity checks is worse than none, because it looks fine.
        ch_cat_checked = CHECK_DB(ch_cat_db).summary
            .combine(ch_cat_db)
            .map { _summary, db -> db }

        // Fan out over drugs x read-support thresholds. The published sweep
        // varies FRS at a fixed background and p, so the threshold belongs to
        // the extraction and each setting is its own comparable unit.
        ch_cat_frs = params.catalogue_frs
            ? channel.fromList(params.catalogue_frs.toString().split(',').collect { it.trim() })
            : channel.of('none')

        // "RIF:rpoB;BDQ:Rv0678,atpE,pepQ" -> [RIF: 'rpoB', BDQ: 'Rv0678,atpE,pepQ']
        def geneMap = (params.catalogue_genes ?: '').toString().split(';')
            .findAll { it.trim() }
            .collectEntries { entry ->
                def (d, g) = entry.trim().split(':', 2)
                [(d.trim()): g.trim()]
            }

        ch_cat_units = channel.fromList(
                params.catalogue_drugs.toString().split(',').collect { it.trim() })
            .combine(ch_cat_frs)
            .combine(ch_cat_checked)
            .map { drug, frs, db ->
                tuple(drug, geneMap[drug], params.catalogue_dataset_tag, frs, db) }

        ch_cohort_files = EXTRACT_COHORT(ch_cat_units).cohort

        // The wildcards travel WITH the cohort: they are catomatic's default
        // rules expanded over the genes this cohort actually contains, so they
        // are derived from the same extraction rather than from a panel
        // maintained somewhere else. Override with --catalogue_wildcards when
        // the rules themselves are the subject rather than a default.
        ch_cat_ready = params.catalogue_wildcards
            ? ch_cohort_files.combine(
                  channel.fromPath(params.catalogue_wildcards, checkIfExists: true))
                .map { drug, tag, frs, sm, mu, co, _generated, supplied ->
                       tuple(drug, tag, frs, sm, mu, co, supplied) }
            : ch_cohort_files

        ch_catalogues = RUN_CATOMATIC(ch_cat_ready).catalogue

        ch_comparisons = params.catalogue_reference
            ? COMPARE_CATALOGUE(
                ch_catalogues,
                channel.fromPath(params.catalogue_reference, checkIfExists: true).first()
              ).comparison
            : channel.empty()

        // The stages this path skipped publish nothing; they were not run.
        ch_checked     = channel.empty()
        ch_cohort      = channel.empty()
        ch_marts       = channel.empty()
        ch_sel_fold    = channel.empty()
        ch_eval        = channel.empty()
        ch_tiers       = channel.empty()
        ch_models      = channel.empty()
        ch_models_fold = channel.empty()
    }
    else if( params.from_duckdb ) {
        ch_catalogues  = channel.empty()
        ch_comparisons = channel.empty()

        ch_unpacked = UNPACK_DUCKDB(
            channel.fromPath(params.from_duckdb, checkIfExists: true)).index

        // Only `.name` is ever read downstream, and it is a grouping label. A
        // deposit rerun does NOT redo feature engineering, so advertising this
        // run's fe_* knobs as though they had produced the mart would be a
        // fabrication -- the mart came from the deposit. Carry the label alone.
        def feLabel = [name: params.fe_name]

        ch_units = ch_unpacked.splitJson()
            .map { u -> tuple(u.drug, feLabel, u.held_out,
                              file(u.mart), file(u.meta),
                              file(u.concordance), file(u.selection)) }

        // The deposit holds BOTH arms, and they are not interchangeable. Units
        // with a held-out lineage carry a selection made inside that training
        // fold and are the only honest input to evaluation; the `none` unit
        // carries a selection refit on all data and is the only legitimate
        // input to a shipped model. Feeding either arm to the other stage --
        // which the first version of this branch did -- silently produces
        // deployment models fitted on fold-restricted features and a "fold"
        // score computed on rows its features already saw.
        ch_by_arm  = ch_units.branch { _d, _fe, held_out, _m, _mt, _c, _s ->
            deployment: held_out == 'none'
            evaluation: true
        }

        ch_models      = TRAIN_H2O(ch_by_arm.deployment)
        ch_models_fold = TRAIN_H2O_FOLD(ch_by_arm.evaluation)
        ch_eval   = EVALUATE_CV(ch_by_arm.evaluation)
        ch_tiers  = TIER_REPORT(
            ch_eval.map { drug, fe, _held_out, manifest -> tuple(drug, fe.name, manifest) }
                   .groupTuple(by: [0, 1]))

        // The stages this path skipped publish nothing; they were not run.
        ch_checked  = channel.empty()
        ch_cohort   = channel.empty()
        ch_marts    = channel.empty()
        ch_sel_fold = channel.empty()
    }
    else {

    ch_catalogues  = channel.empty()
    ch_comparisons = channel.empty()

    /* ---- A. the slim database -------------------------------------------
     * Thirteen stages that all mutate ONE 3.5 GB DuckDB file in place, so they
     * cannot fan out and staging the file between thirteen tasks would move
     * ~45 GB to save ~14 minutes. They run as one cached unit instead. Pass
     * --db to skip this entirely and reuse a database you already have, which
     * is the common case: the DB changes on a CRyPTIC release, the FE changes
     * every afternoon.
     */
    if( params.db ) {
        ch_db = channel.fromPath(params.db, checkIfExists: true)
    }
    else if( params.cryptic_src ) {
        ch_db = BUILD_SLIM_DB(channel.fromPath(params.cryptic_src, checkIfExists: true)).db
    }
    else {
        error "Provide --db <cryptic-slim.duckdb> or --cryptic_src <cryptic-tables-vX.Y.Z/>"
    }

    // Referential integrity is a GATE, not a report. A mart built on a database
    // that fails its own RI checks is worse than no mart, because it looks fine.
    //
    // The gate no longer passes the database through. Re-emitting an input as an
    // output cannot work on an S3 work dir, where staged inputs are swept after
    // the task: the file is gone before the outputs are collected. Ordering is
    // instead enforced by joining the ORIGINAL channel to the gate's summary,
    // which is a real output -- downstream still cannot start until the check
    // has passed, and nothing has to survive the sweep.
    if( params.skip_db_check ) {
        ch_checked = ch_db
    }
    else {
        ch_checked = ch_db.combine( CHECK_DB(ch_db).summary ).map { db, _summary -> db }
    }

    /* ---- B. cohort ------------------------------------------------------ */
    ch_cohort = BUILD_COHORT(ch_checked)

    /* ---- C. feature marts, one per drug × FE configuration --------------
     * This is the fan-out that makes FE a parameter. params.fe is a single map
     * by default; pass a list of maps to sweep several configurations against
     * identical downstream folds.
     */
    ch_marts = BUILD_MART(
        channel.fromList(params.drugs.tokenize(','))
               .combine(channel.fromList(feConfigs()))
               .combine(ch_checked)
               .combine(ch_cohort)
    )

    /* ---- C2. label-free FE layers ---------------------------------------
     * Universe filters, annotation joins and derived columns. None of these can
     * see the phenotype, so they run ONCE per mart and the result is valid for
     * every fold.
     */
    ch_annot = channel.fromPath(params.gene_annotations, checkIfExists: true)
    ch_pre = FE_PRE(ch_marts.map { d, fe, n, mart, meta -> tuple(d, fe, n, mart, meta, 'none') },
                    'pre-mart', params.fe_pre_steps, ch_annot).mart
    ch_ready = params.fe_pre_steps?.trim()
        ? ch_pre.map { d, fe, n, mart, meta, ho, notes -> tuple(d, fe, n, mart, meta) }
        : ch_marts

    /* ---- D. concordance selection ---------------------------------------
     * Two arms, deliberately distinguished.
     *
     * The EVALUATION arm runs selection inside each training fold, so the
     * feature set never sees the held-out lineage. Without this, a lineage-LOO
     * AUC is not honest: the features were chosen while looking at the test
     * rows. This is the reason the pipeline exists at all — it is 11 drugs ×
     * 4 lineages = 44 independent EBM + causal-forest + refutation runs, which
     * is unaffordable in a loop on one machine and unremarkable distributed.
     *
     * The DEPLOYMENT arm refits selection on all data, which is legitimate for
     * a shipped model and is labelled as such so the two are never conflated.
     */
    ch_fold_inputs = ch_ready.flatMap { d, fe, n, mart, meta ->
        params.lineages.tokenize(',').collect { lin -> tuple(d, fe, n, mart, meta, lin) } }
    // Label-using FE runs INSIDE the fold, so the ranking never sees the lineage
    // it will be tested on.
    ch_fold_fe = params.fe_fold_steps?.trim()
        ? FE_FOLD(ch_fold_inputs, 'per-fold', params.fe_fold_steps, ch_annot).mart
              .map { d, fe, n, mart, meta, ho, notes -> tuple(d, fe, n, mart, meta, ho) }
        : ch_fold_inputs

    ch_sel_fold = SELECT_CONCORDANT(ch_fold_fe, ch_checked.first())
    ch_sel_full = SELECT_FULL(ch_ready.map { d, fe, n, mart, meta ->
        tuple(d, fe, n, mart, meta, 'none') }, ch_checked.first())

    /* ---- D2. the publishable artefact (optional) -------------------------
     * The state of the data after feature engineering and before training, as
     * one self-contained DuckDB. Nothing downstream reads it: TRAIN_H2O still
     * takes the mart and the selection directly, so enabling this changes what
     * is deposited, never what is computed.
     *
     * Gated on params.create_duckdb, so the default run does not pay for it.
     */
    if (params.create_duckdb) {
        CREATE_DUCKDB(
            ch_ready.map { d, fe, n, mart, meta -> mart }.collect(),
            ch_ready.map { d, fe, n, mart, meta -> meta }.collect(),
            ch_sel_fold.map { _d, _fe, _ho, _mart, _meta, conc, _sel -> conc }
                       .mix(ch_sel_full.map { _d, _fe, _ho, _mart, _meta, conc, _sel -> conc })
                       .collect(),
            // The manifests say HOW each selection was made. Without them the
            // deposit shows which features were chosen but not on what basis,
            // which is the part a reader cannot reconstruct.
            ch_sel_fold.map { _d, _fe, _ho, _mart, _meta, _conc, sel -> sel }
                       .mix(ch_sel_full.map { _d, _fe, _ho, _mart, _meta, _conc, sel -> sel })
                       .collect()
        )
    }

    /* ---- E. honest evaluation ------------------------------------------- */
    ch_eval = EVALUATE_CV(ch_sel_fold)

    // One operating range per (drug, FE configuration) -- grouping on the FE
    // name too, so a sweep does not silently average two configurations into
    // one tier.
    ch_tiers = TIER_REPORT(
        ch_eval.map { drug, fe, held_out, manifest -> tuple(drug, fe.name, manifest) }
               .groupTuple(by: [0, 1])
    )

    /* ---- F. training ----------------------------------------------------
     * H2O boots INSIDE the task and is torn down with it. A long-lived server
     * would be shared mutable state across tasks: one run's memory pressure
     * becomes another's failure, and -resume stops meaning anything. Per-task
     * means the JVM heap is bounded by the Nomad allocation that was actually
     * reserved for it.
     */
    ch_models = TRAIN_H2O(ch_sel_full)

    /* ---- F2. the same candidates, on a lineage no fold contained ---------
     * The leaderboard of F ranks models by out-of-fold AUC within the training
     * lineages, and EVALUATE_CV reports a held-out lineage using one fixed
     * gradient-boosted baseline. Those are different models on different rows,
     * so neither answers whether the H2O apparatus earns its place. This arm
     * trains the same stack per fold and scores every candidate on the withheld
     * lineage, which is the axis EVALUATE_CV already reports on.
     */
    ch_models_fold = TRAIN_H2O_FOLD(ch_sel_fold)

    }

    publish:
    database = ch_checked
    cohort   = ch_cohort
    marts    = ch_marts
    folds    = ch_sel_fold
    evals    = ch_eval
    tiers    = ch_tiers
    models   = ch_models.mix(ch_models_fold)
    catalogues  = ch_catalogues
    comparisons = ch_comparisons
}


/* These channels emit TUPLES, so a path closure receives a List. `r.drug` on a
 * List is a Groovy spread -- it asks every element for a `.drug`, hits the drug
 * String first and throws MissingPropertyException. It cost a full campaign:
 * every task ran to completion and the run then died publishing its results.
 * (`r.fe_hash` never existed either; the FE map carries `name`.)
 *
 * Index positionally against each process's declared output tuple:
 *   marts   BUILD_MART        (drug, fe, fe.name, mart, meta)
 *   folds   SELECT_*          (drug, fe, held_out, mart, meta, conc, sel)
 *   evals   EVALUATE_CV       (drug, fe, held_out, manifest)
 *   models  TRAIN_H2O         (drug, fe, h2o/**, manifest)
 */
output {
    database { path 'database' }
    cohort   { path 'cohort' }
    marts    { path { r -> "marts/${r[0]}/${r[2]}" } }
    folds    { path { r -> "selection/${r[0]}/held-out-${r[2]}" } }
    evals    { path { r -> "evaluation/${r[0]}" } }
    tiers    { path 'tiers' }
    // held_out is r[2]: the deployment refit and four fold models are the same
    // drug and would otherwise publish to one path and overwrite each other.
    models   { path { r -> "models/${r[0]}/held-out-${r[2]}" } }
    // catalogues  RUN_CATOMATIC     (drug, dataset_tag, frs, csv, json, build)
    // comparisons COMPARE_CATALOGUE (drug, dataset_tag, frs, comparison)
    // One directory per cohort tag, drug and read-support setting: the sweep
    // publishes many catalogues for one drug and they must not overwrite.
    catalogues  { path { r -> "catalogues/${r[1]}/${r[0]}/frs-${r[2]}" } }
    comparisons { path { r -> "catalogues/${r[1]}/${r[0]}/frs-${r[2]}" } }
}
