// Unpack a published DuckDB into the files the training stages already expect.
//
// The entry point for somebody who downloaded the artefact and wants to train
// from it. It reconstructs feature_mart_<DRUG>.parquet, its metadata, and one
// concordance.parquet + selection.json per fold, so TRAIN_H2O is fed the same
// tuple whether the selection ran here or came out of a deposit.

process UNPACK_DUCKDB {
    tag "${db.name}"
    label 'mart'

    input:
    path db

    output:
    path "unpacked/units.json", emit: index
    path "unpacked/**",         emit: files

    script:
    def drugs = params.drugs ? "--drugs ${params.drugs}" : ''
    def arm   = params.from_duckdb_arm ? "--arm ${params.from_duckdb_arm}" : ''
    """
    cd ${params.project_root} && \\
    ${params.python} -m analysis.scripts.feature_mart.unpack_duckdb \\
        --db \$OLDPWD/${db} \\
        --out \$OLDPWD/unpacked ${drugs} ${arm}
    """
}
