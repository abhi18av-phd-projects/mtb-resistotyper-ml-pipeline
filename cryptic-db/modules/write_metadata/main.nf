/*
 * Stage 12 — build metadata, plus the source manifest folded in.
 *
 * The database carries its own provenance so a mart built from it can quote the
 * exact release and input checksums without the operator having to remember
 * which directory was on disk that week.
 */
process WRITE_METADATA {
    tag "metadata"
    label 'load'

    input:
    path db
    path source_manifest

    output:
    path db, emit: db

    script:
    """
    bash ${params.scripts_dir}/12_metadata.sh "${params.cryptic_src}" ${db} > metadata.log 2>&1

    ${params.python} - <<'PY'
import duckdb, json, pathlib
manifest = json.loads(pathlib.Path("${source_manifest}").read_text())
con = duckdb.connect("${db}")
con.execute("CREATE OR REPLACE TABLE source_provenance (key VARCHAR, value VARCHAR)")
con.executemany(
    "INSERT INTO source_provenance VALUES (?, ?)",
    [("release", str(manifest.get("release"))),
     ("source_dir", str(manifest.get("source_dir"))),
     ("n_files", str(manifest.get("n_files"))),
     ("files_sha256", json.dumps(manifest.get("files", {})))],
)
con.close()
PY
    """
}
