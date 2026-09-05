# cryptic-db

CRyPTIC release tables into a reproducible, checksummed DuckDB. The output is the input to the
[training pipeline](../pipeline), which consumes it with `--db`.

**The artefact is the full database.** The thirteen shell stages load every table in the release —
about 19 GB, 80M mutation rows, and an `effects` table carrying `EVIDENCE` and `PREDICTION_VALUES`.
Feature engineering runs on that. `cryptic-slim.duckdb` was a speed hack for iterating on the
analysis, and it is not free: it drops those two columns and the plate-level tables, which blocks
three catalogued techniques outright —

| technique | needs |
|---|---|
| T2 RFUS / SOLO grading | `effects.PREDICTION_VALUES` |
| T3 EVIDENCE JSON — observed samples, confidence grading, WHO HGVS | `effects.EVIDENCE` |
| T4 real WHO tier grading (`cat__n_tier1/tier2/interim`) | `effects.EVIDENCE` |

`build_mart` currently substitutes per-mutation R/S/U/F counts for T4 and says why: *"the slim DB's
effects table is trimmed … so tier grading isn't recoverable."* On the full database it is.

Pass `--make_slim` for a reduced copy when the wait is what you are optimising. It records what it
dropped, in a `slim_provenance` table and a manifest, so a mart built from it cannot silently claim
capabilities it does not have.

```bash
nextflow run . --cryptic_src /path/to/cryptic-tables-v3.4.0
```

## Why it is a separate pipeline

The database changes when CRyPTIC cuts a release. The feature engineering changes every
afternoon. Rebuilding a 3.5 GB database to try a different top-N is the wrong unit of work, so
this emits a versioned artefact and the training pipeline consumes it.

## What it does that the shell build did not

**Identifies the inputs by content.** `VERIFY_SOURCE` checksums every input parquet and CSV
before anything reads them, and the manifest is written into the database as a
`source_provenance` table. A release directory edited in place is otherwise indistinguishable
from one that has not been, and every downstream provenance claim inherits that uncertainty.

**Runs the loads concurrently.** Stages 01–09 each read only the source parquet they load and
the tables they themselves create — verified against the scripts, not assumed — so they run as
independent shards, each writing its own DuckDB, and are attached together afterwards. DuckDB
permits one writer per file, which is why sharding is the way to parallelise it at all.

**Fails on referential integrity.** Stage 11 created the check views but never failed on them,
so a database with orphaned rows could go on to produce marts that look perfectly well formed.
Here a violation stops the run.

**Exports a parquet mirror**, so the release is readable without DuckDB — by another language,
another tool, or a reviewer who should not have to install anything to look at the data.

## Stages

```
VERIFY_SOURCE ──► LOAD_SHARD × 9 ──► ASSEMBLE_DB ──► BUILD_VIEWS ──► WRITE_METADATA
   checksums        concurrent          attach          stage 10        stage 12 +
                                        + copy                          provenance
                                                                             │
                                        EXPORT_PARQUET ◄── CHECK_INTEGRITY ◄──┘
                                          mirror              stage 11, as a GATE
```

The shell stage scripts are reused unmodified. They already take `(src, db)` positionally, so
the pipeline changes where they write, not what they do.

## Outputs

| path | what |
|---|---|
| `cryptic.duckdb` | the full database, carrying its own `source_provenance` table |
| `slim/cryptic-slim.duckdb` | optional reduced copy, with a manifest of what it dropped |
| `provenance/source_manifest.json` | release, and a SHA-256 for every input file |
| `checks/ri_summary.tsv` | the referential-integrity result |
| `parquet/*.parquet` | one file per table |

## Next: end-to-end validation on the cluster

Both pipelines have been linted, previewed and exercised stage by stage on a laptop
against the already-built slim database. What has **not** been done is a cold run from the
CRyPTIC release parquet files through to feature engineering, on cluster hardware.

That is the next validation, on the cluster:

1. Stage the full CRyPTIC release parquet files (v3.4.0, 16 files, 1.3 GB).
2. Run `cryptic-db` cold — `VERIFY_SOURCE` through the parquet mirror — and confirm the
   nine load shards genuinely run concurrently, that `ASSEMBLE_DB` reproduces the table
   set the sequential shell build produced, and that `CHECK_INTEGRITY` passes on a
   database this pipeline built rather than one it inherited.
3. Feed the resulting database to the training pipeline and run the FE layers, then the
   44 fold-internal selections, under `-profile nomad`.

Two things only a cold run can tell us: whether the sharded assembly is byte-equivalent
to the sequential build, and what the concurrency actually buys once the large loads
(`08_mutations`, `09_variants`) compete for I/O.
