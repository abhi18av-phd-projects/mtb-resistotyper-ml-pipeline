# Building the CRyPTIC DuckDB database

This directory contains numbered scripts that turn the [CRyPTIC consortium](https://www.crypticproject.org/) parquet release into a single, queryable DuckDB database with proper keys, indexes, and integrity checks.

> **Output:** `cryptic.duckdb` (~19 GB), 18 tables, 175 M+ rows, ~14 min wall time.

---

## Quick start

```bash
# One-shot build
./run_all.sh \
    /path/to/cryptic-tables-v3.4.0 \
    /path/to/cryptic.duckdb

# or via env vars
export MTB_CRYPTIC_SRC=/path/to/cryptic-tables-v3.4.0
export MTB_CRYPTIC_DB=/path/to/cryptic.duckdb
./run_all.sh

# Run a single stage (any number, any order, but earlier stages must have run first)
./08_mutations.sh "$MTB_CRYPTIC_SRC" "$MTB_CRYPTIC_DB"
```

The scripts use only the `duckdb` CLI — no Python deps. Install with `brew install duckdb` or download from [duckdb.org/docs/installation](https://duckdb.org/docs/installation/).

---

## Script index

| # | Script | What it loads | Rows | Wall time |
|---|---|---|---:|---:|
| 00 | `00_init_db.sh` | (clears + creates empty DB) | — | <1 s |
| 01 | `01_reference.sh` | drugs, countries, sites, plate_layout | ~640 | <1 s |
| 02 | `02_sample_anchors.sh` | wgs_samples, dst_samples | 120 k | 1 s |
| 03 | `03_one_to_one.sh` | genomes, ukmyc_plates | 81 k | 1 s |
| 04 | `04_effects_predictions.sh` | effects, predictions | ~2 M | 13 s |
| 05 | `05_phenotypes.sh` | dst_measurements, ukmyc_phenotypes | ~1 M | 7 s |
| 06 | `06_bashthebug.sh` | bashthebug, bashthebug_classifications | ~5 M | 16 s |
| 07 | `07_ukmyc_growth.sh` | ukmyc_growth | 3 M | 10 s |
| 08 | `08_mutations.sh` | **mutations** | 80 M | ~7 min |
| 09 | `09_variants.sh` | **variants** | 84 M | ~4 min |
| 10 | `10_views.sh` | v_catalogues, v_sample_coverage, v_sample_fanout | — | <1 s |
| 11 | `11_ri_checks.sh` | 14 × `v_ri_*` integrity views | — | <1 s |
| 12 | `12_metadata.sh` | _database_metadata | 5 | <1 s |
| — | `_lib.sh` | shared helpers (sourced by every stage) | — | — |
| — | `run_all.sh` | sequential orchestrator | — | — |

---

## What this database represents — the conceptual model

The CRyPTIC release is a network of related observations about *Mycobacterium tuberculosis* isolates. The shape that links them all is **`UNIQUEID`**: a string like `site.06.subj.06TB_0269.lab.06MIL1372.iso.1` that identifies one bacterial isolate. Every fact table is a fan-out from `UNIQUEID`.

### The two anchor populations

There are **two** parallel anchor tables, and they overlap but neither contains the other:

```
                    UNIQUEID
                  /          \
        wgs_samples           dst_samples
        (54,565)              (65,842)
        sequenced isolate     clinically-tested isolate
        intersection: 54,354 isolates appear in both
```

- **`wgs_samples`** — every isolate that went through Pathogena/Clockwork to produce a VCF. Parent of all *genetic* tables.
- **`dst_samples`** — every isolate that has clinical drug-susceptibility test data. Parent of `dst_measurements`.

ML pipelines that need genotype + phenotype implicitly take the intersection (`INNER JOIN` on `UNIQUEID`).

### The fan-out (rows per isolate)

| Child table | Rows / sample | What each row represents |
|---|---:|---|
| `genomes` | 1 | species + lineage call |
| `ukmyc_plates` | 0–1 | provenance of the microtitre plate (only for samples tested on UKMYC) |
| `predictions` | ~15 | one per drug — the catalogue's S/R/U call |
| `ukmyc_phenotypes` | ~13 | one per drug — the experimentally measured MIC + binary call |
| `dst_measurements` | ~10 | one per (drug, method, source) — clinical MIC measurements |
| `effects` | ~21 | per (drug, gene, mutation) row that contributed evidence to a `predictions` call |
| `bashthebug` | varies | one per plate image read by citizen scientists |
| `ukmyc_growth` | ~150 | one per well-dilution-day reading |
| **`mutations`** | **~1,483** | one per gene-level mutation called by gnomonicus |
| **`variants`** | **~1,553** | one per VCF record |

The two large tables (`mutations`, `variants`) hold > 97 % of the rows. They are the *features* in any resistotyping ML problem.

---

## The relationship graph

```
                      ┌──────────────┐
                      │     drugs    │ ◄── effects, predictions,
                      │  (39 rows)   │     dst_measurements,
                      └──────┬───────┘     ukmyc_phenotypes,
                             │             ukmyc_growth,
                             │             bashthebug
                             ▼
   ┌────────────┐    ┌──────────────────┐    ┌──────────────┐
   │  countries │ ◄─ │    wgs_samples   │    │ dst_samples  │ ─► countries
   └────────────┘    │   (PK UNIQUEID)  │    │ (PK UNIQUEID)│
                     └─────────┬────────┘    └──────┬───────┘
                               │                    │
       ┌───────────────────────┼────────────┐       │
       │           │           │            │       │
       ▼           ▼           ▼            ▼       ▼
  ┌─────────┐ ┌─────────┐ ┌────────┐ ┌──────────────┐ ┌──────────────────┐
  │ genomes │ │mutations│ │variants│ │ ukmyc_plates │ │ dst_measurements │
  │ (1:1)   │ │ (1:N)   │ │ (1:N)  │ │   (1:0..1)   │ │     (1:N)        │
  └─────────┘ └─────────┘ └────────┘ └──────┬───────┘ └──────────────────┘
                                            │
                  ┌─────────────────────────┼─────────────────┐
                  │                         │                 │
                  ▼                         ▼                 ▼
            ┌───────────────┐       ┌──────────────┐    ┌────────────┐
            │ukmyc_phenoty- │       │ ukmyc_growth │    │ bashthebug │
            │     pes       │       │              │    │            │
            └───────────────┘       └──────┬───────┘    └─────┬──────┘
                                           │                  │
                                           ▼                  ▼
                                   ┌─────────────┐    ┌────────────────┐
                                   │plate_layout │    │bashthebug_     │
                                   │             │    │classifications │
                                   └─────────────┘    └────────────────┘

  effects, predictions ─── catalogue-versioned: every row carries
                            (CATALOGUE_NAME, CATALOGUE_VERSION) — filter
                            explicitly before joining into the model frame.
```

---

## Why these particular joins

### 1. The genetics → label join (the core ML query)

```sql
SELECT m.UNIQUEID, m.GENE, m.MUTATION, p.PREDICTION
FROM mutations m
INNER JOIN predictions p USING (UNIQUEID)
INNER JOIN wgs_samples w  USING (UNIQUEID)
WHERE p.DRUG = 'RIF'
  AND p.CATALOGUE_NAME    = 'WHO-UCN-GTB-PCI-2023.5'
  AND p.CATALOGUE_VERSION = 2.0
  AND p.PREDICTION IN ('R', 'S')
  AND w.status = 'complete'
```

Why each clause is there:

- **`mutations` is the feature side, `predictions` is the label side.** This is the simplest sample-level genotype → drug-resistance join.
- **`INNER JOIN wgs_samples`** is a quality guard: not every UNIQUEID in `predictions` reaches a "complete" Pathogena run; some are `cannot_assemble` or `cannot_speciate`. Those samples have unreliable mutations and should not be used as training data.
- **Filter on `(CATALOGUE_NAME, CATALOGUE_VERSION)`** even when only one is present today. Future releases re-add WHO1, and an unfiltered join would silently double-count.
- **Filter `PREDICTION IN ('R', 'S')`** drops `U` (unknown) and `F` (failed). For binary classification the model needs unambiguous labels.

### 2. The drug filter belongs early

Only ~22 of the 39 drugs in `drugs.drug` actually have predictions. Filtering `WHERE p.DRUG = 'RIF'` at the top of the CTE avoids materialising the whole 810k-row predictions table into the join.

### 3. The catalogue join (when comparing catalogues)

```sql
SELECT * FROM v_catalogues;
-- Pick which (CATALOGUE_NAME, CATALOGUE_VERSION) to use, then drive everything off that.
```

`v_catalogues` is the right place to discover what's loaded — never assume.

### 4. The genotype + clinical-DST join (gold-standard label)

```sql
-- "All samples with a sequencing-derived prediction AND a clinical DST result"
SELECT p.UNIQUEID, p.DRUG,
       p.PREDICTION                          AS catalogue_call,
       d.PHENOTYPE                           AS clinical_call,
       d.METHOD_1
FROM predictions p
INNER JOIN dst_measurements d USING (UNIQUEID, DRUG)
WHERE p.CATALOGUE_NAME = 'WHO-UCN-GTB-PCI-2023.5'
  AND p.CATALOGUE_VERSION = 2.0;
```

This is the join you use to *evaluate* the catalogue against ground truth. Note the `USING (UNIQUEID, DRUG)` — both columns are part of the natural join condition.

### 5. Per-well-level data (rare but useful)

```sql
SELECT g.UNIQUEID, g.DRUG, g.READINGDAY, g.DILUTION,
       g.GROWTH,
       l.CONC, l.BINARY_PHENOTYPE
FROM ukmyc_growth g
INNER JOIN plate_layout l USING (PLATEDESIGN, DRUG, DILUTION);
```

The composite `(PLATEDESIGN, DRUG, DILUTION)` is the natural FK. `plate_layout` resolves it to the actual concentration and the design's binary R/S call.

### 6. Why `mutations` + `variants` are kept separate

They look similar but are computed at different layers:

- **`variants`** = "what's in the VCF" (raw genotype, base-pair level).
- **`mutations`** = "what gnomonicus interprets that VCF as in terms of named gene mutations" — the *biological* layer.

For catalogue-aware features, join on `mutations`. For raw-genotype experiments (e.g., training your own catalogue), join on `variants`.

---

## Schema decisions worth understanding

### Surrogate `row_id` PKs on four tables

`mutations`, `variants`, `effects`, `dst_measurements`, `ukmyc_growth` carry a generated `row_id` PK rather than the natural composite. Why:

| Table | Natural near-key | Why it doesn't work as PK |
|---|---|---|
| `mutations` | (UNIQUEID, GENE, MUTATION, IS_MINOR) | ~469 rows in 80 M differ only in quality fields |
| `variants` | (UNIQUEID, GENE, VARIANT, MINOR_VARIANT) | similar |
| `effects` | (UNIQUEID, DRUG, GENE, MUTATION, CAT, VER) | `GENE` can be NULL for multi-gene mutations like `fabG1@G241G&fabG1@c723t` |
| `dst_measurements` | (UNIQUEID, DRUG, SOURCE) | repeated measurements per source |
| `ukmyc_growth` | (UNIQUEID, DRUG, READINGDAY, DILUTION, SITEID, PLATEDESIGN) | plates re-imaged on multiple days |

Surrogate PK keeps every row addressable; secondary indexes on the natural-key columns restore join speed.

### Catalogue versioning is real even when latent

`effects` and `predictions` carry `(CATALOGUE_NAME, CATALOGUE_VERSION)` even though v3.4.0 currently has only WHO2 (`WHO-UCN-GTB-PCI-2023.5 v2.0`). Past releases held WHO1 + WHO2 simultaneously and future releases likely will again. **Filter on the catalogue tuple in every query** — it's a one-line guard against silent double-counting on the next data refresh.

### Foreign keys: declared as integrity views, not constraints

DuckDB v1.x parses FK constraints but does not reliably enforce them across bulk-loaded fact tables. Instead of declaring FKs that might be silently ignored, this build creates 14 `v_ri_*` views that *materially compute* orphan rows. After a build, run:

```sql
SELECT * FROM v_ri_mutations_to_wgs LIMIT 10;
-- (should return zero rows on a clean cryptic-tables-v3.4.0 build)
```

A clean genetic side has 0 orphans. The clinical side has known issues:

- `v_ri_ukmyc_growth_to_plates` ⇒ 43,872 rows / 416 distinct UIDs
- `v_ri_bashthebug_to_plates` ⇒ 1,190 rows / 214 distinct UIDs

These are pre-existing data-quality issues in the source release; downstream code that joins through `ukmyc_plates` should `LEFT JOIN` and tolerate `NULL`s.

### Indexing strategy

Each fact table has indexes on the columns you'll *actually* join on:

- `(UNIQUEID)` — for "all rows for this sample"
- `(GENE, MUTATION)` or `(GENE, VARIANT)` — for "all samples carrying X"
- `(UNIQUEID, DRUG)` — for the most common ML join
- `(PLATEDESIGN, DRUG, DILUTION)` — for plate-layout joins

Indexing the 80 M-row tables takes longer than loading them. That's fine — it's a one-time cost. Without these indexes, a single feature-engineering query takes minutes; with them it takes seconds.

---

## Using the database

### Recommended client setup

```python
import duckdb
con = duckdb.connect("cryptic.duckdb", read_only=True)

# Discover what's loaded
con.sql("SELECT * FROM v_catalogues").show()
con.sql("SELECT * FROM _database_metadata ORDER BY key").show()

# Check coverage
con.sql("""
    SELECT SUM(has_wgs::INT) AS n_wgs,
           SUM(has_dst::INT) AS n_dst,
           SUM((has_wgs AND has_dst)::INT) AS intersection
    FROM v_sample_coverage
""").show()
```

### Pulling a feature frame for ML

```python
sql = """
WITH lab AS (
    SELECT p.UNIQUEID,
           CASE p.PREDICTION WHEN 'R' THEN 1 WHEN 'S' THEN 0 END AS y
    FROM predictions p
    INNER JOIN wgs_samples w USING (UNIQUEID)
    WHERE p.DRUG = 'RIF'
      AND p.CATALOGUE_NAME    = 'WHO-UCN-GTB-PCI-2023.5'
      AND p.CATALOGUE_VERSION = 2.0
      AND p.PREDICTION IN ('R', 'S')
      AND w.status = 'complete'
),
top_mut AS (
    SELECT GENE || '_' || MUTATION AS gm
    FROM mutations m INNER JOIN lab USING (UNIQUEID)
    GROUP BY 1
    ORDER BY COUNT(DISTINCT UNIQUEID) DESC
    LIMIT 1000
)
PIVOT (SELECT m.UNIQUEID, m.GENE || '_' || m.MUTATION AS gm
       FROM mutations m
       INNER JOIN lab USING (UNIQUEID)
       INNER JOIN top_mut USING (gm))
ON gm USING COUNT(*) > 0
GROUP BY UNIQUEID
"""
features = con.sql(sql).pl()   # zero-copy to polars; .df() for pandas
```

This single SQL replaces the entire Clojure feature-engineering stack the project used previously.

---

## Nextflow integration

The scripts are designed to drop into a Nextflow pipeline as separate processes. The simplest pattern:

```nextflow
nextflow.enable.dsl = 2

params.cryptic_src = "/data/cryptic-tables-v3.4.0"

workflow {
    Channel.fromPath(params.cryptic_src) | build_cryptic_db
}

process build_cryptic_db {
    publishDir "results", mode: "copy"
    cpus 4
    memory "32 GB"

    input:
        path src_dir

    output:
        path "cryptic.duckdb"

    script:
    """
    bash ${projectDir}/analysis/scripts/build_cryptic_db/run_all.sh ${src_dir} cryptic.duckdb
    """
}
```

If you want each stage as its own process (independently cacheable, fail-isolated, resumable):

```nextflow
process stage_load {
    tag "${stage}"
    input:
        tuple val(stage), path(src_dir), path(db_in)
    output:
        path db_in, includeInputs: true   // mutates the DB in-place
    script:
    """
    bash ${projectDir}/analysis/scripts/build_cryptic_db/${stage} ${src_dir} ${db_in}
    """
}

workflow {
    db = init_db(params.cryptic_src)
    stages = Channel.of(
        "01_reference.sh",
        "02_sample_anchors.sh",
        "03_one_to_one.sh",
        "04_effects_predictions.sh",
        "05_phenotypes.sh",
        "06_bashthebug.sh",
        "07_ukmyc_growth.sh",
        "08_mutations.sh",
        "09_variants.sh",
        "10_views.sh",
        "11_ri_checks.sh",
        "12_metadata.sh"
    )
    // Sequential reduction — each stage produces the DB consumed by the next
    stages.reduce(db) { acc, s ->
        stage_load(tuple(s, params.cryptic_src, acc))
    }
}
```

Note that the DuckDB *file* is the artifact passing between stages, not separate per-stage outputs — every stage reads and mutates the same file. Caching at the stage level still works (Nextflow caches by input hash), but you get less out of it because the input file changes after every stage. The single-process version above is usually the right choice; the per-stage version helps when running on a cluster where you want to parallelise the indexing of the two big tables (08, 09) on a high-memory node.

### Resource recommendations

- `08_mutations.sh` and `09_variants.sh` are the only stages that need real resources. Give them ≥16 GB RAM and 4–8 cores.
- All other stages run comfortably in <2 GB RAM.
- Disk: needs ~25 GB free (the output is 19 GB; intermediates push it to ~22 GB at peak).

---

## Reproducibility

The build is fully reproducible from the parquet/CSV inputs. There is no random sampling, no Python kernel state, no order-dependent decisions. Re-running `run_all.sh` with the same source directory always produces the same database byte-for-byte (modulo `build_timestamp_utc` in `_database_metadata`).

To upgrade to a new CRyPTIC release:

1. Drop the new `cryptic-tables-vX.Y.Z` directory next to the old one.
2. Re-run `./run_all.sh /path/to/cryptic-tables-vX.Y.Z /path/to/cryptic.duckdb`.
3. Inspect `_database_metadata` to confirm the new version landed.
4. Inspect `v_catalogues` to see if any new (CATALOGUE_NAME, CATALOGUE_VERSION) rows appeared — if so, downstream filters need updating.

---

## Troubleshooting

| Symptom | Cause | Fix |
|---|---|---|
| `Constraint Error: Data contains duplicates on indexed column(s)` | New release added rows that break a declared composite PK. | Investigate which table; if the dups are real, switch to surrogate `row_id` PK following the same pattern as `mutations`. |
| `NOT NULL constraint failed` on a PK | New release added rows with NULL in a PK column. | Same — switch to surrogate PK. |
| Stage 08 / 09 OOMs | Default DuckDB memory limit too low. | Edit the stage to add `SET memory_limit='32GB';` at the top of the heredoc, or set `DUCKDB_MEMORY_LIMIT` if you've patched DuckDB to honour env vars. |
| `Catalog Error: Table 'X' does not exist` in a stage | Earlier stage was skipped. | Stages have implicit dependencies — run them in numbered order or use `run_all.sh`. |
| RI check has new orphans on a fresh release | Source data quality changed. | Inspect the offending `v_ri_*` view; add a documented exception in `12_metadata.sh` if the orphan is intentional. |
