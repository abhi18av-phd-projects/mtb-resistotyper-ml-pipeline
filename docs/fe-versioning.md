# Versioning feature engineering

A trained model is a function of three things: the data, the feature engineering
(FE) that turned the data into features, and the learner. The pipeline's own
version and the CRyPTIC release pin the first and the code; this scheme pins the
second, independently, so that "which features did this model see?" has an
answer that survives the model being copied somewhere else.

## Three identifiers

| identifier | what it names | how it is made | changes when |
|---|---|---|---|
| `fe_version` | the FE **technique**, as released | declared in `FE_VERSION`, semantic | a person bumps it; the lock makes that mandatory when FE code changes |
| `fe_id` | the FE **configuration** actually used | `fe-` + a hash of version, settings, determinant genes, FE steps, selection settings | any of those differ |
| `mart_id` | this FE **on this database** | `mart-` + a hash of `fe_id` and the database checksum | the FE or the database changes |

`fe_id` deliberately excludes the database. The same technique applied to
CRyPTIC v2.1.2 and v3.4.0 carries the **same** `fe_id` and a different
`mart_id` -- so a comparison across releases can require equal FE, and a
difference between two models with equal `fe_id` is a difference in the data.

`fe_id` also excludes the source code's hash. The code is versioned by
`fe_version` instead: hashing it would give two builds that produce identical
features different ids after any edit, even to a comment. The code hashes are
still recorded, in `fe_config`, for provenance.

## What counts as a MAJOR, MINOR or PATCH change

- **MAJOR** -- features are encoded or named differently, or the fold scheme
  changes. A model trained under one cannot read the other's input.
- **MINOR** -- the feature space changes but the encoding does not: a new
  covariate or FE step, a changed selection rule, a changed default such as
  `top_n`.
- **PATCH** -- output is unchanged on a fixed reference release (CRyPTIC
  v3.4.0). Verify by rebuilding a mart with the old and new code and comparing.

Each bump gets an entry in `analysis/scripts/feature_mart/FE_CHANGELOG.md`
saying what changed and how it was verified.

## Keeping the version honest

`analysis/scripts/feature_mart/FE_LOCK.json` records the FE version together with
a hash of every FE source file (the mart builder, the selection stage, every FE
step). The check

    python -m analysis.scripts.feature_mart.fe_identity --check

fails when that code has changed but `FE_VERSION` has not. After bumping and
writing the changelog entry:

    python -m analysis.scripts.feature_mart.fe_identity --lock

The previous, informal FE version (`v1.0.2`) was carried as the mart filename
label, which every real run overwrote with `vFULL` -- so it never moved and the
shipped models recorded nothing about their FE. The lock exists so that cannot
recur silently.

## Where the identity is recorded

The stages that apply FE record their own part, so a consumer never re-derives
it from code it may not share:

- **mart sidecar** (`feature_mart_<DRUG>_*.metadata.json`) -- `fe_version` and an
  `fe` block: the mart-builder settings, the complete determinant-gene map, and
  the builder's code hash.
- **selection manifest** -- an `fe` block with the concordance settings and the
  selection stage's code hash. The realised candidate-pool size is data, not
  configuration, and is left out.

Consumers assemble `fe_id` from those blocks and carry it on:

- **the DuckDB deposit** (`--create_duckdb`) -- `fe_version`, `fe_id` and
  `mart_id` columns in `marts_index`, and an `fe_identity` table holding each
  technique's full configuration once.
- **H2O models** -- in `h2o_manifest.json`, as MLflow tags (`fe_version`,
  `fe_id`, `mart_id`, `cryptic_version`), and **inside every MOJO** as
  `mtb/fe_manifest.json`. H2O's loader ignores the extra entry (a MOJO carrying
  it predicts bit-identically on H2O 3.46.0.7), and it reads with a plain unzip:

      unzip -p model.zip mtb/fe_manifest.json

- **model releases** -- in each `model_card.json`, in the release `MANIFEST.json`
  with the full `fe_config`, and in the release name: a release is a
  (database, FE) pair, `mtb-resistotyper-ml-models-v3.4.0+fe1.0.3`.

## Artefacts built before the scheme

Marts and models built before FE was recorded carry no `fe` block. The packager
can reconstruct the configuration from the older sidecar fields and the source
tree at the revision that built them (`--fe-code-root`), and assign an FE version
only where the output has been verified equivalent (`--fe-version`). Where
neither is possible the identity is recorded as unknown, never guessed.
