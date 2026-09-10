# Feature-engineering changelog

The FE version (`fe_identity.FE_VERSION`) versions the technique that turns a
compendium database into the features a model sees: the mart builder, the
determinant-gene seeding, the optional FE steps and the concordance selection.
It is independent of the pipeline's release version and of the CRyPTIC release.

- MAJOR: features are encoded or named differently, or the fold scheme changes.
- MINOR: the feature space changes but the encoding does not.
- PATCH: output is unchanged on a fixed reference release (CRyPTIC v3.4.0).

`python -m analysis.scripts.feature_mart.fe_identity --check` fails when the FE
code has changed but the version has not. Bump, add an entry here, then `--lock`.

## 1.0.3 -- 2026-09-10

First version carried as `fe_version`, alongside the measured `fe_id`.

- Release label read from the database (`v3.4.0`, `v2.1.2`, `v3.4.0-slim`)
  instead of the hard-coded `slim-2026.05`. Metadata only.
- Optional genome and cohort covariates selected by intersection; the minor
  allele column mapped across releases (`IS_MINOR_ALLELE` -> `IS_MINOR`);
  numeric covariates a release does not publish are omitted, not NaN. Changes
  nothing on releases that publish every column.
- Concordance selection chooses confounders by intersection.

PATCH, verified: the v3.4.0 RIF mart built by this code and by the 1.0.2 code
is identical -- 14,490 x 1,817, same columns, same values.

## 1.0.2

Co-resistance confounder panel added (design/15). Carried until now only as the
mart filename label (`MART_VERSION_DEFAULT`), which every real run overrode with
`vFULL`; that label is a scope marker for file matching, not an FE version.
