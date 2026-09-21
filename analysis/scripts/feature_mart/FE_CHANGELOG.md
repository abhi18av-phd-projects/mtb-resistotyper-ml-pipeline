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

## 1.0.4 -- 2026-09-11

Causal-concordance selection made actually deterministic under its seed, not
just declared so.

- `causal.py` pins `n_jobs=1` on `CausalForestDML` and `ExplainableBoosting-
  Classifier` (previously the unset defaults, `-1` and `-2`), and pins
  `OMP_NUM_THREADS`/`OPENBLAS_NUM_THREADS`/`MKL_NUM_THREADS`/
  `NUMEXPR_NUM_THREADS`/`VECLIB_MAXIMUM_THREADS` to `1` before numpy/scipy are
  imported. Every voter's own randomness (per-tree seeds, subsample indices,
  cross-fit folds) was already seeded and generated sequentially before any
  parallel work starts; what was not pinned was thread count, and
  CausalForestDML's first-stage nuisance-model selection (`model_y='auto'`/
  `model_t='auto'`, choosing between 'forest' and 'linear' by cross-validated
  score) can land on either side of a near-tied score depending on it. A
  same-CPU-count rerun already matched to 1e-8 (float noise, same family every
  time); a different-CPU-count rerun previously moved the CATE estimate by
  whole units -- a different family winning, not noise.
- Does not touch `build_mart.py`; marts are unaffected.

PATCH, verified: not by rebuilding a mart (this change is in the selection
stage, not the mart builder) but the natural equivalent -- the same
selection stage, on the same mart (RIF, lineage1 held out), run three times,
two at 8 CPUs and one at 2. All three `causal_concordance_*.parquet` outputs
are byte-identical (same MD5), not merely close: every vote, every logged
statistic, to `0.000e+00`.

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
