"""Three-layer H2O modeling — AutoML discovery -> curated base models ->
hand-built stacked ensembles (design/03, design/15; Paper 1 deliverable).

Rationale (why this supersedes train_h2o.py's AutoML-black-box): AutoML alone
hands back a single leader and discards the leaderboard, the individual base
models, and all control over the metalearner. Here the three layers are
explicit:

  LAYER 1 — AutoML as DISCOVERY. Run it to read the leaderboard (which
            families clear the bar, the ceiling AUC) and to keep AutoML's own
            stacked ensembles as a BASELINE TO BEAT. Not the deliverable.

  LAYER 2 — a curated portfolio of INDIVIDUAL base models, each trained
            deliberately and interpreted on its own:
              GLM   elastic-net logistic (lambda-search) -> signed coefficient
                    catalogue over mutations, directly comparable to WHO-2023.5
              GBM   tuned gradient boosting (usual top performer)
              XGB   second boosting flavor (diversity; guarded — H2O XGBoost is
                    frequently unavailable on macOS, skipped if so)
              DRF   distributed random forest (bagged trees)
              XRT   extremely randomized trees (histogram_type=Random)
              DL    feedforward net (diversity; slow, uninterpretable)

  LAYER 3 — stacked ensembles WE build over the Layer-2 models via
            H2OStackedEnsembleEstimator, with the metalearner treated as an
            experiment (GLM metalearner -> readable per-base-model weights;
            GBM metalearner -> nonlinear blend). We read se.metalearner() so
            "the stack trusts GBM 0.6 / GLM 0.25 / DRF 0.15" is a REPORTED
            result, not a hidden default.

THE UNIFYING INVARIANT — honest CV and valid stacking are the SAME constraint:
H2O stacking is only leakage-free if the base models' cross-validation
predictions are leakage-free, and the metalearner trains on exactly those
level-1 predictions. So every base model AND the metalearner use the mart's
`__fold__` (GroupKFold-by-isolate / shared stratified folds) as fold_column.
The ensemble is therefore honest BY CONSTRUCTION, not by a bolt-on afterward.

Two feature sets (as in train_h2o.py): `denovo` (all raw__ + covariates, no
catalogue/co-resistance leakage) and `concordant` (causal.py's
qf__causal_concordant mutations + covariates) — so the concordance value
proposition is measured per-algorithm AND at the ensemble level.

Outputs per feature set under <out>/<feature_set>/:
  automl_leaderboard.csv       Layer-1 leaderboard
  base_models.csv              Layer-2 per-family xval AUC / AUCPR
  glm_catalogue.csv            Layer-2 GLM nonzero signed coefficients
  varimp_<family>.csv          Layer-2 tree/DL variable importance
  ensembles.csv                Layer-3 SE xval AUC + metalearner weights
  SE_<drug>_<fs>_best.zip       best model's MOJO (deployable, design/07)
  surrogate_terms.csv          surrogate-EBM readout of the best model
  manifest.json                everything above -> tracking.py -> MLflow

Run (pixi env: h2o 3.46 + openjdk 21 + interpret):
    JAVA_HOME=$PWD/.pixi/envs/default/lib/jvm \
    pixi run python -m analysis.scripts.feature_mart.train_h2o_stacked \
        --mart analysis/results/feature_mart/feature_mart_RIF_cryptic-slim-2026.05_v1.0.2.parquet \
        --concordance analysis/results/causal/RIF_mutation/causal_concordance_RIF_mutation_level.parquet \
        --out analysis/results/h2o_stacked --automl-secs 120
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from pathlib import Path

import pandas as pd

# Local, and deliberately importable without mlflow installed: the module
# degrades to no-ops so the untracked and tracked code paths cannot diverge.
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
from analysis.scripts.feature_mart import fe_identity, mlflow_tracking

RANDOM_STATE = 42
FOLD_COL = "__fold__"


# --------------------------------------------------------------------------- #
# feature-set plumbing (shared shape with train_h2o.py)
# --------------------------------------------------------------------------- #
LINEAGE_COV = ["cov__lineage_L1", "cov__lineage_L2", "cov__lineage_L3", "cov__lineage_L4"]


def _covariates(df: pd.DataFrame, drop_lineage: bool = False) -> list[str]:
    """Model covariates, optionally without the lineage one-hots.

    They must be dropped whenever a whole lineage is held out. Within such a
    split the one-hot is constant on each side, so it identifies the held-out
    group exactly: the model is handed the grouping variable as a feature, which
    is the circularity evaluate_cv.py excludes for the same reason.
    """
    cov = LINEAGE_COV + ["cov__median_coverage", "cov__tb_breadth"]
    if drop_lineage:
        cov = [c for c in cov if c not in LINEAGE_COV]
    return [c for c in cov if c in df.columns]


def _feature_columns(df: pd.DataFrame, feature_set: str, concordant: list[str] | None,
                     drop_lineage: bool = False) -> list[str]:
    cov = _covariates(df, drop_lineage=drop_lineage)
    if feature_set == "denovo":
        return [c for c in df.columns if c.startswith("raw__")] + cov
    if feature_set == "concordant":
        if not concordant:
            raise SystemExit("concordant feature set requested but no concordant mutations passed")
        return [c for c in concordant if c in df.columns] + cov
    raise ValueError(f"unknown feature_set {feature_set!r}")


def _load_concordant(concordance_path: Path | None) -> list[str]:
    if concordance_path is None or not concordance_path.exists():
        return []
    res = pd.read_parquet(concordance_path)
    return list(res.loc[res["qf__causal_concordant"], "feature_column"])


def _infer_drug(mart_path: Path) -> str:
    name = mart_path.name
    if name.startswith("feature_mart_"):
        return name[len("feature_mart_"):].split("_")[0]
    return "UNKNOWN"


# --------------------------------------------------------------------------- #
# LAYER 1 — AutoML discovery
# --------------------------------------------------------------------------- #
def _layer1_automl(hf, features, y, sort_metric, secs, out_dir) -> dict:
    """Run AutoML to map the achievable ceiling + keep its SEs as a baseline."""
    from h2o.automl import H2OAutoML

    aml = H2OAutoML(
        max_runtime_secs=secs, seed=RANDOM_STATE,
        keep_cross_validation_predictions=True, sort_metric=sort_metric,
    )
    aml.train(x=features, y=y, training_frame=hf, fold_column=FOLD_COL)
    lb = aml.leaderboard.as_data_frame()
    lb.to_csv(out_dir / "automl_leaderboard.csv", index=False)

    automl_se = lb[lb["model_id"].str.contains("StackedEnsemble")]
    return {
        "n_models": int(lb.shape[0]),
        "leader_model_id": aml.leader.model_id,
        "leader_auc": float(lb.iloc[0]["auc"]) if "auc" in lb.columns else None,
        "best_automl_se_auc": float(automl_se["auc"].max()) if not automl_se.empty else None,
        "families_seen": sorted({m.split("_")[0] for m in lb["model_id"]}),
    }


# --------------------------------------------------------------------------- #
# LAYER 2 — curated individual base models
# --------------------------------------------------------------------------- #
def _base_model_specs(balance_classes: bool):
    """One deliberate config per family. Common kwargs (fold_column,
    keep_cross_validation_predictions, seed) are injected in _train_base so
    every model is stack-compatible AND honest on the shared folds."""
    from h2o.estimators import (
        H2OGeneralizedLinearEstimator, H2OGradientBoostingEstimator,
        H2ORandomForestEstimator, H2ODeepLearningEstimator,
    )
    bc = {"balance_classes": balance_classes}
    specs = [
        ("GLM", H2OGeneralizedLinearEstimator, dict(
            family="binomial", alpha=0.5, lambda_search=True, standardize=True)),
        ("GBM", H2OGradientBoostingEstimator, dict(
            ntrees=250, max_depth=6, learn_rate=0.05, sample_rate=0.8,
            col_sample_rate=0.8, **bc)),
        ("DRF", H2ORandomForestEstimator, dict(ntrees=250, max_depth=20, **bc)),
        ("XRT", H2ORandomForestEstimator, dict(
            ntrees=250, max_depth=20, histogram_type="Random", **bc)),
        ("DL", H2ODeepLearningEstimator, dict(
            hidden=[64, 64], epochs=50, activation="RectifierWithDropout")),
    ]
    return specs


def _maybe_xgboost(balance_classes: bool):
    """H2O XGBoost is frequently unavailable on macOS — guard it."""
    try:
        from h2o.estimators import H2OXGBoostEstimator
        if not H2OXGBoostEstimator.available():
            return None
        return ("XGB", H2OXGBoostEstimator, dict(
            ntrees=250, max_depth=6, learn_rate=0.05, subsample=0.8,
            col_sample_rate=0.8))
    except Exception:
        return None


def _train_base(hf, features, y, specs, out_dir) -> list[dict]:
    models = []
    for family, Est, params in specs:
        t = time.time()
        # fold_column is passed to train() (not the ctor); keep_cross_validation
        # _predictions=True retains the level-1 frame the metalearner needs.
        est = Est(
            model_id=f"base_{family}_{os.getpid()}",
            keep_cross_validation_predictions=True,
            seed=RANDOM_STATE,
            **params,
        )
        est.train(x=features, y=y, training_frame=hf, fold_column=FOLD_COL)
        perf = est.model_performance(xval=True)
        rec = {
            "family": family,
            "model": est,
            "model_id": est.model_id,
            "cv_auc": float(perf.auc()),
            "cv_aucpr": float(perf.aucpr()) if perf.aucpr() is not None else None,
            "cv_logloss": float(perf.logloss()),
            "train_secs": round(time.time() - t, 1),
        }
        _dump_interpretability(family, est, out_dir)
        models.append(rec)
        print(f"  base {family:<4} xval AUC={rec['cv_auc']:.4f} "
              f"AUCPR={rec['cv_aucpr'] and round(rec['cv_aucpr'],4)} ({rec['train_secs']}s)")
    return models


def _dump_interpretability(family: str, est, out_dir: Path) -> None:
    """GLM -> signed coefficient catalogue; tree/DL -> variable importance."""
    try:
        if family == "GLM":
            coef = est.coef()  # {feature: coefficient} on the standardized scale
            rows = [(k, float(v)) for k, v in coef.items() if k != "Intercept" and abs(v) > 1e-8]
            rows.sort(key=lambda r: abs(r[1]), reverse=True)
            pd.DataFrame(rows, columns=["feature", "coefficient"]).to_csv(
                out_dir / "glm_catalogue.csv", index=False)
        else:
            vi = est.varimp(use_pandas=True)
            if vi is not None:
                vi.head(30).to_csv(out_dir / f"varimp_{family}.csv", index=False)
    except Exception as e:  # interpretability is best-effort, never fatal
        print(f"    [warn] interpretability dump for {family} failed: {e}")


# --------------------------------------------------------------------------- #
# LAYER 3 — hand-built stacked ensembles + metalearner experiment
# --------------------------------------------------------------------------- #
def _metalearner_weights(se) -> dict | None:
    """Read the metalearner so base-model trust is a reported result."""
    try:
        ml = se.metalearner()
        if ml is None:
            return None
        # GLM metalearner -> coef dict keyed by base model_id.
        if hasattr(ml, "coef"):
            try:
                return {k: float(v) for k, v in ml.coef().items()}
            except Exception:
                pass
        vi = ml.varimp(use_pandas=True)
        if vi is not None:
            return {r["variable"]: float(r["relative_importance"]) for _, r in vi.iterrows()}
    except Exception:
        return None
    return None


def _train_ensembles(hf, features, y, base, out_dir) -> list[dict]:
    from h2o.estimators import H2OStackedEnsembleEstimator

    all_ids = [b["model_id"] for b in base]
    # pruned = top-4 base models by xval AUC (diversity/parsimony)
    topk = [b["model_id"] for b in sorted(base, key=lambda b: b["cv_auc"], reverse=True)[:4]]

    configs = [
        ("SE_all_glm", all_ids, "glm"),   # full portfolio, interpretable blend
        ("SE_all_gbm", all_ids, "gbm"),   # metalearner experiment: nonlinear blend
        ("SE_top4_glm", topk, "glm"),     # pruned portfolio
    ]
    out = []
    for name, ids, meta_algo in configs:
        t = time.time()
        se = H2OStackedEnsembleEstimator(
            model_id=f"{name}_{os.getpid()}",
            base_models=ids,
            metalearner_algorithm=meta_algo,
            # metalearner cross-validated on the SAME folds -> honest SE metric
            metalearner_fold_column=FOLD_COL,
            seed=RANDOM_STATE,
        )
        # x MUST be passed explicitly: with metalearner_fold_column set, H2O
        # removes the fold column from `ignored_columns`, which is only
        # populated when x is given (else it is None and .remove() crashes).
        se.train(x=features, y=y, training_frame=hf)
        perf = se.model_performance(xval=True)
        rec = {
            "name": name,
            "metalearner": meta_algo,
            "n_base_models": len(ids),
            "base_models": ids,
            "cv_auc": float(perf.auc()),
            "cv_aucpr": float(perf.aucpr()) if perf.aucpr() is not None else None,
            "metalearner_weights": _metalearner_weights(se),
            "model": se,
            "train_secs": round(time.time() - t, 1),
        }
        out.append(rec)
        print(f"  ens  {name:<12} xval AUC={rec['cv_auc']:.4f} (meta={meta_algo}, {rec['train_secs']}s)")
    return out


# --------------------------------------------------------------------------- #
# surrogate EBM on the winning model (glass-box global readout)
# --------------------------------------------------------------------------- #
def _fit_surrogate(df, features, model, hf, out_dir) -> dict:
    """Glass-box global readout of the winning model.

    The surrogate explains the model; it is not the model. Losing it should cost
    the explanation, not the fit. An unguarded import here destroyed ten
    completed trainings that had already been logged to the tracking server,
    because the image was built without interpret. The manifest records the
    absence so a reader can tell "no surrogate was fitted" from "the surrogate
    found nothing".
    """
    try:
        from interpret.glassbox import ExplainableBoostingRegressor
    except ImportError as exc:
        print(f"[surrogate] SKIPPED: {exc}. The models are unaffected; the "
              f"glass-box readout is absent from this run.", file=sys.stderr)
        return {"fitted": False, "reason": f"interpret unavailable: {exc}"}

    pred = model.predict(hf).as_data_frame()
    p1 = pred["p1"].to_numpy() if "p1" in pred.columns else pred.iloc[:, -1].to_numpy()
    X = df[features].to_numpy()
    ebm = ExplainableBoostingRegressor(interactions=0, random_state=RANDOM_STATE, feature_names=features)
    ebm.fit(X, p1)
    fidelity = float(ebm.score(X, p1))
    terms = sorted(zip(ebm.term_names_, ebm.term_importances()), key=lambda t: t[1], reverse=True)
    pd.DataFrame([(n, float(i)) for n, i in terms[:15]], columns=["term", "importance"]).to_csv(
        out_dir / "surrogate_terms.csv", index=False)
    return {"fidelity_r2": fidelity, "top_terms": [n for n, _ in terms[:8]]}


# --------------------------------------------------------------------------- #
# orchestration
# --------------------------------------------------------------------------- #
def _score_holdout(recs, hf_test, label):
    """AUC of each candidate on rows it never saw. Non-fatal: a scoring failure
    must not lose a trained model, so the key is simply absent."""
    for r in recs:
        try:
            r["holdout_auc"] = float(r["model"].model_performance(test_data=hf_test).auc())
        except Exception as exc:                      # noqa: BLE001
            print(f"  [warn] holdout scoring failed for {r.get(label, '?')}: {exc}")
    ranked = [r for r in recs if r.get("holdout_auc") is not None]
    for r in sorted(ranked, key=lambda r: -r["holdout_auc"]):
        print(f"  holdout {r.get(label, '?'):<14} AUC={r['holdout_auc']:.4f}")


def train_one(df, drug, feature_set, features, out_dir, automl_secs, sort_metric,
              balance_classes, holdout_df=None) -> dict:
    import h2o

    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    y = "y__binary"

    hf = h2o.H2OFrame(df[features + [y, FOLD_COL]])
    hf[y] = hf[y].asfactor()
    hf[FOLD_COL] = hf[FOLD_COL].asfactor()

    print(f"--- LAYER 1 (AutoML discovery, {automl_secs}s) ---")
    layer1 = _layer1_automl(hf, features, y, sort_metric, automl_secs, out_dir)

    print("--- LAYER 2 (curated base models) ---")
    specs = _base_model_specs(balance_classes)
    xgb = _maybe_xgboost(balance_classes)
    if xgb:
        specs.insert(2, xgb)
    else:
        print("  [info] H2O XGBoost unavailable on this platform — skipping XGB")
    base = _train_base(hf, features, y, specs, out_dir)
    pd.DataFrame([{k: b[k] for k in ("family", "model_id", "cv_auc", "cv_aucpr", "cv_logloss", "train_secs")}
                  for b in base]).to_csv(out_dir / "base_models.csv", index=False)

    print("--- LAYER 3 (hand-built stacked ensembles) ---")
    ens = _train_ensembles(hf, features, y, base, out_dir)
    holdout = {}
    if holdout_df is not None and len(holdout_df):
        # The comparison the leaderboard cannot make. cv_auc is out-of-fold within
        # the training lineages; this is the score on a lineage no fold contained,
        # which is the axis evaluate_cv.py reports on. Same rows, same features,
        # so an H2O candidate and the fixed baseline are finally comparable.
        print(f"--- HOLD-OUT ({len(holdout_df):,} rows the model never saw) ---")
        hf_test = h2o.H2OFrame(holdout_df[features + [y]])
        hf_test[y] = hf_test[y].asfactor()
        _score_holdout(base, hf_test, "family")
        _score_holdout(ens, hf_test, "name")
        scored = [r for r in base + ens if r.get("holdout_auc") is not None]
        if scored:
            best = max(scored, key=lambda r: r["holdout_auc"])
            holdout = {
                "n_rows": int(len(holdout_df)),
                "n_resistant": int(holdout_df[y].sum()),
                "best": best.get("family") or best.get("name"),
                "best_auc": best["holdout_auc"],
                "by_model": {(r.get("family") or r.get("name")): r["holdout_auc"] for r in scored},
            }

    pd.DataFrame([{k: e[k] for k in ("name", "metalearner", "n_base_models", "cv_auc", "cv_aucpr", "train_secs")}
                  for e in ens]).to_csv(out_dir / "ensembles.csv", index=False)

    # best model overall = best of {base models, our ensembles}
    best_base = max(base, key=lambda b: b["cv_auc"])
    best_ens = max(ens, key=lambda e: e["cv_auc"])
    best = best_ens if best_ens["cv_auc"] >= best_base["cv_auc"] else best_base
    stacking_lift = round(best_ens["cv_auc"] - best_base["cv_auc"], 4)

    mojo_path = best["model"].download_mojo(path=str(out_dir), get_genmodel_jar=False)
    surrogate = _fit_surrogate(df, features, best["model"], hf, out_dir)

    manifest = {
        "drug": drug,
        "stage": "h2o-3layer (AutoML discovery -> base models -> stacked ensembles)",
        "feature_set": feature_set,
        "n_features": len(features),
        "random_state": RANDOM_STATE,
        "sort_metric": sort_metric,
        "balance_classes": balance_classes,
        "layer1_automl": layer1,
        "layer2_base_models": [
            {k: b[k] for k in ("family", "model_id", "cv_auc", "cv_aucpr", "cv_logloss", "train_secs")}
            for b in base
        ],
        "layer3_ensembles": [
            {k: e[k] for k in ("name", "metalearner", "n_base_models", "base_models",
                               "cv_auc", "cv_aucpr", "metalearner_weights", "train_secs")}
            for e in ens
        ],
        "best_model": {
            "which": best["name"] if "name" in best else best["family"],
            "cv_auc": best["cv_auc"],
            "is_ensemble": "name" in best,
        },
        "stacking_lift_over_best_base": stacking_lift,
        "best_base_family": best_base["family"],
        "best_base_auc": best_base["cv_auc"],
        "holdout": holdout,
        "surrogate": surrogate,
        "outputs": {
            "mojo": Path(mojo_path).name,
            "automl_leaderboard": "automl_leaderboard.csv",
            "base_models": "base_models.csv",
            "ensembles": "ensembles.csv",
            "glm_catalogue": "glm_catalogue.csv",
        },
        "wall_time_seconds": round(time.time() - t0, 2),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def _embed_fe(out_dir: Path, fe: dict, drug: str, feature_set: str) -> None:
    """Write mtb/fe_manifest.json INSIDE every MOJO in out_dir.

    A MOJO is a zip, and H2O's loader ignores entries it does not know: on
    3.46.0.7 a MOJO carrying this file loads and predicts bit-identically to the
    original. The manifest therefore travels with the model file itself and is
    readable with a plain unzip, no H2O needed. Re-check on an H2O upgrade.
    """
    import zipfile
    body = json.dumps({"drug": drug, "feature_set": feature_set, **fe}, indent=2)
    for z in sorted(out_dir.glob("*.zip")):
        with zipfile.ZipFile(z, "a") as zf:
            if "mtb/fe_manifest.json" not in zf.namelist():
                zf.writestr("mtb/fe_manifest.json", body)


def main() -> None:
    p = argparse.ArgumentParser(description="Three-layer H2O modeling (AutoML -> base models -> stacked ensembles).")
    p.add_argument("--mart", type=Path, required=True)
    p.add_argument("--concordance", type=Path, default=None)
    p.add_argument("--out", type=Path, default=Path("analysis/results/h2o_stacked"))
    p.add_argument("--feature-set", choices=["denovo", "concordant", "all"], default="all")
    p.add_argument("--drug", default=None)
    # --max-runtime-secs is the name train_h2o.py uses and the name the Nextflow
    # module passes. Without the alias the module's default entrypoint dies at
    # argument parsing, which is a configuration error that only shows up once a
    # mart has already been built.
    p.add_argument("--automl-secs", "--max-runtime-secs", dest="automl_secs",
                   type=int, default=120)
    p.add_argument("--sort-metric", default="AUC", help="AUC (balanced) or AUCPR (low-prevalence drugs).")
    p.add_argument("--balance-classes", action="store_true",
                   help="enable for the ~1-5%%R drugs (BDQ/LZD/DLM/CFZ); off for RIF/INH.")
    p.add_argument("--h2o-mem", default="6G")
    p.add_argument("--mart-meta", type=Path, default=None,
                   help="the mart's metadata sidecar (default: beside --mart)")
    p.add_argument("--selection", type=Path, default=None,
                   help="the selection manifest that chose the concordant features")
    p.add_argument("--fe-pre-steps", default="", help="label-free FE steps that ran")
    p.add_argument("--fe-fold-steps", default="", help="label-using FE steps that ran")
    p.add_argument("--held-out-lineage", default=None,
                   help="train without this lineage and score every candidate on "
                        "it. Turns the run into the EVALUATION arm: comparable to "
                        "evaluate_cv.py because it is the same held-out rows. "
                        "Omit for the deployment refit on all data.")
    args = p.parse_args()

    os.environ.setdefault("JAVA_HOME", os.path.abspath(".pixi/envs/default/lib/jvm"))
    import h2o
    h2o.init(nthreads=-1, max_mem_size=args.h2o_mem, name=f"mtb-h2o-stack-{os.getpid()}")

    drug = args.drug or _infer_drug(args.mart)
    df = pd.read_parquet(args.mart)
    concordant = _load_concordant(args.concordance)

    # The FE identity of the features these models are trained on. Computed once
    # and carried into the manifest, the MLflow tags, and INSIDE every MOJO, so a
    # model file separated from everything else still says what FE made it. A
    # mart that predates the recorded FE block yields no identity; training goes
    # ahead and the manifest says the identity is unknown rather than guessing.
    fe = None
    meta_path = args.mart_meta or args.mart.with_suffix(".metadata.json")
    try:
        meta = json.loads(meta_path.read_text()) if meta_path.exists() else {}
        sel = json.loads(args.selection.read_text()) if args.selection and args.selection.exists() else {}
        fe = fe_identity.identity(meta, sel, args.fe_pre_steps, args.fe_fold_steps)
        fe["mart_id"] = fe_identity.mart_id(
            fe["fe_id"], (meta.get("source_db") or {}).get("checksum_sha256"))
        fe["cryptic_version"] = meta.get("cryptic_version")
        print(f"FE {fe['fe_version']} ({fe['fe_id']}), {fe['mart_id']}")
    except SystemExit as exc:
        print(f"no FE identity for this mart: {exc}")

    # Hold a whole lineage out, or refit on everything. 'none' is accepted so the
    # workflow can pass the fold value straight through without a conditional.
    held = (args.held_out_lineage or "").strip()
    holdout_df = None
    if held and held.lower() != "none":
        if "cov__lineage_raw" not in df.columns:
            raise SystemExit("--held-out-lineage needs cov__lineage_raw in the mart")
        mask = df["cov__lineage_raw"].astype(str) == held
        if not mask.any():
            raise SystemExit(f"--held-out-lineage {held!r} matches no rows")
        holdout_df, df = df[mask].copy(), df[~mask].copy()
        print(f"held out {held}: {len(holdout_df):,} rows withheld, {len(df):,} for training")

    sets = ["denovo", "concordant"] if args.feature_set == "all" else [args.feature_set]
    results = {}
    for fs in sets:
        if fs == "concordant" and not concordant:
            print("skipping concordant set — no --concordance results provided")
            continue
        feats = _feature_columns(df, fs, concordant, drop_lineage=bool(holdout_df is not None))
        print(f"\n=== {drug} / feature_set={fs} ({len(feats)} features) ===")
        # Tracking wraps training rather than following it: a run that dies
        # halfway is the one most worth inspecting, and a backfill from
        # manifest.json can only ever record the runs that finished. The handle
        # is a no-op when MLFLOW_TRACKING_URI is unset or the server is
        # unreachable, so an hour of cluster time is never lost to a tracking
        # failure.
        out_dir = args.out / fs
        fe_tags = {k: fe[k] for k in ("fe_version", "fe_id", "mart_id", "cryptic_version")} if fe else {}
        with mlflow_tracking.run(drug, fs, tags=fe_tags, params={
                "n_features": len(feats),
                "automl_secs": args.automl_secs,
                "sort_metric": args.sort_metric,
                "balance_classes": args.balance_classes,
                "h2o_mem": args.h2o_mem,
                "mart": args.mart.name,
                "random_state": RANDOM_STATE,
                "held_out_lineage": held or "none",
                "arm": "evaluation" if holdout_df is not None else "deployment",
        }) as tracked:
            results[fs] = train_one(df, drug, fs, feats, out_dir,
                                    args.automl_secs, args.sort_metric, args.balance_classes,
                                    holdout_df=holdout_df)
            # Before anything is logged or checksummed: the MOJO a tracking
            # server or a deposit receives is the one that carries its FE.
            results[fs]["fe"] = fe
            (out_dir / "manifest.json").write_text(json.dumps(results[fs], indent=2))
            if fe:
                _embed_fe(out_dir, fe, drug, fs)
            mlflow_tracking.log_manifest(tracked, results[fs], out_dir)

    h2o.cluster().shutdown()

    print("\n=== summary: best xval AUC by feature set ===")
    for fs, m in results.items():
        print(f"  {fs:<11} best={m['best_model']['which']:<12} AUC={m['best_model']['cv_auc']:.4f} "
              f"(stacking lift over best base {m['stacking_lift_over_best_base']:+.4f})")


if __name__ == "__main__":
    main()
