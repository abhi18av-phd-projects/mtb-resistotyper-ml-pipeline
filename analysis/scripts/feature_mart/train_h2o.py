"""H2O AutoML -> Stacked Ensemble training on the RIF feature mart (design/03).

Completes the cycle: mart (build_mart.py) -> feature selection (causal.py)
-> MODEL (this). Produces the per-drug stacked-ensemble predictor — the
Paper-1 deliverable (see abc-universe
brainstorms/mtb-resistotyper-ml/2026-07-08-two-publication-split-idea.md).

Two feature sets, so the FE->metric dashboard gets real AUC deltas:
  denovo     — the full de-novo genomic set: every raw__ column (top-1000
               mutations + gene-level flags + quality summaries) + the
               lineage/coverage covariates. NO cat__ (catalogue leakage),
               NO cov__coresist_* (co-resistance is outcome-adjacent).
               This is the honest genotype->phenotype model.
  concordant — only the causal-concordant mutations from causal.py's
               results parquet (qf__causal_concordant) + the same
               covariates. Tests whether the tiny selected set matches the
               full set — the concordance-selection validation.

Honest CV: uses the mart's own __fold__ column as H2O's fold_column, so the
AutoML leaderboard + stacked-ensemble metrics are cross-validated on the
shared stratified folds (design/09 fold protocol).

Outputs per feature set, under <out>/<feature_set>/:
  leaderboard.csv                 the AutoML leaderboard
  SE_RIF_<feature_set>.zip        the leader MOJO (deployable, design/07)
  surrogate_terms.csv             surrogate-EBM top terms (global readout)
  manifest.json                   FE params + CV metrics (AUC/logloss/aucpr)
                                  -> tracking.py logs these to MLflow

Run (in the pixi env — it has h2o 3.46 + openjdk 21 + interpret):
    JAVA_HOME=$PWD/.pixi/envs/default/lib/jvm \
    pixi run python -m analysis.scripts.feature_mart.train_h2o \
        --mart analysis/results/feature_mart/feature_mart_RIF_cryptic-slim-2026.05_v1.0.2.parquet \
        --concordance analysis/results/causal/RIF_mutation/causal_concordance_RIF_mutation_level.parquet \
        --out analysis/results/h2o \
        --max-runtime-secs 180
"""

from __future__ import annotations

import argparse
import json
import os
import time
from pathlib import Path

import numpy as np
import pandas as pd

RANDOM_STATE = 42


def _feature_columns(df: pd.DataFrame, feature_set: str, concordant: list[str] | None) -> list[str]:
    """Build the model feature list for the chosen feature set.

    Covariates common to both: lineage one-hots + numeric coverage. Excluded
    everywhere: id__/y__/qf__/__fold__ (bookkeeping/labels), cat__ (catalogue
    leakage), cov__coresist_* (outcome-adjacent MDR co-resistance), and the
    string covariates (country/dataset/sublineage/lineage_raw) to keep the
    frame numeric for a clean first model.
    """
    covariates = [
        "cov__lineage_L1", "cov__lineage_L2", "cov__lineage_L3", "cov__lineage_L4",
        "cov__median_coverage", "cov__tb_breadth",
    ]
    covariates = [c for c in covariates if c in df.columns]

    if feature_set == "denovo":
        raw = [c for c in df.columns if c.startswith("raw__")]
        return raw + covariates
    if feature_set == "concordant":
        if not concordant:
            raise SystemExit("concordant feature set requested but no concordant mutations passed")
        present = [c for c in concordant if c in df.columns]
        return present + covariates
    raise ValueError(f"unknown feature_set {feature_set!r}")


def _load_concordant(concordance_path: Path | None) -> list[str]:
    if concordance_path is None or not concordance_path.exists():
        return []
    res = pd.read_parquet(concordance_path)
    return list(res.loc[res["qf__causal_concordant"], "feature_column"])


def train_one(
    df: pd.DataFrame,
    feature_set: str,
    features: list[str],
    out_dir: Path,
    max_runtime_secs: int,
) -> dict:
    import h2o
    from h2o.automl import H2OAutoML

    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    model_cols = features + ["y__binary", "__fold__"]
    hf = h2o.H2OFrame(df[model_cols])
    hf["y__binary"] = hf["y__binary"].asfactor()
    hf["__fold__"] = hf["__fold__"].asfactor()

    aml = H2OAutoML(
        max_runtime_secs=max_runtime_secs,
        seed=RANDOM_STATE,
        keep_cross_validation_predictions=True,
        sort_metric="AUC",
    )
    aml.train(x=features, y="y__binary", training_frame=hf, fold_column="__fold__")

    leader = aml.leader
    lb = aml.leaderboard.as_data_frame()
    lb.to_csv(out_dir / "leaderboard.csv", index=False)

    # Cross-validated metrics of the leader (honest — computed on the shared folds).
    perf = leader.model_performance(xval=True)
    metrics = {
        "cv_auc": float(perf.auc()),
        "cv_aucpr": float(perf.aucpr()) if perf.aucpr() is not None else None,
        "cv_logloss": float(perf.logloss()),
        "leader_model_id": leader.model_id,
        "leader_is_stacked_ensemble": "StackedEnsemble" in leader.model_id,
        "n_leaderboard_models": int(lb.shape[0]),
    }

    # MOJO export (deployable artifact, design/07).
    mojo_path = leader.download_mojo(path=str(out_dir), get_genmodel_jar=False)
    mojo_name = Path(mojo_path).name

    # Surrogate-EBM global readout: approximate the leader with a glass-box
    # EBM on the training features -> leader probabilities (design/03 bridge).
    surrogate = _fit_surrogate(df, features, leader, hf, out_dir)

    manifest = {
        "drug": "RIF",
        "stage": "h2o-stacked-ensemble (design/03, Paper 1)",
        "feature_set": feature_set,
        "n_features": len(features),
        "max_runtime_secs": max_runtime_secs,
        "random_state": RANDOM_STATE,
        "cv_metrics": metrics,
        "surrogate": surrogate,
        "outputs": {
            "mojo": mojo_name,
            "leaderboard": "leaderboard.csv",
        },
        "wall_time_seconds": round(time.time() - t0, 2),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def _fit_surrogate(df, features, leader, hf, out_dir) -> dict:
    """Fit an EBM on features -> leader's predicted P(R); report fidelity +
    top terms. Global glass-box readout of what the SE uses (design/03)."""
    from interpret.glassbox import ExplainableBoostingRegressor

    pred = leader.predict(hf).as_data_frame()
    # H2O binary predict frame: columns [predict, p0, p1]; p1 = P(class '1').
    p1 = pred["p1"].to_numpy() if "p1" in pred.columns else pred.iloc[:, -1].to_numpy()

    X = df[features].to_numpy()
    # Pass feature_names so the surrogate's term names are the real column
    # names (raw__<gene>_<mut>, cov__*) — otherwise interpret falls back to
    # positional feature_0000 names and the "which genes" readout is useless.
    ebm = ExplainableBoostingRegressor(
        interactions=0, random_state=RANDOM_STATE, feature_names=features
    )
    ebm.fit(X, p1)
    fidelity = float(ebm.score(X, p1))  # R^2 of the surrogate approximation

    terms = sorted(zip(ebm.term_names_, ebm.term_importances()), key=lambda t: t[1], reverse=True)
    top = [(name, float(imp)) for name, imp in terms[:15]]
    pd.DataFrame(top, columns=["term", "importance"]).to_csv(out_dir / "surrogate_terms.csv", index=False)
    return {"fidelity_r2": fidelity, "top_terms": [t[0] for t in top[:8]]}


def main() -> None:
    parser = argparse.ArgumentParser(description="H2O stacked-ensemble training on the RIF mart (design/03).")
    parser.add_argument("--mart", type=Path, required=True)
    parser.add_argument("--concordance", type=Path, default=None,
                        help="causal.py results parquet, for the concordant feature set.")
    parser.add_argument("--out", type=Path, default=Path("analysis/results/h2o"))
    parser.add_argument("--feature-set", choices=["denovo", "concordant", "all"], default="all")
    parser.add_argument("--max-runtime-secs", type=int, default=180)
    parser.add_argument("--h2o-mem", default="4G")
    args = parser.parse_args()

    os.environ.setdefault("JAVA_HOME", os.path.abspath(".pixi/envs/default/lib/jvm"))
    import h2o
    h2o.init(nthreads=-1, max_mem_size=args.h2o_mem, name=f"mtb-h2o-{os.getpid()}")

    df = pd.read_parquet(args.mart)
    concordant = _load_concordant(args.concordance)

    sets = ["denovo", "concordant"] if args.feature_set == "all" else [args.feature_set]
    results = {}
    for fs in sets:
        if fs == "concordant" and not concordant:
            print(f"skipping concordant set — no --concordance results provided")
            continue
        features = _feature_columns(df, fs, concordant)
        print(f"=== training feature_set={fs} ({len(features)} features) ===")
        results[fs] = train_one(df, fs, features, args.out / fs, args.max_runtime_secs)

    h2o.cluster().shutdown()

    print(json.dumps({fs: {"cv_auc": m["cv_metrics"]["cv_auc"], "n_features": m["n_features"]}
                      for fs, m in results.items()}, indent=2))


if __name__ == "__main__":
    main()
