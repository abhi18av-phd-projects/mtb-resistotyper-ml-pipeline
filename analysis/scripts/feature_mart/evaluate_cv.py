"""Honest CV evaluation across three regimes (design/03 §Honest evaluation).

The load-bearing experiment (design/03 §"Why lineage-aware CV is
load-bearing"): a random-k-fold AUC on clonal MTB is structurally
optimistic because near-identical isolates (same lineage / transmission
chain) land in both train and test. This script measures how much the AUC
drops when that leakage is removed — and whether the causally-concordant
feature set (docs/15) holds up better than the full de-novo set, which is
the decisive test for the concordance value proposition.

Fixed model across all conditions (LightGBM, seeded) so the ONLY thing
varying is the CV regime × feature set — this isolates the regime effect
(the scientific question), not model architecture. Absolute AUC will differ
from the H2O stacked ensemble (train_h2o.py); the regime *delta* is the
result. A full H2O-SE lineage-CV is a heavier follow-up.

Three regimes:
  random           StratifiedKFold(5) on y — the optimistic baseline
  lineage_loo      leave-one-lineage-out for L1-L4 — the Hunt/Iqbal test;
                   the held-out lineage's AUC is genomic generalisation to
                   an unseen population
  country_loo      leave-one-country-out — worst-case generalisation.
                   COVERAGE-LIMITED: ~85% of isolates have no country
                   (cov__country=''); only Viet Nam/Peru/South Africa/
                   Germany have data, so this is supplementary, not primary.

Two feature sets (matching train_h2o.py, minus the lineage one-hots which
are excluded everywhere because lineage is the CV grouping variable — using
it as a feature would be circular):
  denovo      all raw__ genomic + coverage covariates
  concordant  only qf__causal_concordant mutations + coverage covariates

Reports overall (pooled) AUC + per-lineage AUC per (feature_set, regime).
Manifest -> tracking.py -> MLflow.

Run (causal venv has lightgbm + sklearn; no H2O/Java needed):
    python -m analysis.scripts.feature_mart.evaluate_cv \
        --mart analysis/results/feature_mart/feature_mart_RIF_cryptic-slim-2026.05_v1.0.2.parquet \
        --concordance analysis/results/causal/RIF_mutation/causal_concordance_RIF_mutation_level.parquet \
        --out analysis/results/cv_eval
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COUNTRIES = ["Viet Nam", "Peru", "South Africa", "Germany"]
COVARIATES = ["cov__median_coverage", "cov__tb_breadth"]


def _lgbm() -> LGBMClassifier:
    return LGBMClassifier(
        n_estimators=300, learning_rate=0.05, num_leaves=31,
        subsample=0.8, colsample_bytree=0.6, random_state=RANDOM_STATE,
        n_jobs=-1, verbose=-1,
    )


def _feature_columns(df: pd.DataFrame, feature_set: str, concordant: list[str]) -> list[str]:
    cov = [c for c in COVARIATES if c in df.columns]
    if feature_set == "denovo":
        raw = [c for c in df.columns if c.startswith("raw__")]
        return raw + cov
    if feature_set == "concordant":
        return [c for c in concordant if c in df.columns] + cov
    raise ValueError(feature_set)


def _safe_auc(y_true: np.ndarray, y_score: np.ndarray) -> float | None:
    """AUC only if both classes are present (else undefined)."""
    if len(np.unique(y_true)) < 2:
        return None
    return float(roc_auc_score(y_true, y_score))


def _per_lineage_auc(y: np.ndarray, oof: np.ndarray, lineage: np.ndarray) -> dict:
    out = {}
    for lin in LINEAGES:
        mask = lineage == lin
        if mask.sum() >= 20:
            out[lin] = {"n": int(mask.sum()), "n_R": int(y[mask].sum()), "auc": _safe_auc(y[mask], oof[mask])}
    return out


def _random_regime(X: np.ndarray, y: np.ndarray, lineage: np.ndarray) -> dict:
    """StratifiedKFold(5); pooled OOF AUC + per-lineage AUC on OOF preds."""
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    oof = np.zeros(len(y))
    for tr, te in skf.split(X, y):
        m = _lgbm().fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    return {"overall_auc": _safe_auc(y, oof), "per_lineage": _per_lineage_auc(y, oof, lineage)}


def _leave_one_group_out(X, y, groups, group_values, lineage) -> dict:
    """Train on all-but-one group, test on the held-out group. The held-out
    group's AUC is generalisation to an unseen population."""
    per_group = {}
    pooled_true, pooled_score, pooled_lin = [], [], []
    for g in group_values:
        te = groups == g
        tr = ~te & (groups != "")  # train on labelled-other isolates
        # For lineage LOO, train on the other lineages (drop unlabeled/minor).
        if set(group_values) == set(LINEAGES):
            tr = np.isin(groups, [x for x in LINEAGES if x != g])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            per_group[g] = {"n": int(te.sum()), "n_R": int(y[te].sum()), "auc": None}
            continue
        m = _lgbm().fit(X[tr], y[tr])
        score = m.predict_proba(X[te])[:, 1]
        per_group[g] = {"n": int(te.sum()), "n_R": int(y[te].sum()), "auc": _safe_auc(y[te], score)}
        pooled_true.append(y[te]); pooled_score.append(score); pooled_lin.append(lineage[te])
    aucs = [v["auc"] for v in per_group.values() if v["auc"] is not None]
    pooled = {
        "overall_auc": _safe_auc(np.concatenate(pooled_true), np.concatenate(pooled_score)) if pooled_true else None,
        "mean_group_auc": float(np.mean(aucs)) if aucs else None,
        "per_group": per_group,
    }
    return pooled


def _tracked_metrics(manifest: dict) -> dict[str, float]:
    """Flatten a manifest into MLflow metric keys.

    Keys are namespaced with `/` because MLflow groups charts on the prefix, so
    `lineage_loo/*` becomes one section instead of five unrelated single-bar
    charts. Per-lineage AUCs are lifted onto the PARENT rather than left in the
    manifest blob: four bars in one chart is the comparison this evaluation
    exists to make, and it is unreadable if each lineage needs a drill-down.
    """
    out: dict[str, float] = {}
    res = manifest.get("results", {})

    for fs, byregime in res.items():
        if isinstance(byregime.get("n_features"), int):
            out[f"n_features/{fs}"] = float(byregime["n_features"])
        for regime in ("random", "lineage_loo", "country_loo"):
            r = byregime.get(regime) or {}
            for key, name in (("overall_auc", "overall"), ("mean_group_auc", "mean_group")):
                if isinstance(r.get(key), (int, float)):
                    out[f"{regime}/{fs}_{name}_auc"] = float(r[key])
            for group, g in (r.get("per_group") or {}).items():
                if isinstance(g.get("auc"), (int, float)):
                    out[f"{regime}/{fs}/{group}_auc"] = float(g["auc"])
                    out[f"support/{regime}_{fs}_{group}_n"] = float(g.get("n") or 0)

    # The two numbers this script exists to produce, computed once here rather
    # than eyeballed off two charts.
    #
    # leakage_delta: how much AUC a random k-fold buys you that an unseen
    # lineage does not. It is the optimism a random split hides, and a headline
    # AUC without it is not interpretable.
    #
    # concordance_gain: whether the causally-concordant feature set survives the
    # honest regime BETTER than the de-novo set. That is the decisive test for
    # the concordance value proposition, and it is a difference of differences,
    # so no single chart shows it.
    for fs in res:
        rand = (res[fs].get("random") or {}).get("overall_auc")
        loo = (res[fs].get("lineage_loo") or {}).get("mean_group_auc")
        if isinstance(rand, (int, float)) and isinstance(loo, (int, float)):
            out[f"headline/leakage_delta_{fs}"] = float(rand) - float(loo)

    dn = (res.get("denovo", {}).get("lineage_loo") or {}).get("mean_group_auc")
    co = (res.get("concordant", {}).get("lineage_loo") or {}).get("mean_group_auc")
    if isinstance(dn, (int, float)) and isinstance(co, (int, float)):
        out["headline/concordance_gain_lineage_loo"] = float(co) - float(dn)

    if isinstance(manifest.get("wall_time_seconds"), (int, float)):
        out["wall_time_seconds"] = float(manifest["wall_time_seconds"])
    return out


def _track(manifest: dict, out_dir: Path, held_out: str | None) -> None:
    """Project the manifest onto MLflow. Never fatal: the manifest is the truth.

    Without this the honest evaluation was invisible in MLflow -- the dashboard
    showed H2O's internal CV and nothing about generalisation to an unseen
    lineage, which is the claim the campaign is built to support.
    """
    try:
        from analysis.scripts.feature_mart import mlflow_tracking
    except Exception:
        return
    fs_label = "concordant" if "concordant" in manifest.get("results", {}) else "denovo"
    name = f"{manifest['drug']}-cv-{held_out}" if held_out else f"{manifest['drug']}-cv"
    with mlflow_tracking.run(name, fs_label, {
            "model": manifest.get("model"),
            "random_state": manifest.get("random_state"),
            "held_out_lineage": held_out or "none",
            "stage": "evaluate_cv"}) as h:
        h.metrics(_tracked_metrics(manifest))
        # Per (feature set, regime) as children, so the per-lineage support
        # counts stay browsable without crowding the parent.
        for fs, byregime in manifest.get("results", {}).items():
            for regime in ("random", "lineage_loo", "country_loo"):
                r = byregime.get(regime) or {}
                if not r:
                    continue
                h.child(f"{fs}:{regime}",
                        {"feature_set": fs, "regime": regime,
                         "n_features": byregime.get("n_features")},
                        {"overall_auc": r.get("overall_auc"),
                         "mean_group_auc": r.get("mean_group_auc")})
        h.artifact(out_dir)


def _infer_drug(mart_path: Path) -> str:
    """feature_mart_<DRUG>_cryptic-... -> <DRUG>."""
    name = mart_path.name
    if name.startswith("feature_mart_"):
        return name[len("feature_mart_"):].split("_")[0]
    return "UNKNOWN"


def evaluate(mart_path: Path, concordance_path: Path | None, out_dir: Path,
             drug: str | None = None, held_out: str | None = None) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    drug = drug or _infer_drug(mart_path)
    df = pd.read_parquet(mart_path)
    y = df["y__binary"].to_numpy()
    lineage = df["cov__lineage_raw"].to_numpy()
    country = df["cov__country"].fillna("").to_numpy()

    concordant = []
    if concordance_path and concordance_path.exists():
        res = pd.read_parquet(concordance_path)
        concordant = list(res.loc[res["qf__causal_concordant"], "feature_column"])

    results = {}
    for fs in ["denovo", "concordant"]:
        if fs == "concordant" and not concordant:
            continue
        feats = _feature_columns(df, fs, concordant)
        X = df[feats].to_numpy()
        results[fs] = {
            "n_features": len(feats),
            "random": _random_regime(X, y, lineage),
            "lineage_loo": _leave_one_group_out(X, y, lineage, LINEAGES, lineage),
            "country_loo": _leave_one_group_out(X, y, country, COUNTRIES, lineage),
        }

    manifest = {
        "drug": drug,
        "stage": "cv-evaluation (design/03 honest evaluation)",
        "model": "LightGBM (fixed across regimes to isolate the regime effect)",
        "random_state": RANDOM_STATE,
        "country_coverage_caveat": "~85% of isolates have no country; country_loo is supplementary",
        "results": results,
        "wall_time_seconds": round(time.time() - t0, 2),
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    _track(manifest, out_dir, held_out)
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description="Three-regime honest CV evaluation (design/03).")
    p.add_argument("--mart", type=Path, required=True)
    p.add_argument("--concordance", type=Path, default=None)
    p.add_argument("--out", type=Path, default=Path("analysis/results/cv_eval"))
    p.add_argument("--drug", default=None, help="Drug label (inferred from the mart filename if omitted).")
    p.add_argument("--held-out", default=None,
                   help="Held-out lineage this fold's feature selection was made without.")
    args = p.parse_args()

    m = evaluate(args.mart, args.concordance, args.out, drug=args.drug,
                 held_out=args.held_out)

    print("\n=== overall AUC by feature set × regime ===")
    print(f"{'feature_set':<12}{'n_feat':>7}{'random':>9}{'lineage_loo':>13}{'country_loo':>13}")
    for fs, r in m["results"].items():
        rnd = r["random"]["overall_auc"]
        lin = r["lineage_loo"]["mean_group_auc"]
        cty = r["country_loo"]["mean_group_auc"]
        print(f"{fs:<12}{r['n_features']:>7}{rnd:>9.4f}"
              f"{(lin if lin is not None else float('nan')):>13.4f}"
              f"{(cty if cty is not None else float('nan')):>13.4f}")


if __name__ == "__main__":
    main()
