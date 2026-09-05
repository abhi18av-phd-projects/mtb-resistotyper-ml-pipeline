"""Step 1 — deploy the usable/moderate per-drug models as tuned L2-logistic.

The model-comparison + DL-sweep results settled it: under lineage-LOO a regularised
linear model beats gradient-boosted trees and tuned deep nets, with tighter per-lineage
stability and full interpretability. This trains the production per-drug models on that
basis:

- **Honest C-tuning:** nested lineage-LOO — the outer loop holds out a lineage; the
  inner loop tunes `C` by lineage-LOO over the *training* lineages only, so `C` never
  sees the held-out lineage. The reported AUC is therefore a clean generalisation
  estimate (no tuning leakage).
- **Deployable artefact:** the final model is fit on all data with the C picked by
  full lineage-LOO. Because it is StandardScaler + L2-logistic, the whole model is a
  tiny JSON (feature list + scaler mean/scale + coefficients + intercept) — inference
  is pure numpy, no H2O/JVM. A **signed-coefficient catalogue** falls out directly
  (coefficient sign = resistance/susceptibility direction), which is the interpretable
  replacement for the GLM surrogate.

Usable + moderate tiers only (RIF/INH/EMB/MXF/LEV + ETH/KAN/AMI); the at-chance tier
(BDQ/CFZ/LZD/DLM) stays at the data floor and is not deployed.

In-memory, small JSON outputs (storage-constrained machine).

Run:
    python -m analysis.scripts.feature_mart.train_logistic \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --baseline analysis/results/fulldata_ml/model_comparison/manifest.json \
        --out analysis/results/fulldata_ml/logistic_models
"""

from __future__ import annotations

import argparse
import glob
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.preprocessing import StandardScaler

RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COV = ["cov__median_coverage", "cov__tb_breadth"]
DEPLOY = ["RIF", "INH", "EMB", "MXF", "LEV", "ETH", "KAN", "AMI"]  # usable + moderate
C_GRID = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0, 10.0]


def _fit(X, y, C):
    sc = StandardScaler().fit(X)
    lr = LogisticRegression(C=C, max_iter=2000).fit(sc.transform(X), y)
    return sc, lr


def _proba(sc, lr, X):
    return lr.predict_proba(sc.transform(X))[:, 1]


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _loo_auc_at_C(X, y, lineage, C, folds):
    """lineage-LOO mean AUC at a fixed C over the given lineage folds."""
    oof = np.full(len(y), np.nan)
    for L in folds:
        te = lineage == L
        tr = np.isin(lineage, [x for x in folds if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        sc, lr = _fit(X[tr], y[tr], C)
        oof[te] = _proba(sc, lr, X[te])
    per = {L: _auc(y[(lineage == L) & ~np.isnan(oof)], oof[(lineage == L) & ~np.isnan(oof)])
           for L in folds if ((lineage == L) & ~np.isnan(oof)).sum() >= 20}
    aucs = [v for v in per.values() if v is not None]
    return (np.mean(aucs) if aucs else None), per


def _best_C(X, y, lineage, folds):
    scores = {C: _loo_auc_at_C(X, y, lineage, C, folds)[0] for C in C_GRID}
    scores = {C: s for C, s in scores.items() if s is not None}
    return max(scores, key=scores.get) if scores else 1.0


def _nested_loo(X, y, lineage):
    """Outer lineage-LOO; inner lineage-LOO on training lineages tunes C (no leakage)."""
    oof = np.full(len(y), np.nan)
    chosen = {}
    for L in LINEAGES:
        te = lineage == L
        train_folds = [x for x in LINEAGES if x != L]
        tr = np.isin(lineage, train_folds)
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        C = _best_C(X[tr], y[tr], lineage[tr], train_folds)  # inner: train lineages only
        chosen[L] = C
        sc, lr = _fit(X[tr], y[tr], C)
        oof[te] = _proba(sc, lr, X[te])
    per = {L: _auc(y[(lineage == L) & ~np.isnan(oof)], oof[(lineage == L) & ~np.isnan(oof)])
           for L in LINEAGES if ((lineage == L) & ~np.isnan(oof)).sum() >= 20}
    aucs = [v for v in per.values() if v is not None]
    return {"honest_auc": round(float(np.mean(aucs)), 4) if aucs else None,
            "stability_std": round(float(np.std(aucs)), 4) if len(aucs) > 1 else None,
            "per_lineage": {k: (round(v, 3) if v else None) for k, v in per.items()},
            "C_per_outer_fold": chosen}


def train(marts_dir, causal_dir, baseline_path, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    base = {}
    if baseline_path and Path(baseline_path).exists():
        base = json.loads(Path(baseline_path).read_text()).get("per_drug", {})
    summary = {}
    for drug in DEPLOY:
        marts = glob.glob(str(marts_dir / f"feature_mart_{drug}_*vFULL.parquet"))
        cpaths = glob.glob(str(causal_dir / drug / f"causal_concordance_{drug}_mutation_level.parquet"))
        if not marts or not cpaths:
            continue
        df = pd.read_parquet(marts[0])
        y = df["y__binary"].to_numpy()
        lineage = df["cov__lineage_raw"].to_numpy()
        conc = [c for c in pd.read_parquet(cpaths[0]).query("qf__causal_concordant")["feature_column"] if c in df.columns]
        cov = [c for c in COV if c in df.columns]
        feats = conc + cov
        X = df[feats].to_numpy("float32")

        nested = _nested_loo(X, y, lineage)
        # deployable model: C by full lineage-LOO, fit on all data
        final_C = _best_C(X, y, lineage, LINEAGES)
        sc, lr = _fit(X, y, final_C)
        coef = lr.coef_[0]
        catalogue = sorted(
            [{"feature": f, "coef": round(float(c), 4),
              "direction": "R+" if c > 0 else "S-", "abs": abs(float(c)),
              "is_covariate": f in cov} for f, c in zip(feats, coef)],
            key=lambda r: -r["abs"])
        model = {"drug": drug, "model": "standardscaler+l2_logistic", "C": final_C,
                 "features": feats, "scaler_mean": sc.mean_.tolist(),
                 "scaler_scale": sc.scale_.tolist(), "coef": coef.tolist(),
                 "intercept": float(lr.intercept_[0]),
                 "predict": "p = sigmoid(((x - scaler_mean)/scaler_scale) . coef + intercept)"}
        (out_dir / f"{drug}_logistic.json").write_text(json.dumps(model, indent=2))

        b = base.get(drug, {})
        summary[drug] = {
            "n_features": len(feats), "n_concordant": len(conc), "tuned_C": final_C,
            "honest_nested_loo_auc": nested["honest_auc"], "stability_std": nested["stability_std"],
            "per_lineage": nested["per_lineage"], "C_per_outer_fold": nested["C_per_outer_fold"],
            "baseline_logistic_C1": b.get("logistic_l2", {}).get("auc"),
            "baseline_lightgbm": b.get("lightgbm", {}).get("auc"),
            "delta_vs_lightgbm": (round(nested["honest_auc"] - b["lightgbm"]["auc"], 4)
                                  if nested["honest_auc"] and b.get("lightgbm", {}).get("auc") else None),
            "top_signed_coefficients": catalogue[:10],
        }
        print(f"{drug:<5} C={final_C:<5} honest_loo={nested['honest_auc']} stab={nested['stability_std']} "
              f"| vs C=1 {summary[drug]['baseline_logistic_C1']} vs lgbm {summary[drug]['baseline_lightgbm']} "
              f"(Δlgbm {summary[drug]['delta_vs_lightgbm']:+})")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "step1 — deployed L2-logistic (nested-CV C-tuning) for usable/moderate tiers",
         "tiers_deployed": DEPLOY, "C_grid": C_GRID,
         "per_drug": summary, "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return summary


def main():
    p = argparse.ArgumentParser(description="Step 1: deploy usable/moderate models as tuned L2-logistic.")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--baseline", type=Path, default=Path("analysis/results/fulldata_ml/model_comparison/manifest.json"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/logistic_models"))
    a = p.parse_args()
    train(a.marts, a.causal, a.baseline, a.out)


if __name__ == "__main__":
    main()
