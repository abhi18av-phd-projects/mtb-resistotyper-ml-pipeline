"""How far can a tuned deep net (H2O-DL-class) get under HONEST lineage-LOO?

The DL literature (MD-CNN, RGTFormer, the 2026 Nat Commun CNN) reports 90–99%
accuracy — but on random/k-fold CV, the regime our study shows is inflated. The
question here: if we *optimise* a feedforward deep net — which is exactly what H2O's
Deep Learning estimator is (a multilayer perceptron with L1/L2 + dropout) — how close
does it get under lineage-LOO, and does more capacity help or hurt?

We sweep the same algorithm class as H2O-DL (sklearn `MLPClassifier`: feedforward net,
L2 `alpha`, early stopping — H2O-DL adds dropout/L1 but is the same family, so the
lineage-LOO conclusion transfers) over capacity × regularisation, on the *same*
concordant features as the model-comparison panel, and read it against the L2-logistic
and LightGBM baselines from `model_comparison/manifest.json`.

Hypothesis (from our own finding that linear > trees under lineage-LOO): a plain MLP
has *more* capacity than trees, so it should overfit lineage structure at least as
much and not beat the linear model — i.e. the DL papers' edge comes from features +
architecture (CNN locality, evolutionary/biophysical inputs, rare-variant
augmentation), not from depth per se.

In-memory, manifest only (storage-constrained machine).

Run:
    python -m analysis.scripts.feature_mart.evaluate_dl_sweep \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --baseline analysis/results/fulldata_ml/model_comparison/manifest.json \
        --out analysis/results/fulldata_ml/dl_sweep
"""

from __future__ import annotations

import argparse
import glob
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.exceptions import ConvergenceWarning
from sklearn.metrics import roc_auc_score
from sklearn.neural_network import MLPClassifier
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

warnings.filterwarnings("ignore", category=ConvergenceWarning)
RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COV = ["cov__median_coverage", "cov__tb_breadth"]

# capacity × L2 regularisation grid (H2O-DL analogue: hidden layers × l2)
ARCHS = {"mlp_32": (32,), "mlp_128": (128,), "mlp_64_32": (64, 32), "mlp_128_64_32": (128, 64, 32)}
ALPHAS = {"a1e-4": 1e-4, "a1e-2": 1e-2, "a1e0": 1.0}


def _mlp(hidden, alpha):
    return make_pipeline(StandardScaler(), MLPClassifier(
        hidden_layer_sizes=hidden, alpha=alpha, activation="relu", solver="adam",
        early_stopping=True, n_iter_no_change=15, max_iter=500, random_state=RANDOM_STATE))


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _loo(X, y, lineage, model_fn):
    oof = np.full(len(y), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        oof[te] = model_fn().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    per = {L: _auc(y[(lineage == L) & ~np.isnan(oof)], oof[(lineage == L) & ~np.isnan(oof)])
           for L in LINEAGES if ((lineage == L) & ~np.isnan(oof)).sum() >= 20}
    aucs = [v for v in per.values() if v is not None]
    return (round(float(np.mean(aucs)), 4) if aucs else None,
            round(float(np.std(aucs)), 4) if len(aucs) > 1 else None)


def evaluate(marts_dir, causal_dir, baseline_path, out_dir):
    t0 = time.time()
    base = {}
    if baseline_path and Path(baseline_path).exists():
        base = json.loads(Path(baseline_path).read_text()).get("per_drug", {})
    results = {}
    for cpath in sorted(glob.glob(str(causal_dir / "*" / "causal_concordance_*_mutation_level.parquet"))):
        drug = Path(cpath).name.split("_")[2]
        marts = glob.glob(str(marts_dir / f"feature_mart_{drug}_*vFULL.parquet"))
        if not marts:
            continue
        df = pd.read_parquet(marts[0])
        y = df["y__binary"].to_numpy()
        lineage = df["cov__lineage_raw"].to_numpy()
        conc = [c for c in pd.read_parquet(cpath).query("qf__causal_concordant")["feature_column"] if c in df.columns]
        cov = [c for c in COV if c in df.columns]
        if not conc:
            continue
        X = np.hstack([df[conc].to_numpy("float32"), df[cov].to_numpy("float32")])
        grid = {}
        best = (None, -1, None)
        for aname, arch in ARCHS.items():
            for lname, alpha in ALPHAS.items():
                auc, stab = _loo(X, y, lineage, lambda a=arch, al=alpha: _mlp(a, al))
                grid[f"{aname}/{lname}"] = {"auc": auc, "stab": stab}
                if auc is not None and auc > best[1]:
                    best = (f"{aname}/{lname}", auc, stab)
        b = base.get(drug, {})
        results[drug] = {
            "n_concordant": len(conc),
            "dl_best_config": best[0], "dl_best_auc": best[1] if best[1] >= 0 else None,
            "dl_best_stab": best[2],
            "baseline_logistic": b.get("logistic_l2", {}).get("auc"),
            "baseline_lightgbm": b.get("lightgbm", {}).get("auc"),
            "dl_minus_logistic": (round(best[1] - b["logistic_l2"]["auc"], 4)
                                  if best[1] >= 0 and b.get("logistic_l2", {}).get("auc") else None),
            "grid": grid,
        }
        print(f"{drug:<5} conc={len(conc):>3} | DL best={best[0]} {best[1]:.4f} (stab {best[2]}) "
              f"| logit={results[drug]['baseline_logistic']} lgbm={results[drug]['baseline_lightgbm']} "
              f"| DL-logit={results[drug]['dl_minus_logistic']:+}")
    manifest = {"experiment": "dl_sweep — tuned MLP (H2O-DL class) vs logistic/lgbm under lineage-LOO",
                "note": "MLPClassifier is the same algorithm family as H2O Deep Learning; "
                        "best-of-grid is the optimistic DL number.",
                "grid": {"archs": list(ARCHS), "alphas": list(ALPHAS)},
                "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}
    d = out_dir
    d.mkdir(parents=True, exist_ok=True)
    (d / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def main():
    p = argparse.ArgumentParser(description="DL (H2O-DL-class MLP) sweep under lineage-LOO.")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--baseline", type=Path, default=Path("analysis/results/fulldata_ml/model_comparison/manifest.json"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/dl_sweep"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.baseline, a.out)


if __name__ == "__main__":
    main()
