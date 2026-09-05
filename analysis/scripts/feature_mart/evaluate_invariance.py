"""Invariance (C5) — does restricting to lineage-invariant features generalise
more stably? (design/18 Track 4)

The honest-generalisation thesis, at the feature level: a *causal* mutation's
association with resistance is stable across environments (lineages); a
lineage-confounded hitchhiker's association appears only in the lineage where it
co-occurs with R. Invariant Causal Prediction formalises this. Here, per drug:

1. For each concordant feature, compute its within-lineage point-biserial
   correlation with y in each of L1–L4 (where it has ≥ MIN_CARR carriers).
2. **Invariance vote:** the effect is present (|r| ≥ R_MIN) with the SAME sign in
   ≥ MIN_LINEAGES of the lineages where it's computable → invariant; else
   lineage-specific (a co-occurrence artifact).
3. Re-run lineage-LOO (fixed LightGBM) with (a) all concordant features vs
   (b) the invariant-only subset, and report mean AUC AND per-lineage **stability**
   (std across held-out lineages). The bet: invariant-only is more stable and
   holds/gains AUC on the drugs whose concordant set carries lineage-specific
   noise (e.g. BDQ, whose transfer AUC collapsed on L2).

Run:
    python -m analysis.scripts.feature_mart.evaluate_invariance \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --out analysis/results/invariance_eval
"""

from __future__ import annotations

import argparse
import glob
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.metrics import roc_auc_score

RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COV = ["cov__median_coverage", "cov__tb_breadth"]
MIN_CARR = 10       # min carriers of a feature within a lineage to judge its effect
R_MIN = 0.10        # |point-biserial| threshold for "effect present"
MIN_LINEAGES = 3    # same-sign present in ≥ this many lineages → invariant


def _lgbm():
    return LGBMClassifier(n_estimators=300, learning_rate=0.05, num_leaves=31,
                          subsample=0.8, colsample_bytree=0.6, random_state=RANDOM_STATE,
                          n_jobs=-1, verbose=-1)


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _lineage_loo(df, feats, y, lineage):
    """Mean + per-lineage AUC and stability (std) for a feature set."""
    oof = np.full(len(df), np.nan)
    X = df[feats].to_numpy()
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        m = _lgbm().fit(X[tr], y[tr])
        oof[te] = m.predict_proba(X[te])[:, 1]
    per = {}
    for L in LINEAGES:
        mask = (lineage == L) & ~np.isnan(oof)
        if mask.sum() >= 20:
            per[L] = _auc(y[mask], oof[mask])
    aucs = [v for v in per.values() if v is not None]
    return {"mean_auc": round(float(np.mean(aucs)), 4) if aucs else None,
            "stability_std": round(float(np.std(aucs)), 4) if len(aucs) > 1 else None,
            "per_lineage": {k: (round(v, 3) if v else None) for k, v in per.items()}}


def _invariance(df, feat, y, lineage) -> dict:
    """Per-lineage point-biserial corr; sign-consistency → invariant flag."""
    per = {}
    for L in LINEAGES:
        m = lineage == L
        col = df[feat].to_numpy()[m]
        if col.sum() >= MIN_CARR and len(np.unique(y[m])) == 2 and col.std() > 0:
            per[L] = float(np.corrcoef(col, y[m])[0, 1])
    present = {L: r for L, r in per.items() if abs(r) >= R_MIN}
    signs = {np.sign(r) for r in present.values()}
    invariant = len(present) >= MIN_LINEAGES and len(signs) == 1
    return {"invariant": bool(invariant), "n_lineages_present": len(present),
            "per_lineage_r": {k: round(v, 3) for k, v in per.items()}}


def evaluate(marts_dir: Path, causal_dir: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
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
        inv = {f: _invariance(df, f, y, lineage) for f in conc}
        invariant = [f for f in conc if inv[f]["invariant"]]
        all_loo = _lineage_loo(df, conc + cov, y, lineage)
        inv_loo = _lineage_loo(df, invariant + cov, y, lineage) if invariant else None
        results[drug] = {
            "n_concordant": len(conc), "n_invariant": len(invariant),
            "invariant_features": invariant,
            "lineage_loo_all_concordant": all_loo,
            "lineage_loo_invariant_only": inv_loo,
        }
        print(f"{drug:<5} conc={len(conc):>3} inv={len(invariant):>3} | "
              f"all AUC={all_loo['mean_auc']} std={all_loo['stability_std']} | "
              f"inv AUC={inv_loo['mean_auc'] if inv_loo else '-'} std={inv_loo['stability_std'] if inv_loo else '-'}")
    manifest = {"experiment": "invariance (C5) — invariant-only vs all-concordant, lineage-LOO",
                "params": {"R_MIN": R_MIN, "MIN_LINEAGES": MIN_LINEAGES, "MIN_CARR": MIN_CARR},
                "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description="Invariance (C5) feature evaluation (design/18).")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/invariance_eval"))
    args = p.parse_args()
    evaluate(args.marts, args.causal, args.out)


if __name__ == "__main__":
    main()
