"""Multi-drug transfer under lineage-LOO — does a single SE beat per-drug models?

The transfer bet (design/15 Part B + T14 drug-class): train ONE model over the
long-format (isolate, drug) mart with `id__DRUG` + `cov__drug_class` as
categorical features, so same-class drugs share signal (FQ MXF/LEV,
aminoglycoside KAN/AMI; efflux BDQ↔CFZ via the shared Rv0678 feature). The
decisive test is whether this lifts the per-drug **lineage-LOO** AUC — especially
for the moderate/at-chance tiers — over the independent per-drug models.

Fixed LightGBM (as in evaluate_cv.py) to isolate the *transfer* effect, not model
architecture. Lineage-LOO: train on 3 lineages (all drugs), predict the held-out
lineage, then score per drug. Reported per drug = mean over held-out L1–L4,
alongside the per-drug baseline (cv_eval_fulldata) → the transfer delta.

Run:
    python -m analysis.scripts.feature_mart.evaluate_multidrug_cv \
        --mart <multidrug mart.parquet> --baseline analysis/results/cv_eval_fulldata/summary.json \
        --out ./md_eval
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

RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]


def _lgbm() -> LGBMClassifier:
    return LGBMClassifier(n_estimators=400, learning_rate=0.05, num_leaves=63,
                          subsample=0.8, colsample_bytree=0.6, random_state=RANDOM_STATE,
                          n_jobs=-1, verbose=-1)


def _safe_auc(y, p) -> float | None:
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def evaluate(mart_path: Path, baseline_path: Path | None, out_dir: Path,
             drugs: list[str] | None = None) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    df = pd.read_parquet(mart_path)
    if drugs:  # mechanism-grouped ablation: restrict the pool to these drugs only
        df = df[df["id__DRUG"].isin(drugs)].reset_index(drop=True)
    lineage = df["cov__lineage_raw"].to_numpy()
    drug = df["id__DRUG"].to_numpy()

    # features: raw__ genotype + numeric coverage + id__DRUG + drug_class (cats).
    # EXCLUDE lineage one-hots (the grouping variable) — using them is circular.
    raw = [c for c in df.columns if c.startswith("raw__")]
    num = [c for c in ("cov__median_coverage", "cov__tb_breadth") if c in df.columns]
    cats = [c for c in ("id__DRUG", "cov__drug_class") if c in df.columns]
    X = df[raw + num + cats].copy()
    codes = {}
    for c in cats:                       # label-encode → LightGBM categorical
        X[c] = X[c].astype("category")
        codes[c] = list(X[c].cat.categories)
        X[c] = X[c].cat.codes
    y = df["y__binary"].to_numpy()
    cat_idx = [X.columns.get_loc(c) for c in cats]

    # lineage-LOO: OOF prediction on each held-out lineage
    oof = np.full(len(df), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() == 0 or len(np.unique(y[tr])) < 2:
            continue
        m = _lgbm()
        m.fit(X.iloc[tr], y[tr], categorical_feature=cat_idx)
        oof[te] = m.predict_proba(X.iloc[te])[:, 1]

    # per-drug: mean over held-out lineages + per-lineage
    per_drug = {}
    for d in sorted(set(drug)):
        pl = {}
        for L in LINEAGES:
            mask = (drug == d) & (lineage == L) & ~np.isnan(oof)
            if mask.sum() >= 20:
                pl[L] = _safe_auc(y[mask], oof[mask])
        aucs = [v for v in pl.values() if v is not None]
        per_drug[d] = {"multidrug_lineage_loo_auc": round(float(np.mean(aucs)), 4) if aucs else None,
                       "per_lineage": {k: (round(v, 3) if v else None) for k, v in pl.items()},
                       "n_R": int(y[drug == d].sum())}

    # compare to per-drug baseline (concordant lineage-LOO)
    base = {}
    if baseline_path and baseline_path.exists():
        b = json.loads(baseline_path.read_text())
        base = {d: b[d]["concordant"]["lineage_loo"] for d in b if b[d].get("concordant")}
    for d, r in per_drug.items():
        md = r["multidrug_lineage_loo_auc"]
        bl = base.get(d)
        r["baseline_perdrug_lineage_loo"] = bl
        r["transfer_delta"] = round(md - bl, 4) if (md is not None and bl is not None) else None

    manifest = {
        "experiment": "multidrug-transfer (T14 drug-class + id__DRUG single model)",
        "model": "LightGBM (fixed), lineage-LOO",
        "pool": sorted(set(drug)) if drugs else "all",
        "n_rows": len(df), "n_features": len(raw + num + cats),
        "categorical_features": cats,
        "per_drug": per_drug,
        "wall_time_seconds": round(time.time() - t0, 2),
    }
    out_name = "manifest.json" if not drugs else f"manifest_{'_'.join(sorted(drugs))}.json"
    (out_dir / out_name).write_text(json.dumps(manifest, indent=2))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description="Multi-drug transfer lineage-LOO evaluation (design/15 Part B, T14).")
    p.add_argument("--mart", type=Path, required=True)
    p.add_argument("--baseline", type=Path, default=Path("analysis/results/cv_eval_fulldata/summary.json"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/md_eval"))
    p.add_argument("--drugs", nargs="+", default=None,
                   help="restrict the transfer pool to these drugs (mechanism-grouped ablation)")
    args = p.parse_args()
    m = evaluate(args.mart, args.baseline, args.out, args.drugs)
    print(f"\n{'drug':<5}{'md_linLOO':>10}{'baseline':>10}{'Δtransfer':>10}  per-lineage")
    for d, r in sorted(m["per_drug"].items(), key=lambda kv: -(kv[1]['transfer_delta'] or -9)):
        md = r["multidrug_lineage_loo_auc"]; bl = r["baseline_perdrug_lineage_loo"]; dl = r["transfer_delta"]
        print(f"{d:<5}{(md if md else 0):>10.3f}{(bl if bl else 0):>10.3f}{(dl if dl is not None else 0):>+10.3f}"
              f"  {r['per_lineage']}")


if __name__ == "__main__":
    main()
