"""Genome-annotation FE, tier 1 — PE/PPE/ESX/IS region flags (design/20).

The region analogue of the E3 synonymous-hitchhiker test. PE/PPE/ESX/IS genes are
hypervariable, error-prone-mapping regions — a canonical population-structure /
hitchhiker source, and 19% of the mart universe. Parsed from the `GENE` name (zero
download, leakage-free). Per drug, split the full universe by region and evaluate under
BOTH random-CV and lineage-LOO (LightGBM, as in E3 — trees exploit hitchhikers hardest):

  (all)  everything   (core) exclude PE/PPE/ESX/IS   (pe_ppe_only) only those regions

Expectation (mirroring synonymous): the region-only model scores high under random-CV
and collapses to ~chance under lineage-LOO (confirming region hitchhikers), and
excluding the regions is neutral-to-positive for the core model (a cleaner universe).

Run:
    python -m analysis.scripts.feature_mart.evaluate_region_annotation \
        --marts analysis/results/fulldata_ml/marts --out analysis/results/fulldata_ml/region_annotation
"""

from __future__ import annotations

import argparse
import glob
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COV = ["cov__median_coverage", "cov__tb_breadth"]
AGG = {"raw__n_mutations", "raw__mean_frs", "raw__mean_coverage", "raw__n_minor"}
_LGBM = dict(n_estimators=300, learning_rate=0.05, num_leaves=31, subsample=0.8,
             colsample_bytree=0.6, random_state=RANDOM_STATE, n_jobs=-1, verbose=-1)
_PE = re.compile(r"^PE(\d|_)")
_PPE = re.compile(r"^PPE(\d|_)")


def region(col: str) -> str:
    g = col.removeprefix("raw__").split("_")[0]
    if _PPE.match(col.removeprefix("raw__")):
        return "PPE"
    if _PE.match(col.removeprefix("raw__")):
        return "PE"
    if g.lower().startswith("esx"):
        return "ESX"
    if g.startswith("IS") and any(ch.isdigit() for ch in g):
        return "IS"
    return "core"


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _lineage_loo(X, y, lineage):
    oof = np.full(len(y), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        oof[te] = LGBMClassifier(**_LGBM).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    aucs = [_auc(y[(lineage == L) & ~np.isnan(oof)], oof[(lineage == L) & ~np.isnan(oof)])
            for L in LINEAGES if ((lineage == L) & ~np.isnan(oof)).sum() >= 20]
    aucs = [a for a in aucs if a is not None]
    return round(float(np.mean(aucs)), 4) if aucs else None


def _random_cv(X, y):
    skf = StratifiedKFold(n_splits=5, shuffle=True, random_state=RANDOM_STATE)
    oof = np.full(len(y), np.nan)
    for tr, te in skf.split(X, y):
        oof[te] = LGBMClassifier(**_LGBM).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    return round(_auc(y, oof) or float("nan"), 4)


def evaluate(marts_dir, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    results = {}
    for m in sorted(glob.glob(str(marts_dir / "feature_mart_*_*vFULL.parquet"))):
        drug = Path(m).name.split("_")[2]
        df = pd.read_parquet(m)
        y = df["y__binary"].to_numpy()
        lineage = df["cov__lineage_raw"].to_numpy()
        universe = [c for c in df.columns if c.startswith("raw__") and c not in AGG
                    and not c.startswith("raw__gene_")]
        reg = {c: region(c) for c in universe}
        cov = [c for c in COV if c in df.columns]
        Xcov = df[cov].to_numpy("float32")
        pe_ppe = [c for c in universe if reg[c] in ("PE", "PPE", "ESX", "IS")]
        core = [c for c in universe if reg[c] == "core"]

        def loo(cols):
            X = np.hstack([df[cols].to_numpy("float32"), Xcov]) if cols else Xcov
            return _lineage_loo(X, y, lineage)

        res = {"n_universe": len(universe), "n_pe_ppe": len(pe_ppe), "n_core": len(core),
               "all_lineage_loo": loo(universe), "core_lineage_loo": loo(core)}
        if pe_ppe:
            Xpp = np.hstack([df[pe_ppe].to_numpy("float32"), Xcov])
            res["pe_ppe_only_random_cv"] = _random_cv(Xpp, y)
            res["pe_ppe_only_lineage_loo"] = _lineage_loo(Xpp, y, lineage)
            res["pe_ppe_leakage_gap"] = (round(res["pe_ppe_only_random_cv"] - res["pe_ppe_only_lineage_loo"], 4)
                                         if res["pe_ppe_only_random_cv"] and res["pe_ppe_only_lineage_loo"] else None)
        res["core_minus_all"] = (round(res["core_lineage_loo"] - res["all_lineage_loo"], 4)
                                 if res["core_lineage_loo"] and res["all_lineage_loo"] else None)
        results[drug] = res
        print(f"{drug:<5} pe/ppe={len(pe_ppe):>3}/{len(universe):<4} | all={res['all_lineage_loo']} "
              f"core={res['core_lineage_loo']} (Δ {res['core_minus_all']:+}) | "
              f"PEPPE-only rand={res.get('pe_ppe_only_random_cv')} loo={res.get('pe_ppe_only_lineage_loo')} "
              f"gap={res.get('pe_ppe_leakage_gap')}")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "region annotation (PE/PPE/ESX/IS) — hitchhiker collapse + core-vs-all",
         "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return results


def main():
    p = argparse.ArgumentParser(description="Region annotation (PE/PPE) FE experiment (design/20).")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/region_annotation"))
    a = p.parse_args()
    evaluate(a.marts, a.out)


if __name__ == "__main__":
    main()
