"""Mutation-class FE (E3/E4) — the synonymous-hitchhiker test + class enrichment.

Mutation class is deterministic from the GARC string in each `raw__<gene>_<mut>`
column, so this is leakage-free and needs no new infra (parses the existing marts).

E3  Split the FULL per-drug feature universe by class and evaluate three feature
    sets under BOTH random-CV and lineage-LOO:
      (all) everything   (functional) drop synonymous   (synonymous) synonymous only
    The thesis-level result: synonymous-only should score HIGH under random-CV and
    COLLAPSE to ~chance under lineage-LOO — a direct demonstration that silent
    mutations are lineage hitchhikers, not resistance signal. Functional-only should
    match `all` under lineage-LOO (synonymous safely droppable).

E4  Class enrichment among the causally-concordant + lineage-invariant drivers vs the
    background universe (Fisher exact) — do the survivors avoid synonymous and favour
    functional classes?

Fixed LightGBM (trees exploit hitchhikers hardest → starkest collapse). In-memory,
manifests only (storage-constrained machine).

Run (causal venv):
    python -m analysis.scripts.feature_mart.evaluate_mutation_class \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --invariance analysis/results/fulldata_ml/invariance_eval/manifest.json \
        --out analysis/results/fulldata_ml/mutation_class
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

from analysis.scripts.feature_mart.mutation_class import mut_class  # noqa: F401
from scipy.stats import fisher_exact
from sklearn.metrics import roc_auc_score
from sklearn.model_selection import StratifiedKFold

RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COV = ["cov__median_coverage", "cov__tb_breadth"]
_LGBM = dict(n_estimators=300, learning_rate=0.05, num_leaves=31, subsample=0.8,
             colsample_bytree=0.6, random_state=RANDOM_STATE, n_jobs=-1, verbose=-1)
_AGG = {"raw__n_mutations", "raw__mean_frs", "raw__mean_coverage", "raw__n_minor"}
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


def evaluate(marts_dir, causal_dir, invariance_path, out_dir):
    t0 = time.time()
    inv = {}
    if invariance_path and Path(invariance_path).exists():
        inv = json.loads(Path(invariance_path).read_text()).get("per_drug", {})
    e3, e4 = {}, {}
    for cpath in sorted(glob.glob(str(causal_dir / "*" / "causal_concordance_*_mutation_level.parquet"))):
        drug = Path(cpath).name.split("_")[2]
        marts = glob.glob(str(marts_dir / f"feature_mart_{drug}_*vFULL.parquet"))
        if not marts:
            continue
        df = pd.read_parquet(marts[0])
        y = df["y__binary"].to_numpy()
        lineage = df["cov__lineage_raw"].to_numpy()
        universe = [c for c in df.columns if c.startswith("raw__")
                    and c not in _AGG and not c.startswith("raw__gene_")]
        cls = {c: mut_class(c) for c in universe}
        cov = [c for c in COV if c in df.columns]
        Xcov = df[cov].to_numpy(dtype="float32")

        def feats(keep):
            cols = [c for c in universe if keep(cls[c])]
            X = np.hstack([df[cols].to_numpy(dtype="float32"), Xcov]) if cols else Xcov
            return X, len(cols)

        sets = {
            "all": (lambda k: True),
            "functional_no_synonymous": (lambda k: k != "synonymous"),
            "synonymous_only": (lambda k: k == "synonymous"),
        }
        res = {"class_counts_universe": {k: sum(v == k for v in cls.values())
                                        for k in set(cls.values())}}
        for name, keep in sets.items():
            X, n = feats(keep)
            res[name] = {"n_features": n,
                         "random_cv_auc": _random_cv(X, y) if n else None,
                         "lineage_loo_auc": _lineage_loo(X, y, lineage) if n else None}
        # collapse magnitude for synonymous-only
        so = res["synonymous_only"]
        res["synonymous_leakage_gap"] = (round(so["random_cv_auc"] - so["lineage_loo_auc"], 4)
                                         if so["random_cv_auc"] and so["lineage_loo_auc"] else None)
        e3[drug] = res
        print(f"{drug:<5} | all loo={res['all']['lineage_loo_auc']} "
              f"func={res['functional_no_synonymous']['lineage_loo_auc']} "
              f"(n_syn={res['synonymous_only']['n_features']}) | "
              f"SYN random={so['random_cv_auc']} loo={so['lineage_loo_auc']} gap={res['synonymous_leakage_gap']}")

        # ---- E4: class enrichment among concordant + invariant drivers ----
        conc = set(pd.read_parquet(cpath).query("qf__causal_concordant")["feature_column"]) & set(universe)
        invariant = set(inv.get(drug, {}).get("invariant_features", [])) & set(universe)
        bg = pd.Series(cls)
        def enrich(sel):
            out = {}
            for k in ("nonsynonymous", "promoter_nt", "indel", "synonymous", "nonsense_Z"):
                a = sum(cls[c] == k for c in sel)                 # selected & class k
                b = len(sel) - a                                  # selected & not k
                cc = int((bg == k).sum()) - a                     # bg & class k
                d = len(universe) - len(sel) - cc                 # bg & not k
                if a + b > 0 and cc + d > 0:
                    orr, p = fisher_exact([[a, b], [cc, d]], alternative="two-sided")
                    out[k] = {"n_selected": a, "odds_ratio": round(float(orr), 2) if np.isfinite(orr) else None,
                              "p": round(float(p), 4)}
            return out
        e4[drug] = {"n_concordant": len(conc), "n_invariant": len(invariant),
                    "concordant_enrichment": enrich(conc),
                    "invariant_enrichment": enrich(invariant) if invariant else None}
    wall = round(time.time() - t0, 2)
    for name, res in (("e3_class_feature_sets", e3), ("e4_class_enrichment", e4)):
        d = out_dir / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "manifest.json").write_text(json.dumps(
            {"experiment": name, "judged_by": "random-CV vs lineage-LOO AUC (E3); Fisher enrichment (E4)",
             "per_drug": res, "wall_time_seconds": wall}, indent=2))
    return {"e3": e3, "e4": e4}


def main():
    p = argparse.ArgumentParser(description="Mutation-class FE (E3 synonymous-hitchhiker, E4 enrichment).")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--invariance", type=Path,
                   default=Path("analysis/results/fulldata_ml/invariance_eval/manifest.json"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/mutation_class"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.invariance, a.out)


if __name__ == "__main__":
    main()
