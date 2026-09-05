"""B5 — elastic-net logistic: unify feature selection + modelling.

Our pipeline is two stages: causal-concordance *selects* (EBM + CATE + FWL, heavy
`econml`), then L2-logistic *models*. Elastic-net does both in one fit — L1 zeros
coefficients (selection), L2 shares weight across correlated variants. We fit it on the
**full per-drug universe** with `(C, l1_ratio)` tuned by lineage-LOO, and ask:

  1. Does unified elastic-net **match** the two-stage concordance + L2 under lineage-LOO?
  2. Does it **select the same features** (Jaccard vs the concordant set)?

Match → pipeline simplification + reproducibility win (one deterministic fit, no
stochastic multi-voter selection). Underperform → causal selection earns its keep.

The elastic-net AUC is tuned on the same lineage-LOO it is reported on (mild optimism)
— that can only *favour* elastic-net, so if it still does not beat concordance+L2 the
conclusion is robust. In-memory, manifest only.

Run:
    python -m analysis.scripts.feature_mart.evaluate_elasticnet \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --deployed analysis/results/fulldata_ml/logistic_models/manifest.json \
        --out analysis/results/fulldata_ml/elasticnet
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

LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COV = ["cov__median_coverage", "cov__tb_breadth"]
DEPLOY = ["RIF", "INH", "EMB", "MXF", "LEV", "ETH", "KAN", "AMI"]
AGG = {"raw__n_mutations", "raw__mean_frs", "raw__mean_coverage", "raw__n_minor"}
C_GRID = [0.03, 0.1, 0.3]
L1_GRID = [0.2, 0.5, 0.8]


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _enet(C, l1):
    return LogisticRegression(penalty="elasticnet", solver="saga", C=C, l1_ratio=l1,
                              max_iter=1500, tol=1e-3)


def _loo(X, y, lineage, C, l1):
    oof = np.full(len(y), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        sc = StandardScaler().fit(X[tr])
        m = _enet(C, l1).fit(sc.transform(X[tr]), y[tr])
        oof[te] = m.predict_proba(sc.transform(X[te]))[:, 1]
    per = {L: _auc(y[(lineage == L) & ~np.isnan(oof)], oof[(lineage == L) & ~np.isnan(oof)])
           for L in LINEAGES if ((lineage == L) & ~np.isnan(oof)).sum() >= 20}
    aucs = [v for v in per.values() if v is not None]
    return (np.mean(aucs) if aucs else None,
            {k: (round(v, 3) if v else None) for k, v in per.items()})


def evaluate(marts_dir, causal_dir, deployed_path, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    deployed = {}
    if deployed_path and Path(deployed_path).exists():
        deployed = json.loads(Path(deployed_path).read_text()).get("per_drug", {})
    results = {}
    for drug in DEPLOY:
        marts = glob.glob(str(marts_dir / f"feature_mart_{drug}_*vFULL.parquet"))
        cpaths = glob.glob(str(causal_dir / drug / f"causal_concordance_{drug}_mutation_level.parquet"))
        if not marts or not cpaths:
            continue
        df = pd.read_parquet(marts[0])
        y = df["y__binary"].to_numpy()
        lineage = df["cov__lineage_raw"].to_numpy()
        universe = [c for c in df.columns if c.startswith("raw__") and c not in AGG
                    and not c.startswith("raw__gene_")]
        # prefilter the wide universe to the carrier-windowed top-K by |corr with y|
        # (saga is O(iters x nnz); the full ~1,800-feature fit is intractable). Elastic-net
        # still performs L1 selection WITHIN this pool; the concordant set is unioned in so
        # the two-stage pipeline's picks are always available for a fair comparison.
        carr = df[universe].sum()
        elig = carr[carr >= 10].index
        corr = df[elig].corrwith(df["y__binary"]).abs().sort_values(ascending=False)
        conc = set(pd.read_parquet(cpaths[0]).query("qf__causal_concordant")["feature_column"]) & set(universe)
        pool = sorted(set(corr.head(400).index) | conc)
        cov = [c for c in COV if c in df.columns]
        feats = pool + cov
        X = df[feats].to_numpy("float32")

        # tune (C, l1) by lineage-LOO
        best = (None, None, -1, None)
        for C in C_GRID:
            for l1 in L1_GRID:
                auc, per = _loo(X, y, lineage, C, l1)
                if auc is not None and auc > best[2]:
                    best = (C, l1, auc, per)
        C, l1, auc, per = best
        # final full-data fit at best params → selected features
        sc = StandardScaler().fit(X)
        final = _enet(C, l1).fit(sc.transform(X), y)
        nz = {feats[i] for i in np.nonzero(final.coef_[0])[0]}
        nz_mut = nz & set(universe)
        inter = nz_mut & conc
        jac = len(inter) / len(nz_mut | conc) if (nz_mut | conc) else None

        dep = deployed.get(drug, {})
        dep_auc = dep.get("honest_nested_loo_auc")
        results[drug] = {
            "enet_auc": round(float(auc), 4), "enet_C": C, "enet_l1_ratio": l1,
            "enet_per_lineage": per, "n_selected": len(nz_mut), "n_concordant": len(conc),
            "n_overlap": len(inter), "jaccard_vs_concordant": round(jac, 3) if jac is not None else None,
            "concordant_L2_auc": dep_auc,
            "enet_minus_concordantL2": round(auc - dep_auc, 4) if dep_auc else None,
        }
        print(f"{drug:<5} enet={auc:.4f} (C={C} l1={l1}) vs conc+L2 {dep_auc} "
              f"(Δ {results[drug]['enet_minus_concordantL2']:+}) | "
              f"selected {len(nz_mut)} vs conc {len(conc)} overlap {len(inter)} J={results[drug]['jaccard_vs_concordant']}")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "B5 elastic-net (unify selection+model) vs concordance+L2, lineage-LOO",
         "C_grid": C_GRID, "l1_grid": L1_GRID, "per_drug": results,
         "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return results


def main():
    p = argparse.ArgumentParser(description="B5 elastic-net logistic vs concordance+L2.")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--deployed", type=Path, default=Path("analysis/results/fulldata_ml/logistic_models/manifest.json"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/elasticnet"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.deployed, a.out)


if __name__ == "__main__":
    main()
