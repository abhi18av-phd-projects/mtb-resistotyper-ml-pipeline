"""Per-mutation quality (E1/E2) + model/technique comparison — storage-light.

Everything is built IN MEMORY from the existing full DB + vFULL marts; only small
JSON manifests are written (the machine is storage-constrained — no new marts).

E1  FRS-weighted presence   binary vs FRS-valued vs binary+FRS   (per-mutation quality as dose)
E2  quality-gating sweep    presence gated at FRS >= tau, + non-minor variant
MODELS  fixed concordant binary features, lineage-LOO across a learner panel
        (LightGBM, HistGBM, RandomForest, ExtraTrees, L2-Logistic) — model comparison
        for the tiers where we underperform.

All judged by lineage-LOO AUC + per-lineage equity (the honest yardstick).

Run (causal venv):
    python -m analysis.scripts.feature_mart.evaluate_quality_and_models \
        --db /tmp/cryptic-full.duckdb --marts analysis/results/fulldata_ml/marts \
        --causal analysis/results/fulldata_ml/causal --out analysis/results/fulldata_ml
"""

from __future__ import annotations

import argparse
import glob
import json
import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.ensemble import ExtraTreesClassifier, HistGradientBoostingClassifier, RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COV = ["cov__median_coverage", "cov__tb_breadth"]
_LGBM = dict(n_estimators=300, learning_rate=0.05, num_leaves=31, subsample=0.8,
             colsample_bytree=0.6, random_state=RANDOM_STATE, n_jobs=-1, verbose=-1)


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _parse(col):
    gene, _, mut = col.removeprefix("raw__").partition("_")
    return gene, mut


def _quality_matrices(con, uids, conc):
    """In-memory per-determinant FRS / IS_MINOR aligned to `uids` order.
    Returns (frs[n,k], minor[n,k]) with 0 where the isolate lacks the call."""
    gm = pd.DataFrame([_parse(c) for c in conc], columns=["GENE", "MUTATION"])
    gm["feat"] = conc
    con.register("qgm", gm)
    con.register("quids", pd.DataFrame({"UNIQUEID": uids}))
    q = con.execute("""
        SELECT m.UNIQUEID, g.feat,
               MAX(CASE WHEN m.FRS = 'inf' THEN 1.0 ELSE m.FRS END) AS frs,
               MAX(m.IS_MINOR::INT) AS minor
        FROM mutations m JOIN qgm g ON g.GENE=m.GENE AND g.MUTATION=m.MUTATION
        WHERE m.UNIQUEID IN (SELECT UNIQUEID FROM quids)
        GROUP BY m.UNIQUEID, g.feat
    """).fetchdf()
    idx = {u: i for i, u in enumerate(uids)}
    fpos = {f: j for j, f in enumerate(conc)}
    frs = np.zeros((len(uids), len(conc)), dtype="float32")
    minor = np.zeros((len(uids), len(conc)), dtype="float32")
    for u, f, fr, mi in q.itertuples(index=False):
        i, j = idx.get(u), fpos.get(f)
        if i is not None and j is not None:
            frs[i, j] = fr if fr is not None else 0.0
            minor[i, j] = mi or 0.0
    con.unregister("qgm"); con.unregister("quids")
    return frs, minor


def _loo(X, y, lineage, model_fn):
    oof = np.full(len(y), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        oof[te] = model_fn().fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    per = {}
    for L in LINEAGES:
        m = (lineage == L) & ~np.isnan(oof)
        if m.sum() >= 20:
            per[L] = _auc(y[m], oof[m])
    aucs = [v for v in per.values() if v is not None]
    return (round(float(np.mean(aucs)), 4) if aucs else None,
            round(float(np.std(aucs)), 4) if len(aucs) > 1 else None,
            {k: (round(v, 3) if v else None) for k, v in per.items()})


MODELS = {
    "lightgbm": lambda: LGBMClassifier(**_LGBM),
    "histgbm": lambda: HistGradientBoostingClassifier(max_iter=300, learning_rate=0.05,
                                                      max_leaf_nodes=31, random_state=RANDOM_STATE),
    "randomforest": lambda: RandomForestClassifier(n_estimators=400, n_jobs=-1,
                                                   random_state=RANDOM_STATE),
    "extratrees": lambda: ExtraTreesClassifier(n_estimators=400, n_jobs=-1,
                                              random_state=RANDOM_STATE),
    "logistic_l2": lambda: make_pipeline(StandardScaler(),
                                        LogisticRegression(max_iter=1000, C=1.0)),
}


def evaluate(db, marts_dir, causal_dir, out_dir):
    t0 = time.time()
    con = duckdb.connect(str(db), read_only=True)
    e1, e2, mc = {}, {}, {}
    for cpath in sorted(glob.glob(str(causal_dir / "*" / "causal_concordance_*_mutation_level.parquet"))):
        drug = Path(cpath).name.split("_")[2]
        marts = glob.glob(str(marts_dir / f"feature_mart_{drug}_*vFULL.parquet"))
        if not marts:
            continue
        df = pd.read_parquet(marts[0])
        y = df["y__binary"].to_numpy()
        lineage = df["cov__lineage_raw"].to_numpy()
        uids = df["id__UNIQUEID"].tolist()
        conc = [c for c in pd.read_parquet(cpath).query("qf__causal_concordant")["feature_column"] if c in df.columns]
        cov = [c for c in COV if c in df.columns]
        if not conc:
            continue
        Xbin = df[conc].to_numpy(dtype="float32")
        Xcov = df[cov].to_numpy(dtype="float32")
        frs, minor = _quality_matrices(con, uids, conc)

        def stack(*mats):
            return np.hstack([m for m in mats])

        lg = lambda: LGBMClassifier(**_LGBM)  # noqa: E731
        # ---- E1: encoding of the concordant features (quality as dose) ----
        a_bin = _loo(stack(Xbin, Xcov), y, lineage, lg)
        a_frs = _loo(stack(frs, Xcov), y, lineage, lg)
        a_both = _loo(stack(Xbin, frs, Xcov), y, lineage, lg)
        e1[drug] = {"n_concordant": len(conc),
                    "binary": {"auc": a_bin[0], "stab": a_bin[1], "per_lineage": a_bin[2]},
                    "frs_valued": {"auc": a_frs[0], "stab": a_frs[1], "per_lineage": a_frs[2]},
                    "binary_plus_frs": {"auc": a_both[0], "stab": a_both[1], "per_lineage": a_both[2]},
                    "delta_frs_minus_binary": (round(a_frs[0] - a_bin[0], 4) if a_bin[0] and a_frs[0] else None)}
        # ---- E2: quality gating sweep ----
        gates = {}
        for tau in (0.0, 0.25, 0.5, 0.75, 0.9):
            Xg = (frs >= tau).astype("float32") if tau > 0 else Xbin
            r = _loo(stack(Xg, Xcov), y, lineage, lg)
            gates[f"tau_{tau}"] = {"auc": r[0], "per_lineage": r[2]}
        Xnm = (Xbin * (1 - minor)).astype("float32")  # drop minor calls
        r_nm = _loo(stack(Xnm, Xcov), y, lineage, lg)
        gates["non_minor_only"] = {"auc": r_nm[0], "per_lineage": r_nm[2]}
        best = max((k for k in gates if gates[k]["auc"] is not None),
                   key=lambda k: gates[k]["auc"], default=None)
        e2[drug] = {"n_concordant": len(conc), "gates": gates,
                    "best_gate": best, "best_auc": gates[best]["auc"] if best else None,
                    "baseline_tau0": gates["tau_0.0"]["auc"]}
        # ---- MODELS: learner panel on binary concordant features ----
        mc[drug] = {"n_concordant": len(conc)}
        Xm = stack(Xbin, Xcov)
        for name, fn in MODELS.items():
            r = _loo(Xm, y, lineage, fn)
            mc[drug][name] = {"auc": r[0], "stab": r[1], "per_lineage": r[2]}
        print(f"{drug:<5} conc={len(conc):>3} | E1 bin={a_bin[0]} frs={a_frs[0]} both={a_both[0]} "
              f"| E2 best={best}({e2[drug]['best_auc']}) | "
              f"models lgbm={mc[drug]['lightgbm']['auc']} hgb={mc[drug]['histgbm']['auc']} "
              f"rf={mc[drug]['randomforest']['auc']} logit={mc[drug]['logistic_l2']['auc']}")
    con.close()
    wall = round(time.time() - t0, 2)
    for name, res in (("frs_weighting_eval", e1), ("quality_gating_eval", e2), ("model_comparison", mc)):
        d = out_dir / name
        d.mkdir(parents=True, exist_ok=True)
        (d / "manifest.json").write_text(json.dumps(
            {"experiment": name, "judged_by": "lineage-LOO AUC + per-lineage equity",
             "per_drug": res, "wall_time_seconds": wall}, indent=2))
    return {"e1": e1, "e2": e2, "models": mc, "wall": wall}


def main():
    p = argparse.ArgumentParser(description="Per-mutation quality (E1/E2) + model comparison, lineage-LOO.")
    p.add_argument("--db", type=Path, default=Path("/tmp/cryptic-full.duckdb"))
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml"))
    a = p.parse_args()
    evaluate(a.db, a.marts, a.causal, a.out)


if __name__ == "__main__":
    main()
