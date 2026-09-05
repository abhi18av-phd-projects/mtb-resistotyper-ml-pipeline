"""MIC regression (T8/T9) — does the graded LOG2MIC phenotype recover more
generalisable signal than binary R/S? (design/15 Part D / design/18 Track 4)

The binary phenotype thresholds MIC at the ECOFF/breakpoint, discarding the
gradient and making near-breakpoint isolates noisy coin-flips. CRyPTIC measured
the actual plate MIC, so we can regress on LOG2MIC directly. Two questions,
both under honest lineage-LOO on the same concordant features + rows:

1. **Discrimination** — does an LGBM *regressor* on LOG2MIC, used as an R-score,
   match or beat the LGBM *classifier* on y__binary (AUC vs the binary label)?
   If yes, the gradient carries extra, generalisable resistance information.
2. **Graded accuracy** — Spearman(pred, true LOG2MIC) and RMSE on the held-out
   lineage: a clinically richer output (how resistant, not just yes/no) that the
   binary model cannot produce.

Same fixed-hyperparameter LGBM for classifier and regressor to isolate the
*target* (binary vs graded), not architecture.

Run:
    python -m analysis.scripts.feature_mart.evaluate_mic_regression \
        --db /tmp/cryptic-full.duckdb --marts analysis/results/fulldata_ml/marts \
        --causal analysis/results/fulldata_ml/causal --out analysis/results/fulldata_ml/mic_regression
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
from lightgbm import LGBMClassifier, LGBMRegressor
from scipy.stats import spearmanr
from sklearn.metrics import roc_auc_score

RANDOM_STATE = 42
LINEAGES = ["lineage1", "lineage2", "lineage3", "lineage4"]
COV = ["cov__median_coverage", "cov__tb_breadth"]
_PARAMS = dict(n_estimators=300, learning_rate=0.05, num_leaves=31,
               subsample=0.8, colsample_bytree=0.6, random_state=RANDOM_STATE,
               n_jobs=-1, verbose=-1)


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _load_log2mic(db: Path, drug: str) -> pd.DataFrame:
    con = duckdb.connect(str(db), read_only=True)
    df = con.execute(
        "SELECT UNIQUEID AS id__UNIQUEID, LOG2MIC FROM ukmyc_phenotypes "
        "WHERE DRUG = ? AND LOG2MIC IS NOT NULL", [drug]).fetchdf()
    con.close()
    # one isolate can be plated more than once; keep the median LOG2MIC
    return df.groupby("id__UNIQUEID", as_index=False)["LOG2MIC"].median()


def _lineage_loo(df, feats, drug_db, out_dir):  # noqa: ARG001
    y = df["y__binary"].to_numpy()
    z = df["LOG2MIC"].to_numpy()
    lineage = df["cov__lineage_raw"].to_numpy()
    X = df[feats].to_numpy()
    clf_oof = np.full(len(df), np.nan)
    reg_oof = np.full(len(df), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        clf_oof[te] = LGBMClassifier(**_PARAMS).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
        reg_oof[te] = LGBMRegressor(**_PARAMS).fit(X[tr], z[tr]).predict(X[te])
    per = {}
    for L in LINEAGES:
        m = (lineage == L) & ~np.isnan(clf_oof)
        if m.sum() < 20:
            continue
        sp = spearmanr(reg_oof[m], z[m]).correlation
        per[L] = {
            "clf_auc": round(_auc(y[m], clf_oof[m]) or float("nan"), 3),
            "reg_auc": round(_auc(y[m], reg_oof[m]) or float("nan"), 3),
            "reg_spearman": round(float(sp), 3) if sp == sp else None,
            "reg_rmse": round(float(np.sqrt(np.mean((reg_oof[m] - z[m]) ** 2))), 3),
        }

    def _mean(k):
        vals = [v[k] for v in per.values() if v[k] is not None and v[k] == v[k]]
        return round(float(np.mean(vals)), 4) if vals else None
    return {"clf_auc": _mean("clf_auc"), "reg_auc": _mean("reg_auc"),
            "reg_spearman": _mean("reg_spearman"), "reg_rmse": _mean("reg_rmse"),
            "per_lineage": per}


def evaluate(db: Path, marts_dir: Path, causal_dir: Path, out_dir: Path) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    results = {}
    for cpath in sorted(glob.glob(str(causal_dir / "*" / "causal_concordance_*_mutation_level.parquet"))):
        drug = Path(cpath).name.split("_")[2]
        marts = glob.glob(str(marts_dir / f"feature_mart_{drug}_*vFULL.parquet"))
        if not marts:
            continue
        df = pd.read_parquet(marts[0])
        mic = _load_log2mic(db, drug)
        df = df.merge(mic, on="id__UNIQUEID", how="inner")
        if len(df) < 200:
            continue
        conc = [c for c in pd.read_parquet(cpath).query("qf__causal_concordant")["feature_column"] if c in df.columns]
        feats = conc + [c for c in COV if c in df.columns]
        r = _lineage_loo(df, feats, drug, out_dir)
        r["n"] = len(df)
        r["n_concordant"] = len(conc)
        r["auc_delta_reg_minus_clf"] = (round(r["reg_auc"] - r["clf_auc"], 4)
                                        if r["reg_auc"] is not None and r["clf_auc"] is not None else None)
        results[drug] = r
        print(f"{drug:<5} n={len(df):>6} conc={len(conc):>3} | clf_auc={r['clf_auc']} "
              f"reg_auc={r['reg_auc']} Δ={r['auc_delta_reg_minus_clf']:+} | "
              f"spearman={r['reg_spearman']} rmse={r['reg_rmse']}")
    manifest = {"experiment": "MIC regression (T8/T9) — LOG2MIC regressor vs binary classifier, lineage-LOO",
                "model": "LightGBM (fixed hyperparams), concordant features + coverage",
                "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description="MIC regression (T8/T9) lineage-LOO evaluation.")
    p.add_argument("--db", type=Path, default=Path("/tmp/cryptic-full.duckdb"))
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/mic_regression"))
    args = p.parse_args()
    evaluate(args.db, args.marts, args.causal, args.out)


if __name__ == "__main__":
    main()
