"""B4 Firth/rare-event logistic + C7 site-LOO rigor.

Two cheap, high-value techniques on top of the deployed L2-logistic:

B4 — **Firth penalized logistic** corrects the small-sample / separation bias that
     standard logistic suffers when the positive class is rare (BDQ 1.4% R, CFZ, LZD).
     Firth adds the Jeffreys-prior penalty 0.5·log|I(β)|, implemented here as IRLS with
     the hat-diagonal score adjustment (no new dependency). Compared under lineage-LOO
     against the deployed L2-logistic and a class-weight-balanced L2, on concordant
     features. Question: does principled rare-event correction lift the data-poor tier?

C7 — **site-LOO** is the third, stricter CV regime. CRyPTIC's collection site is
     embedded in every UNIQUEID (`site.NN.…`); country metadata is 53% missing, so site
     is the fully-populated geographic/protocol-effect grouping. Leave-one-site-out
     (10 usable sites), L2-logistic, vs the lineage-LOO baseline — does the honest AUC
     survive an even harder split?

In-memory, manifest only.

Run:
    python -m analysis.scripts.feature_mart.evaluate_firth_siteloo \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --out analysis/results/fulldata_ml/firth_siteloo
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
DATA_POOR = {"BDQ", "CFZ", "LZD", "DLM"}
MIN_SITE = 100
MIN_CLASS = 15


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _firth_fit(X, y, max_iter=50, tol=1e-6):
    """Firth-penalized logistic via IRLS with hat-diagonal score adjustment.
    X already standardised; a bias column is appended here."""
    Xb = np.hstack([np.ones((len(X), 1)), X])
    n, p = Xb.shape
    beta = np.zeros(p)
    for _ in range(max_iter):
        eta = np.clip(Xb @ beta, -30, 30)
        pi = 1.0 / (1.0 + np.exp(-eta))
        W = np.clip(pi * (1 - pi), 1e-9, None)
        XtWX = Xb.T @ (W[:, None] * Xb) + 1e-6 * np.eye(p)
        inv = np.linalg.inv(XtWX)
        h = W * np.einsum("ij,jk,ik->i", Xb, inv, Xb)          # hat diagonal
        U = Xb.T @ (y - pi + h * (0.5 - pi))                    # Firth-adjusted score
        step = inv @ U
        beta = beta + step
        if np.max(np.abs(step)) < tol:
            break
    return beta


def _firth_proba(beta, X):
    eta = np.clip(np.hstack([np.ones((len(X), 1)), X]) @ beta, -30, 30)
    return 1.0 / (1.0 + np.exp(-eta))


MODELS = {
    "l2_logistic": lambda: ("sk", LogisticRegression(C=0.1, max_iter=2000)),
    "l2_balanced": lambda: ("sk", LogisticRegression(C=0.1, max_iter=2000, class_weight="balanced")),
    "firth": lambda: ("firth", None),
}


def _fit_predict(kind, est, Xtr, ytr, Xte):
    sc = StandardScaler().fit(Xtr)
    Ztr, Zte = sc.transform(Xtr), sc.transform(Xte)
    if kind == "firth":
        beta = _firth_fit(Ztr, ytr)
        return _firth_proba(beta, Zte)
    return est.fit(Ztr, ytr).predict_proba(Zte)[:, 1]


def _grouped_loo(X, y, groups, group_vals, model_key):
    oof = np.full(len(y), np.nan)
    for g in group_vals:
        te = groups == g
        tr = ~te
        if te.sum() < 20 or len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
            continue
        kind, est = MODELS[model_key]()
        oof[te] = _fit_predict(kind, est, X[tr], y[tr], X[te])
    per = {}
    for g in group_vals:
        m = (groups == g) & ~np.isnan(oof)
        if m.sum() >= 20 and len(np.unique(y[m])) == 2:
            per[str(g)] = _auc(y[m], oof[m])
    aucs = [v for v in per.values() if v is not None]
    return {"auc": round(float(np.mean(aucs)), 4) if aucs else None,
            "stability_std": round(float(np.std(aucs)), 4) if len(aucs) > 1 else None,
            "n_groups": len(aucs), "per_group": {k: round(v, 3) for k, v in per.items()}}


def evaluate(marts_dir, causal_dir, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    firth_res, site_res = {}, {}
    for cpath in sorted(glob.glob(str(causal_dir / "*" / "causal_concordance_*_mutation_level.parquet"))):
        drug = Path(cpath).name.split("_")[2]
        marts = glob.glob(str(marts_dir / f"feature_mart_{drug}_*vFULL.parquet"))
        if not marts:
            continue
        df = pd.read_parquet(marts[0])
        y = df["y__binary"].to_numpy()
        lineage = df["cov__lineage_raw"].to_numpy()
        site = df["id__UNIQUEID"].str.extract(r"^site\.(\d+)\.")[0].fillna("?").to_numpy()
        conc = [c for c in pd.read_parquet(cpath).query("qf__causal_concordant")["feature_column"] if c in df.columns]
        if not conc:
            continue
        cov = [c for c in COV if c in df.columns]
        X = df[conc + cov].to_numpy("float32")

        # B4 — Firth vs L2 vs balanced, lineage-LOO
        arms = {k: _grouped_loo(X, y, lineage, LINEAGES, k) for k in MODELS}
        firth_res[drug] = {
            "n_concordant": len(conc), "is_data_poor": drug in DATA_POOR,
            **{k: {"auc": v["auc"], "stab": v["stability_std"]} for k, v in arms.items()},
            "firth_minus_l2": (round(arms["firth"]["auc"] - arms["l2_logistic"]["auc"], 4)
                               if arms["firth"]["auc"] and arms["l2_logistic"]["auc"] else None),
            "firth_per_lineage": arms["firth"]["per_group"],
        }
        # C7 — site-LOO (L2-logistic), vs lineage-LOO
        sc = pd.Series(site).value_counts()
        usable = [s for s in sc.index if s != "?" and sc[s] >= MIN_SITE
                  and (y[site == s].sum() >= MIN_CLASS) and ((site == s).sum() - y[site == s].sum() >= MIN_CLASS)]
        siteloo = _grouped_loo(X, y, site, usable, "l2_logistic")
        lin = arms["l2_logistic"]
        site_res[drug] = {"site_loo_auc": siteloo["auc"], "site_loo_stab": siteloo["stability_std"],
                          "n_sites": siteloo["n_groups"], "lineage_loo_auc": lin["auc"],
                          "site_minus_lineage": (round(siteloo["auc"] - lin["auc"], 4)
                                                 if siteloo["auc"] and lin["auc"] else None),
                          "weakest_sites": dict(sorted(siteloo["per_group"].items(), key=lambda kv: kv[1])[:3])}
        print(f"{drug:<5} | Firth {arms['firth']['auc']} vs L2 {arms['l2_logistic']['auc']} "
              f"bal {arms['l2_balanced']['auc']} (Δfirth {firth_res[drug]['firth_minus_l2']:+}) "
              f"| site-LOO {siteloo['auc']} vs lin {lin['auc']} (Δ {site_res[drug]['site_minus_lineage']:+})")
    for name, res in (("firth_eval", firth_res), ("site_loo_eval", site_res)):
        (out_dir / f"{name}.json").write_text(json.dumps(
            {"experiment": name, "per_drug": res, "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return firth_res, site_res


def main():
    p = argparse.ArgumentParser(description="B4 Firth logistic + C7 site-LOO.")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/firth_siteloo"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.out)


if __name__ == "__main__":
    main()
