"""Causal columns AS features (C1/C2/C4) — do causal-adjusted quantities help as
model inputs, or only as selectors? (design/15 Part D / design/18 Track 4)

So far the causal machinery is used purely as a *selector* (concordance votes pick
which mutations enter the model as raw binary features). C1/C2/C4 ask the opposite:
does feeding a causal-*adjusted* quantity in as a feature improve generalisation?

We test the most defensible, leakage-free variant — a **fold-internal causal risk
score** (C2). For each held-out lineage, on the TRAIN lineages only we estimate each
concordant mutation's log-odds weight w_m = logit(P(R | carrier)) − logit(P(R | base))
(Haldane-corrected), then score each isolate as Σ_m w_m · x_{i,m}. The weights never
see the held-out lineage, so this is honest. A lineage-residualised variant (C1) uses
weights estimated after regressing y on the train lineage dummies (removing the
lineage-confounded component of each effect).

Three feature sets, same fixed LGBM, lineage-LOO:
  (A) raw concordant binary + coverage         — the existing baseline
  (B) causal risk score ALONE + coverage       — a compact deconfounded feature
  (C) raw concordant + risk score + coverage   — augmentation

The bet, and the likely-negative worth recording: a nonlinear tree already extracts
per-mutation effects from the raw binaries, so (C) probably ≈ (A); (B) tests whether
a compact causal score regularises better on data-poor drugs. Either way it settles
whether causal quantities earn their place as features or stay selectors.

Run:
    python -m analysis.scripts.feature_mart.evaluate_causal_features \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --out analysis/results/fulldata_ml/causal_features
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
_PARAMS = dict(n_estimators=300, learning_rate=0.05, num_leaves=31,
               subsample=0.8, colsample_bytree=0.6, random_state=RANDOM_STATE,
               n_jobs=-1, verbose=-1)


def _auc(y, p):
    return float(roc_auc_score(y, p)) if len(np.unique(y)) == 2 else None


def _logit(p):  # Haldane-corrected
    p = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(p / (1 - p))


def _risk_weights(Xtr, ytr, resid_lineage=None):
    """Per-mutation log-odds weight from TRAIN only. If resid_lineage given,
    subtract the per-lineage base rate (C1 lineage-residual) instead of the
    grand base rate (C2)."""
    base = _logit(ytr.mean())
    w = np.zeros(Xtr.shape[1])
    for j in range(Xtr.shape[1]):
        carr = Xtr[:, j] > 0
        if carr.sum() >= 10 and (~carr).sum() >= 10:
            if resid_lineage is None:
                w[j] = _logit(ytr[carr].mean()) - base
            else:  # C1: average within-lineage carrier log-odds lift
                lifts = []
                for L in np.unique(resid_lineage):
                    m = resid_lineage == L
                    c = carr & m
                    if c.sum() >= 5 and m.sum() - c.sum() >= 5 and len(np.unique(ytr[m])) == 2:
                        lifts.append(_logit(ytr[c].mean()) - _logit(ytr[m].mean()))
                w[j] = float(np.mean(lifts)) if lifts else 0.0
    return w


def _loo_auc(df, build_X, y, lineage, resid=False):
    """build_X(Xtr_raw, tr_idx, te_idx, w_c2, w_c1) -> (Xtr, Xte); lineage-LOO AUC."""
    oof = np.full(len(df), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        Xtr, Xte = build_X(tr, te)
        oof[te] = LGBMClassifier(**_PARAMS).fit(Xtr, y[tr]).predict_proba(Xte)[:, 1]
    per = {}
    for L in LINEAGES:
        m = (lineage == L) & ~np.isnan(oof)
        if m.sum() >= 20:
            per[L] = _auc(y[m], oof[m])
    aucs = [v for v in per.values() if v is not None]
    return round(float(np.mean(aucs)), 4) if aucs else None


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
        if not conc:
            continue
        Xc = df[conc].to_numpy(dtype=float)
        Xcov = df[cov].to_numpy(dtype=float)

        def score(tr, te, resid):
            w = _risk_weights(Xc[tr], y[tr], resid_lineage=lineage[tr] if resid else None)
            return Xc[tr] @ w, Xc[te] @ w

        # (A) raw concordant + cov
        auc_A = _loo_auc(df, lambda tr, te: (np.hstack([Xc[tr], Xcov[tr]]),
                                             np.hstack([Xc[te], Xcov[te]])), y, lineage)
        # (B) C2 risk score alone + cov
        def bB(tr, te):
            s_tr, s_te = score(tr, te, resid=False)
            return np.hstack([s_tr[:, None], Xcov[tr]]), np.hstack([s_te[:, None], Xcov[te]])
        auc_B = _loo_auc(df, bB, y, lineage)
        # (C) raw + C2 score + cov
        def bC(tr, te):
            s_tr, s_te = score(tr, te, resid=False)
            return (np.hstack([Xc[tr], s_tr[:, None], Xcov[tr]]),
                    np.hstack([Xc[te], s_te[:, None], Xcov[te]]))
        auc_C = _loo_auc(df, bC, y, lineage)
        # (D) C1 lineage-residual score alone + cov
        def bD(tr, te):
            s_tr, s_te = score(tr, te, resid=True)
            return np.hstack([s_tr[:, None], Xcov[tr]]), np.hstack([s_te[:, None], Xcov[te]])
        auc_D = _loo_auc(df, bD, y, lineage)

        results[drug] = {
            "n_concordant": len(conc),
            "auc_A_raw_concordant": auc_A,
            "auc_B_c2_score_only": auc_B,
            "auc_C_raw_plus_c2": auc_C,
            "auc_D_c1_residual_score_only": auc_D,
            "delta_C_minus_A": round(auc_C - auc_A, 4) if auc_A and auc_C else None,
        }
        print(f"{drug:<5} conc={len(conc):>3} | A(raw)={auc_A} B(C2only)={auc_B} "
              f"C(raw+C2)={auc_C} D(C1only)={auc_D} | ΔC-A={results[drug]['delta_C_minus_A']:+}")
    manifest = {"experiment": "causal columns as features (C1/C2/C4) — fold-internal risk score, lineage-LOO",
                "feature_sets": {"A": "raw concordant + cov", "B": "C2 risk score only + cov",
                                 "C": "raw + C2 score + cov", "D": "C1 lineage-residual score only + cov"},
                "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description="Causal-columns-as-features (C1/C2/C4) lineage-LOO evaluation.")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/causal_features"))
    args = p.parse_args()
    evaluate(args.marts, args.causal, args.out)


if __name__ == "__main__":
    main()
