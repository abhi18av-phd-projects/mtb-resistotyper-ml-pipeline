"""ESM zero-shot variant-effect (A2) — validation + feature test (design/19 §9).

Uses precomputed ESM-2 wild-type-marginal LLR scores (esm_scores.parquet; negative =
evolutionarily disfavoured / likely disruptive). Two questions, mirroring the homoplasy
validation:
  1. **Validation** — do the causally-concordant missense determinants get MORE negative
     (more disruptive) ESM scores than the non-selected universe? (ties the protein
     language model to the causal selection.)
  2. **Feature** — does an ESM-weighted missense burden (Σ over carried variants of
     −LLR, so disruptive variants count more) added to the concordant + L2-logistic beat
     the concordant baseline under lineage-LOO?

Protein-only: applies to protein-coding determinants (BDQ `Rv0678`/`atpE`, DLM
`ddn`/`fgd1`/`fbi`, and the protein first-line drugs); not `rrs`/`rrl`/promoters.
In-memory, manifest only.

Run:
    python -m analysis.scripts.feature_mart.evaluate_esm_features \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --scores /tmp/esmwork/esm_scores.parquet --out analysis/results/fulldata_ml/esm_features
"""

from __future__ import annotations

import argparse
import glob
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.scripts.feature_mart.train_logistic import _best_C, _fit, _proba, _auc, LINEAGES

COV = ["cov__median_coverage", "cov__tb_breadth"]


def _loo(X, y, lineage):
    C = _best_C(X, y, lineage, LINEAGES)
    oof = np.full(len(y), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        sc, lr = _fit(X[tr], y[tr], C)
        oof[te] = _proba(sc, lr, X[te])
    per = [_auc(y[(lineage == L) & ~np.isnan(oof)], oof[(lineage == L) & ~np.isnan(oof)])
           for L in LINEAGES if ((lineage == L) & ~np.isnan(oof)).sum() >= 20]
    per = [v for v in per if v is not None]
    return round(float(np.mean(per)), 4) if per else None


def evaluate(marts_dir, causal_dir, scores_path, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    esm = pd.read_parquet(scores_path).set_index("feature")["esm_llr"].to_dict()
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
        if not conc:
            continue
        cov = [c for c in COV if c in df.columns]
        scored = [c for c in df.columns if c.startswith("raw__") and c in esm]  # missense with an ESM score
        # (1) validation: ESM LLR of concordant vs non-concordant scored variants
        cvals = [esm[c] for c in scored if c in conc]
        nvals = [esm[c] for c in scored if c not in conc]
        # (2) feature: ESM-weighted missense burden (disruptive = -LLR) over determinant genes
        det_genes = {c.removeprefix("raw__").split("_")[0] for c in conc}
        burden_cols = [c for c in scored if c.removeprefix("raw__").split("_")[0] in det_genes]
        w = np.array([-esm[c] for c in burden_cols], dtype="float32")   # more disruptive -> larger
        burden = (df[burden_cols].to_numpy("float32") * w).sum(axis=1, keepdims=True) if burden_cols else np.zeros((len(df), 1), "float32")
        Xcov = df[cov].to_numpy("float32")
        Xc = df[conc].to_numpy("float32")
        A = _loo(np.hstack([Xc, Xcov]), y, lineage)
        B = _loo(np.hstack([Xc, burden, Xcov]), y, lineage)
        results[drug] = {
            "n_concordant": len(conc), "n_scored_missense": len(scored),
            "n_concordant_scored": len(cvals), "n_burden_cols": len(burden_cols),
            "esm_llr_concordant_mean": round(float(np.mean(cvals)), 3) if cvals else None,
            "esm_llr_nonconcordant_mean": round(float(np.mean(nvals)), 3) if nvals else None,
            "auc_concordant": A, "auc_plus_esm_burden": B,
            "delta": round(B - A, 4) if A and B else None,
        }
        r = results[drug]
        print(f"{drug:<5} scored={len(scored):>3} concScored={len(cvals):>2} | "
              f"ESM_LLR conc={r['esm_llr_concordant_mean']} non-conc={r['esm_llr_nonconcordant_mean']} | "
              f"AUC A={A} B={B} (Δ {r['delta']})")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "ESM zero-shot variant-effect (A2) — validation + ESM-weighted burden",
         "model": "esm2_t12_35M_UR50D, wild-type-marginal LLR",
         "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return results


def main():
    p = argparse.ArgumentParser(description="ESM zero-shot variant-effect FE (A2).")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--scores", type=Path, default=Path("/tmp/esmwork/esm_scores.parquet"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/esm_features"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.scores, a.out)


if __name__ == "__main__":
    main()
