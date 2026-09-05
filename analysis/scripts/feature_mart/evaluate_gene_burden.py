"""Step 2 — function-aware gene-burden / LoF collapsing features.

For the data-poor drugs resistance is a *union of many individually-rare variants in
one gene* (BDQ/CFZ `Rv0678`, LZD `rrl`/`rplC`, DLM's six activation genes). No single
variant recurs enough to be selected, so the concordant set misses them. The standard
rare-variant move — and the literature's lead (design/19 §7-C, Nat Commun 2026
rare-variant augmentation) — is to **collapse** per determinant gene:

  burden_func__<gene> = # functional variants carried (drops synonymous — E3 showed
                        those are lineage hitchhikers)
  burden_lof__<gene>  = # loss-of-function variants (indel/nonsense/frameshift), the
                        class most likely causal for efflux/activation LoF mechanisms

Leakage-free (genotype + GARC class, deterministic). Evaluated on the deployed model
(L2-logistic, C tuned by lineage-LOO), honest lineage-LOO, three feature sets:
  (A) concordant + cov            (the Step-1 baseline)
  (B) concordant + burden + cov   (augmented)
  (C) burden + cov only           (does collapsing alone carry signal?)

Focus: the data-poor / hard tier (BDQ/CFZ/LZD/DLM), but run all determinant drugs.
In-memory, manifest only.

Run:
    python -m analysis.scripts.feature_mart.evaluate_gene_burden \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --out analysis/results/fulldata_ml/gene_burden
"""

from __future__ import annotations

import argparse
import glob
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.scripts.feature_mart.build_mart import DETERMINANT_GENES
from analysis.scripts.feature_mart.evaluate_mutation_class import mut_class
from analysis.scripts.feature_mart.train_logistic import _best_C, _fit, _proba, _auc, LINEAGES, C_GRID  # noqa

COV = ["cov__median_coverage", "cov__tb_breadth"]
LOF_CLASSES = {"indel", "nonsense_Z", "minor_indel"}


def _gene_of(col: str) -> str:
    return col.removeprefix("raw__").partition("_")[0]


def _burden(df: pd.DataFrame, drug: str) -> pd.DataFrame:
    """Per-determinant-gene functional + LoF burden columns (in-memory)."""
    genes = DETERMINANT_GENES.get(drug, set())
    raw = [c for c in df.columns if c.startswith("raw__") and not c.startswith("raw__gene_")
           and _gene_of(c) in genes]
    out = {}
    for g in sorted(genes):
        gcols = [c for c in raw if _gene_of(c) == g]
        if not gcols:
            continue
        func = [c for c in gcols if mut_class(c) != "synonymous"]
        lof = [c for c in gcols if mut_class(c) in LOF_CLASSES]
        if func:
            out[f"burden_func__{g}"] = df[func].sum(axis=1).astype("float32")
        if lof:
            out[f"burden_lof__{g}"] = df[lof].sum(axis=1).astype("float32")
    return pd.DataFrame(out, index=df.index)


def _loo(X, y, lineage):
    """C tuned by lineage-LOO (applied equally to every arm → fair delta), honest AUC."""
    C = _best_C(X, y, lineage, LINEAGES)
    oof = np.full(len(y), np.nan)
    for L in LINEAGES:
        te = lineage == L
        tr = np.isin(lineage, [x for x in LINEAGES if x != L])
        if te.sum() < 20 or len(np.unique(y[tr])) < 2:
            continue
        sc, lr = _fit(X[tr], y[tr], C)
        oof[te] = _proba(sc, lr, X[te])
    per = {L: _auc(y[(lineage == L) & ~np.isnan(oof)], oof[(lineage == L) & ~np.isnan(oof)])
           for L in LINEAGES if ((lineage == L) & ~np.isnan(oof)).sum() >= 20}
    aucs = [v for v in per.values() if v is not None]
    return {"auc": round(float(np.mean(aucs)), 4) if aucs else None, "C": C,
            "per_lineage": {k: (round(v, 3) if v else None) for k, v in per.items()}}


def evaluate(marts_dir, causal_dir, out_dir):
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
        burden = _burden(df, drug)
        Xcov = df[cov].to_numpy("float32")
        B = burden.to_numpy("float32") if burden.shape[1] else np.empty((len(df), 0), "float32")

        arms = {}
        if conc:
            arms["A_concordant"] = _loo(np.hstack([df[conc].to_numpy("float32"), Xcov]), y, lineage)
            if B.shape[1]:
                arms["B_concordant_plus_burden"] = _loo(np.hstack([df[conc].to_numpy("float32"), B, Xcov]), y, lineage)
        if B.shape[1]:
            arms["C_burden_only"] = _loo(np.hstack([B, Xcov]), y, lineage)
        a, b = arms.get("A_concordant", {}).get("auc"), arms.get("B_concordant_plus_burden", {}).get("auc")
        results[drug] = {"n_concordant": len(conc), "burden_features": list(burden.columns),
                         "arms": arms, "delta_B_minus_A": round(b - a, 4) if a and b else None}
        print(f"{drug:<5} conc={len(conc):>3} burden={burden.shape[1]:>2} | "
              f"A={a} B={b} C={arms.get('C_burden_only',{}).get('auc')} | Δ(B-A)={results[drug]['delta_B_minus_A']}")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "step2 — function-aware gene-burden/LoF collapsing (L2-logistic, lineage-LOO)",
         "lof_classes": sorted(LOF_CLASSES), "per_drug": results,
         "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return results


def main():
    p = argparse.ArgumentParser(description="Step 2: gene-burden / LoF collapsing features.")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/gene_burden"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.out)


if __name__ == "__main__":
    main()
