"""AA property-change burden — a cheap transferable "radicalness" proxy (design/19 §7).

A poor-man's ESM: for each missense variant (wt→mut residue, parsed from the GARC), how
radical is the substitution physicochemically? Encoded from confident residue scales —
Kyte-Doolittle hydrophobicity, net charge, side-chain volume — as
  radicalness(wt,mut) = |Δhydro| + 2·|Δcharge| + |Δvolume|/50 .
Per isolate we sum the radicalness of carried missense variants (in the concordant
determinant genes) into `burden_radical` and add it to the concordant + L2-logistic
model. Transferable (a property of the substitution, lineage-independent) and additive.
Tested vs the concordant baseline under lineage-LOO; must beat it to be adopted.

Run:
    python -m analysis.scripts.feature_mart.evaluate_aa_property \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --out analysis/results/fulldata_ml/aa_property
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

from analysis.scripts.feature_mart.train_logistic import _best_C, _fit, _proba, _auc, LINEAGES

COV = ["cov__median_coverage", "cov__tb_breadth"]
HYDRO = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
         "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
         "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}
CHARGE = {"D": -1, "E": -1, "K": 1, "R": 1, "H": 0.5}
VOL = {"A": 88.6, "R": 173.4, "N": 114.1, "D": 111.1, "C": 108.5, "Q": 143.8, "E": 138.4,
       "G": 60.1, "H": 153.2, "I": 166.7, "L": 166.7, "K": 168.6, "M": 162.9, "F": 189.9,
       "P": 112.7, "S": 89.0, "T": 116.1, "W": 227.8, "Y": 193.6, "V": 140.0}
_MISSENSE = re.compile(r"^([A-Z])(\d+)([A-Z])$")


def radicalness(col: str) -> float:
    mut = col.removeprefix("raw__").partition("_")[2]
    m = _MISSENSE.match(mut)
    if not m:
        return 0.0
    wt, _, mt = m.groups()
    if wt not in HYDRO or mt not in HYDRO or wt == mt:
        return 0.0
    dh = abs(HYDRO[wt] - HYDRO[mt])
    dc = abs(CHARGE.get(wt, 0) - CHARGE.get(mt, 0))
    dv = abs(VOL[wt] - VOL[mt]) / 50.0
    return dh + 2 * dc + dv


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
    per = {L: _auc(y[(lineage == L) & ~np.isnan(oof)], oof[(lineage == L) & ~np.isnan(oof)])
           for L in LINEAGES if ((lineage == L) & ~np.isnan(oof)).sum() >= 20}
    aucs = [v for v in per.values() if v is not None]
    return round(float(np.mean(aucs)), 4) if aucs else None


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
        if not conc:
            continue
        cov = [c for c in COV if c in df.columns]
        # radicalness-weighted missense burden over the concordant genes' universe
        det_genes = {c.removeprefix("raw__").split("_")[0] for c in conc}
        miss = [c for c in df.columns if c.startswith("raw__") and c.removeprefix("raw__").split("_")[0] in det_genes
                and radicalness(c) > 0]
        w = np.array([radicalness(c) for c in miss], dtype="float32")
        burden = (df[miss].to_numpy("float32") * w).sum(axis=1, keepdims=True) if miss else np.zeros((len(df), 1), "float32")
        Xcov = df[cov].to_numpy("float32")
        Xc = df[conc].to_numpy("float32")
        A = _loo(np.hstack([Xc, Xcov]), y, lineage)
        B = _loo(np.hstack([Xc, burden, Xcov]), y, lineage)
        results[drug] = {"n_concordant": len(conc), "n_missense_in_det_genes": len(miss),
                         "auc_concordant": A, "auc_plus_radicalness": B,
                         "delta": round(B - A, 4) if A and B else None}
        print(f"{drug:<5} conc={len(conc):>3} miss={len(miss):>3} | A={A} B={B} (Δ {results[drug]['delta']:+})")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "AA property-change (radicalness) burden vs concordant, L2-logistic lineage-LOO",
         "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return results


def main():
    p = argparse.ArgumentParser(description="AA property-change radicalness burden (design/19 §7).")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/aa_property"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.out)


if __name__ == "__main__":
    main()
