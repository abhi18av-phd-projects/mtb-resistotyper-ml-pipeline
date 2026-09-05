"""Homoplasy / convergence feature (the on-thesis transferable signal, design/19 §7).

A mutation that arose *independently* many times across the phylogeny — under drug
pressure — is causal by construction and **lineage-invariant by construction**: it is
exactly what should survive lineage-LOO, and the opposite of a single-clade hitchhiker.
True homoplasy needs ancestral-state reconstruction on a tree; we use a tree-free proxy:
**sublineage spread** — the number of distinct sublineages (and major lineages) carrying
the mutation. A convergent determinant scatters across many sublineages; a lineage
marker stays within one clade despite high prevalence.

Two questions:
  1. **Validation** — do the causally-concordant drivers have *higher* homoplasy than the
     non-selected universe? (ties phylogenetic convergence to the causal selection.)
  2. **Feature/selector** — does a homoplasy-selected feature set generalise under
     lineage-LOO (L2-logistic) comparably to the concordant set, and does adding a
     per-isolate homoplasy-weighted burden help?

In-memory, manifest only.

Run:
    python -m analysis.scripts.feature_mart.evaluate_homoplasy \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --out analysis/results/fulldata_ml/homoplasy
"""

from __future__ import annotations

import argparse
import glob
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.scripts.feature_mart.evaluate_mutation_class import mut_class
from analysis.scripts.feature_mart.train_logistic import _best_C, _fit, _proba, _auc, LINEAGES

COV = ["cov__median_coverage", "cov__tb_breadth"]
AGG = {"raw__n_mutations", "raw__mean_frs", "raw__mean_coverage", "raw__n_minor"}
MIN_CARR = 10


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
        subl = df["cov__sublineage"].astype(str).to_numpy()
        universe = [c for c in df.columns if c.startswith("raw__") and c not in AGG
                    and not c.startswith("raw__gene_") and mut_class(c) != "synonymous"]
        conc = set(pd.read_parquet(cpath).query("qf__causal_concordant")["feature_column"]) & set(universe)

        # per-mutation homoplasy proxy: distinct sublineages / lineages among carriers
        homo = {}
        for c in universe:
            carr = df[c].to_numpy() > 0
            n = int(carr.sum())
            if n >= MIN_CARR:
                homo[c] = {"n_carriers": n, "n_sublineages": int(len(set(subl[carr]))),
                           "n_lineages": int(len(set(lineage[carr])))}
        eligible = list(homo)
        if not eligible or not conc:
            continue

        # (1) validation: homoplasy of concordant vs non-concordant
        cvals = [homo[c]["n_sublineages"] for c in eligible if c in conc]
        nvals = [homo[c]["n_sublineages"] for c in eligible if c not in conc]
        # (2) selector: top-N by homoplasy (N = |concordant|), lineage-LOO vs concordant
        N = len(conc)
        top_homo = sorted(eligible, key=lambda c: -homo[c]["n_sublineages"])[:N]
        cov = [c for c in COV if c in df.columns]
        Xcov = df[cov].to_numpy("float32")
        auc_conc = _loo(np.hstack([df[sorted(conc)].to_numpy("float32"), Xcov]), y, lineage)
        auc_homo = _loo(np.hstack([df[top_homo].to_numpy("float32"), Xcov]), y, lineage)
        # homoplasy-weighted burden added to concordant
        w = np.array([homo[c]["n_sublineages"] for c in eligible], dtype="float32")
        burden = (df[eligible].to_numpy("float32") * w).sum(axis=1, keepdims=True)
        auc_aug = _loo(np.hstack([df[sorted(conc)].to_numpy("float32"), burden, Xcov]), y, lineage)

        results[drug] = {
            "n_concordant": len(conc), "n_eligible": len(eligible),
            "homoplasy_concordant_mean_sublineages": round(float(np.mean(cvals)), 2) if cvals else None,
            "homoplasy_nonconcordant_mean_sublineages": round(float(np.mean(nvals)), 2) if nvals else None,
            "auc_concordant": auc_conc, "auc_homoplasy_selected": auc_homo,
            "auc_concordant_plus_homoburden": auc_aug,
            "delta_homoburden": round(auc_aug - auc_conc, 4) if auc_conc and auc_aug else None,
        }
        r = results[drug]
        print(f"{drug:<5} homoplasy(subl): conc={r['homoplasy_concordant_mean_sublineages']} "
              f"non-conc={r['homoplasy_nonconcordant_mean_sublineages']} | "
              f"AUC conc={auc_conc} homo-selected={auc_homo} +homoburden={auc_aug} ({r['delta_homoburden']:+})")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "homoplasy/convergence (sublineage-spread proxy) — validation + selector + burden",
         "note": "tree-free proxy; true homoplasy needs ancestral-state reconstruction",
         "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return results


def main():
    p = argparse.ArgumentParser(description="Homoplasy/convergence feature (design/19 §7).")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/homoplasy"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.out)


if __name__ == "__main__":
    main()
