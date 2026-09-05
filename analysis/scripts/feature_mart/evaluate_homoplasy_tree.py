"""Proper phylogenetic homoplasy via parsimony on the SNP-barcode tree (design/15 Part D).

The tree-free spread proxy failed (confounded by prevalence). The fix: count a mutation's
**independent origins** on the phylogeny — prevalence-independent by construction. The
CRyPTIC sublineage labels ARE the SNP-barcode hierarchy (`lineage4.1.2.1` = the nested
tree), so we build that tree from the labels and run **Fitch parsimony**, rooted at the
ancestral wild-type (absent), counting 0→1 transitions = independent gains = homoplasy.
A convergent determinant (arising repeatedly under drug pressure) has many gains; a
single-origin lineage marker has one — regardless of how prevalent it is.

Per mutation: presence on a sublineage tip = carried by ≥ MIN_CARR of that sublineage's
isolates. Then: (1) validation — do concordant drivers have MORE independent origins
than the non-selected universe? (the proxy got this backwards); (2) selector — does a
homoplasy-selected set generalise under lineage-LOO?

In-memory, manifest only.

Run:
    python -m analysis.scripts.feature_mart.evaluate_homoplasy_tree \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --out analysis/results/fulldata_ml/homoplasy_tree
"""

from __future__ import annotations

import argparse
import glob
import json
import time
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd

from analysis.scripts.feature_mart.evaluate_mutation_class import mut_class
from analysis.scripts.feature_mart.train_logistic import _best_C, _fit, _proba, _auc, LINEAGES

COV = ["cov__median_coverage", "cov__tb_breadth"]
AGG = {"raw__n_mutations", "raw__mean_frs", "raw__mean_coverage", "raw__n_minor"}
MIN_CARR = 10          # min carriers overall to score a mutation
MIN_TIP = 2            # min carriers within a sublineage to call the mutation present there
ROOT = "__root__"


def _build_tree(labels):
    """Barcode hierarchy from dotted sublineage labels -> children map + tip nodes.
    Each observed sublineage gets a `$self` leaf so its own isolates carry a state."""
    children = defaultdict(set)
    tips = {}
    for lab in labels:
        parts = lab.split(".")
        nodes = [".".join(parts[: i + 1]) for i in range(len(parts))]
        prev = ROOT
        for n in nodes:
            children[prev].add(n)
            prev = n
        leaf = lab + "$self"
        children[lab].add(leaf)
        tips[lab] = leaf
    return children, tips


def _independent_origins(children, tip_state):
    """Fitch parsimony rooted at absent (0); count 0->1 gains = independent origins."""
    up = {}

    def post(node):
        kids = children.get(node)
        if not kids:                      # leaf
            up[node] = {tip_state.get(node, 0)}
            return up[node]
        sets = [post(k) for k in kids]
        inter = set.intersection(*sets)
        up[node] = inter if inter else set.union(*sets)
        return up[node]

    post(ROOT)
    gains = [0]
    root_state = 0 if 0 in up[ROOT] else next(iter(up[ROOT]))

    def pre(node, parent_state):
        s = up[node]
        state = parent_state if parent_state in s else min(s)
        if parent_state == 0 and state == 1:
            gains[0] += 1
        for k in children.get(node, ()):
            pre(k, state)

    pre(ROOT, root_state)
    return gains[0]


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
        clean = np.array([s if (s.startswith("lineage") and "/" not in s) else "" for s in subl])
        keep = clean != ""
        labels = sorted(set(clean[keep]))
        children, tips = _build_tree(labels)

        universe = [c for c in df.columns if c.startswith("raw__") and c not in AGG
                    and not c.startswith("raw__gene_") and mut_class(c) != "synonymous"]
        conc = set(pd.read_parquet(cpath).query("qf__causal_concordant")["feature_column"]) & set(universe)
        # precompute per-sublineage carrier counts, vectorised
        homo = {}
        for c in universe:
            carr = df[c].to_numpy() > 0
            if int((carr & keep).sum()) < MIN_CARR:
                continue
            counts = pd.Series(clean[carr & keep]).value_counts()
            tip_state = {tips[lab]: 1 for lab, n in counts.items() if lab in tips and n >= MIN_TIP}
            homo[c] = _independent_origins(children, tip_state)
        eligible = list(homo)
        if not eligible or not conc:
            continue

        cvals = [homo[c] for c in eligible if c in conc]
        nvals = [homo[c] for c in eligible if c not in conc]
        N = len(conc)
        top_homo = sorted(eligible, key=lambda c: -homo[c])[:N]
        cov = [c for c in COV if c in df.columns]
        Xcov = df[cov].to_numpy("float32")
        auc_conc = _loo(np.hstack([df[sorted(conc)].to_numpy("float32"), Xcov]), y, lineage)
        auc_homo = _loo(np.hstack([df[top_homo].to_numpy("float32"), Xcov]), y, lineage)
        results[drug] = {
            "n_concordant": len(conc), "n_eligible": len(eligible),
            "origins_concordant_mean": round(float(np.mean(cvals)), 2) if cvals else None,
            "origins_nonconcordant_mean": round(float(np.mean(nvals)), 2) if nvals else None,
            "origins_concordant_median": float(np.median(cvals)) if cvals else None,
            "origins_nonconcordant_median": float(np.median(nvals)) if nvals else None,
            "auc_concordant": auc_conc, "auc_homoplasy_selected": auc_homo,
        }
        r = results[drug]
        print(f"{drug:<5} origins: conc={r['origins_concordant_mean']} (med {r['origins_concordant_median']}) "
              f"non-conc={r['origins_nonconcordant_mean']} (med {r['origins_nonconcordant_median']}) | "
              f"AUC conc={auc_conc} homo-sel={auc_homo}")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "phylogenetic homoplasy (Fitch independent-origins on SNP-barcode tree)",
         "params": {"MIN_CARR": MIN_CARR, "MIN_TIP": MIN_TIP},
         "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return results


def main():
    p = argparse.ArgumentParser(description="Phylogenetic homoplasy (Fitch parsimony on barcode tree).")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/homoplasy_tree"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.out)


if __name__ == "__main__":
    main()
