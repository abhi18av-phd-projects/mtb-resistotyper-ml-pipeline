"""Annotation-derived FE — functional-category burden + curated resistance-gene burden.

Uses the Mycobrowser v5 `gene_annotations` map (design/20) to add the gene-group /
functional-context signal the ML-on-MTB attention & graph models exploit, but as
explicit additive features on the deployed L2-logistic (our thesis: simple model +
transferable features). Two annotation features:

- **functional-category burden** `burden_cat__<category>` — count of non-synonymous
  variants (synonymous dropped, per E3) the isolate carries in genes of each TubercuList
  category (cell wall, information pathways, lipid metabolism, …). A structured
  gene-group feature.
- **curated resistance-gene burden** `burden_resgene__<drug>` — count of variants in
  genes the annotation flags as resistance-associated for this drug (`Drug Resistance
  Mutations`), an annotation-driven analogue of determinant seeding. Only for drugs the
  Mycobrowser flags cover (RIF/INH/EMB/ETH/FLQ/AMI-SM; not the newer drugs).

Feature sets under lineage-LOO (L2-logistic, C tuned by lineage-LOO — applied equally,
so deltas are fair): (A) concordant+cov · (B) +category-burden · (C) +resgene-burden.
Must beat (A) to be adopted, same bar as every FE lead. In-memory, manifest only.

Run:
    python -m analysis.scripts.feature_mart.evaluate_annotation_features \
        --marts analysis/results/fulldata_ml/marts --causal analysis/results/fulldata_ml/causal \
        --annotations analysis/results/fulldata_ml/gene_annotations.parquet \
        --out analysis/results/fulldata_ml/annotation_features
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
# our drug code -> annotation Drug-Resistance-Mutations codes
ANNO_CODE = {"RIF": {"RIF"}, "INH": {"INH"}, "EMB": {"EMB"}, "MXF": {"FLQ"}, "LEV": {"FLQ"},
             "ETH": {"ETH"}, "AMI": {"AMI", "SM"}, "KAN": {"AMI", "SM"}}


def _gene_of(col: str) -> str:
    return col.removeprefix("raw__").split("_")[0]


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
    return {"auc": round(float(np.mean(aucs)), 4) if aucs else None, "C": C,
            "per_lineage": {k: (round(v, 3) if v else None) for k, v in per.items()}}


def evaluate(marts_dir, causal_dir, ann_path, out_dir):
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    ann = pd.read_parquet(ann_path)
    gene2cat, gene2res = {}, {}
    for _, r in ann.iterrows():
        cat = r["functional_category"]
        res = set(x.strip() for x in r["drug_resistance"].split(",")) if r["drug_resistance"] else set()
        for key in (r["gene"], r["locus"]):
            if key:
                if cat:
                    gene2cat[key] = cat
                if res:
                    gene2res[key] = res
    categories = sorted(set(gene2cat.values()))

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
        universe = [c for c in df.columns if c.startswith("raw__") and c not in AGG
                    and not c.startswith("raw__gene_") and mut_class(c) != "synonymous"]  # functional only
        Xcov = df[cov].to_numpy("float32")

        # functional-category burden (non-synonymous count per category)
        cat_cols, cat_names = [], []
        for cat in categories:
            cols = [c for c in universe if gene2cat.get(_gene_of(c)) == cat]
            if cols:
                cat_cols.append(df[cols].sum(axis=1).to_numpy("float32"))
                cat_names.append(f"burden_cat__{cat}")
        Bcat = np.column_stack(cat_cols) if cat_cols else np.empty((len(df), 0), "float32")

        # curated resistance-gene burden for this drug
        codes = ANNO_CODE.get(drug, set())
        res_cols = [c for c in universe if codes & gene2res.get(_gene_of(c), set())]
        Bres = df[res_cols].sum(axis=1).to_numpy("float32").reshape(-1, 1) if res_cols else np.empty((len(df), 0), "float32")

        Xc = df[conc].to_numpy("float32")
        A = _loo(np.hstack([Xc, Xcov]), y, lineage)
        B = _loo(np.hstack([Xc, Bcat, Xcov]), y, lineage) if Bcat.shape[1] else None
        Cc = _loo(np.hstack([Xc, Bres, Xcov]), y, lineage) if Bres.shape[1] else None
        results[drug] = {
            "n_concordant": len(conc), "n_categories": Bcat.shape[1],
            "n_resgene_features": Bres.shape[1], "res_genes": sorted({_gene_of(c) for c in res_cols}),
            "A_concordant": A["auc"], "B_plus_category": B["auc"] if B else None,
            "C_plus_resgene": Cc["auc"] if Cc else None,
            "delta_B": round(B["auc"] - A["auc"], 4) if B and A["auc"] and B["auc"] else None,
            "delta_C": round(Cc["auc"] - A["auc"], 4) if Cc and A["auc"] and Cc["auc"] else None,
            "B_per_lineage": B["per_lineage"] if B else None,
        }
        print(f"{drug:<5} conc={len(conc):>3} cats={Bcat.shape[1]} res={Bres.shape[1]} | "
              f"A={A['auc']} B={results[drug]['B_plus_category']}({results[drug]['delta_B']}) "
              f"C={results[drug]['C_plus_resgene']}({results[drug]['delta_C']})")
    (out_dir / "manifest.json").write_text(json.dumps(
        {"experiment": "annotation FE — functional-category + curated-resistance-gene burden (L2-logistic, lineage-LOO)",
         "categories": categories, "per_drug": results, "wall_time_seconds": round(time.time() - t0, 2)}, indent=2))
    return results


def main():
    p = argparse.ArgumentParser(description="Annotation-derived FE (category + resistance-gene burden).")
    p.add_argument("--marts", type=Path, default=Path("analysis/results/fulldata_ml/marts"))
    p.add_argument("--causal", type=Path, default=Path("analysis/results/fulldata_ml/causal"))
    p.add_argument("--annotations", type=Path, default=Path("analysis/results/fulldata_ml/gene_annotations.parquet"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/annotation_features"))
    a = p.parse_args()
    evaluate(a.marts, a.causal, a.annotations, a.out)


if __name__ == "__main__":
    main()
