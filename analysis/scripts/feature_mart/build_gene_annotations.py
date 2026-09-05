"""Build the `gene_annotations` map from the Mycobrowser v5 H37Rv GFF (design/20).

Source: `Mycobacterium_tuberculosis_H37Rv_gff_v5.gff` (Mycobrowser Release v5,
NC-000962-3 / H37Rv), located in the sibling MAGMA pipeline resources. It carries, per
locus: gene Name, Functional_Category (the TubercuList categories), Product/Function/
Comments (DeJesus essentiality references), a curated `Drug Resistance Mutations` flag,
PFAM domain, GO, coordinates, pseudogene flag. This is the source design/02 named.

We derive a compact, committed map (join key = `gene` Name or `locus` Rv-number, which
match the CRyPTIC `mutations.GENE` strings) rather than vendoring the 4.4 MB GFF — it is
a standard public release. 99% of the mart universe genes map.

Run:
    python -m analysis.scripts.feature_mart.build_gene_annotations \
        --gff /path/to/Mycobacterium_tuberculosis_H37Rv_gff_v5.gff \
        --out analysis/results/fulldata_ml/gene_annotations.parquet
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import pandas as pd

# hypervariable / repetitive categories that behave as lineage hitchhikers (design/20)
HITCHHIKER_CATEGORIES = {"PE/PPE", "insertion seqs and phages"}
DEFAULT_GFF = ("/Users/abhi/projects/PHD-pub-magma-pipeline/analysis/pipelines/"
               "torch-magma/_resources/exp_gatk4_sv/Mycobacterium_tuberculosis_H37Rv_gff_v5.gff")


def build(gff: Path, out: Path) -> pd.DataFrame:
    rows = []
    for line in Path(gff).read_text().splitlines():
        if line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 9:
            continue
        attr = dict(kv.split("=", 1) for kv in f[8].split(";") if "=" in kv)
        text = f"{attr.get('Product', '')} {attr.get('Comments', '')}"
        cat = attr.get("Functional_Category", "")
        rows.append({
            "locus": attr.get("Locus", ""), "gene": attr.get("Name", ""),
            "feature_type": f[2], "functional_category": cat,
            "is_hitchhiker_region": cat in HITCHHIKER_CATEGORIES,
            "drug_resistance": attr.get("Drug Resistance Mutations", ""),
            "is_pseudogene": attr.get("Is_Pseudogene", "") == "Yes",
            "essential_hint": bool(re.search(r"\bessential\b", text, re.I)),
            "pfam": attr.get("PFAM", ""), "start": int(f[3]), "end": int(f[4]), "strand": f[6],
        })
    ann = pd.DataFrame(rows)
    out.parent.mkdir(parents=True, exist_ok=True)
    ann.to_parquet(out, index=False)
    return ann


def main():
    p = argparse.ArgumentParser(description="Build gene_annotations map from Mycobrowser v5 GFF.")
    p.add_argument("--gff", type=Path, default=Path(DEFAULT_GFF))
    p.add_argument("--out", type=Path, default=Path("analysis/results/fulldata_ml/gene_annotations.parquet"))
    a = p.parse_args()
    ann = build(a.gff, a.out)
    print(f"{len(ann)} features -> {a.out}")
    print(ann[ann.functional_category != ""].functional_category.value_counts().to_string())


if __name__ == "__main__":
    main()
