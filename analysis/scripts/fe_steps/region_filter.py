"""L0 — drop hypervariable regions from the candidate universe.

PE/PPE, ESX and insertion-sequence/phage genes are the hypervariable,
error-prone-mapping part of the genome. design/20 measured what they contribute:
a PE/PPE-only model scores 0.67-0.83 under random cross-validation and collapses
to roughly chance (0.49-0.61) under lineage-leave-one-out, for every drug, mean
gap 0.22. That is the same signature as synonymous mutations -- a
population-structure proxy, not resistance biology.

Dropping them was measured as neutral-to-positive across the usable and moderate
tiers. The one apparent loss, BDQ -0.056, is diagnostic rather than a cost: BDQ's
inflated de-novo score was partly built on these lineage shortcuts, and removing
them moves it toward its honest concordant value.

Label-free: the region is parsed from the gene name, so this cannot see the
phenotype and runs once, before the mart is shared with any fold.
"""

from __future__ import annotations

import re

import pandas as pd

from analysis.scripts.fe_steps.registry import Step, register

# Parsed from the gene name alone -- no annotation file needed, no download,
# and nothing here can depend on the label.
FAMILIES = {
    "pe_ppe": re.compile(r"^(PE|PPE)", re.I),
    "esx":    re.compile(r"^esx", re.I),
    "is_phage": re.compile(r"^(IS\d|Rv\d+c?_IS)", re.I),
}
DEFAULT_FAMILIES = ("pe_ppe", "esx", "is_phage")


def gene_of(column: str) -> str | None:
    """raw__<gene>_<mutation> -> gene. Returns None for non-mutation columns."""
    if not column.startswith("raw__"):
        return None
    gene, _, _ = column[len("raw__"):].partition("_")
    return gene or None


def run(mart: pd.DataFrame, *, families: str = ",".join(DEFAULT_FAMILIES),
        annotations: pd.DataFrame | None = None, **_) -> tuple[pd.DataFrame, dict]:
    wanted = [f.strip() for f in families.split(",") if f.strip()]
    unknown = [f for f in wanted if f not in FAMILIES]
    if unknown:
        raise ValueError(f"unknown region family {unknown}; known: {sorted(FAMILIES)}")

    hits: dict[str, list[str]] = {f: [] for f in wanted}
    for col in mart.columns:
        gene = gene_of(col)
        if not gene:
            continue
        for f in wanted:
            if FAMILIES[f].match(gene):
                hits[f].append(col)
                break

    # An annotation map, when supplied, catches families the gene name does not
    # spell out -- the curated insertion-sequence and phage list in particular.
    annotated = []
    if annotations is not None and "is_hitchhiker_region" in annotations.columns:
        flagged = set(annotations.loc[annotations["is_hitchhiker_region"].astype(bool),
                                      annotations.columns[0]].astype(str))
        annotated = [c for c in mart.columns
                     if (g := gene_of(c)) and g in flagged
                     and not any(c in v for v in hits.values())]

    drop = sorted({c for v in hits.values() for c in v} | set(annotated))
    kept = mart.drop(columns=drop)

    n_mut = sum(1 for c in mart.columns if c.startswith("raw__"))
    note = {
        "step": "region_filter",
        "families": wanted,
        "dropped": len(drop),
        "dropped_by_family": {f: len(v) for f, v in hits.items()},
        "dropped_by_annotation": len(annotated),
        "mutation_features_before": n_mut,
        "mutation_features_after": n_mut - len(drop),
        "fraction_dropped": round(len(drop) / n_mut, 4) if n_mut else 0.0,
    }
    return kept, note


register(Step(
    name="region_filter",
    layer="universe",
    label_free=True,
    summary="drop PE/PPE, ESX and IS/phage features -- measured population-structure proxies",
    run=run,
))
