"""L0 — keep or drop mutation classes.

design/15 Part D measured what the synonymous class contributes: a
synonymous-only model scores 0.77-0.80 under random cross-validation and loses
0.21 AUC under lineage-leave-one-out, reaching chance for RIF and EMB. Silent
mutations cannot change protein function, so that apparent signal is population
structure and nothing else. Concordance selection already depletes them (odds
ratio 0.14), but dropping them from the universe makes the de-novo arm honest
too, and is one line rather than a fitted model.

The classifier is imported from evaluate_mutation_class rather than reimplemented,
so the filter and the experiment that justified it can never disagree about what
"synonymous" means.

Label-free: the class is parsed from the GARC string.
"""

from __future__ import annotations

import pandas as pd

from analysis.scripts.feature_mart.mutation_class import CLASSES, mut_class
from analysis.scripts.fe_steps.registry import Step, register

def run(mart: pd.DataFrame, *, drop: str = "synonymous", keep: str | None = None,
        **_) -> tuple[pd.DataFrame, dict]:
    if keep and drop:
        drop = ""  # keep wins; naming both is ambiguous rather than additive
    wanted_keep = {c.strip() for c in keep.split(",")} if keep else None
    wanted_drop = {c.strip() for c in drop.split(",") if c.strip()} if drop else set()
    for name in (wanted_keep or set()) | wanted_drop:
        if name not in CLASSES:
            raise ValueError(f"unknown mutation class {name!r}; known: {CLASSES}")

    counts: dict[str, int] = {}
    removed = []
    for col in mart.columns:
        if not col.startswith("raw__"):
            continue
        cls = mut_class(col)
        counts[cls] = counts.get(cls, 0) + 1
        if (wanted_keep is not None and cls not in wanted_keep) or cls in wanted_drop:
            removed.append(col)

    kept = mart.drop(columns=removed)
    n_mut = sum(counts.values())
    note = {
        "step": "class_filter",
        "mode": "keep" if wanted_keep else "drop",
        "classes": sorted(wanted_keep or wanted_drop),
        "class_counts": counts,
        "dropped": len(removed),
        "mutation_features_before": n_mut,
        "mutation_features_after": n_mut - len(removed),
        "fraction_dropped": round(len(removed) / n_mut, 4) if n_mut else 0.0,
    }
    return kept, note


register(Step(
    name="class_filter",
    layer="universe",
    label_free=True,
    summary="keep or drop mutation classes; synonymous is the measured lineage proxy",
    run=run,
))
