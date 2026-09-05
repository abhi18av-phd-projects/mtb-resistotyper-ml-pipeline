#!/usr/bin/env python3
"""Assemble a per-drug operating range from the fold-internal evaluations.

Each input manifest comes from one fold, in which concordance selection was run
WITHOUT a particular lineage. Only that lineage's held-out AUC is honest for
that fold: the other lineages in the same manifest were visible to the selection
and their numbers carry the leakage this pipeline exists to remove. So the
report takes one number from each fold -- the AUC on the lineage that fold's
selection never saw -- and averages those.

Selecting on all data and then holding a lineage out for scoring gives a higher
number, which is why the difference between the two is worth reporting rather
than quietly correcting.

Input manifests are named fold_<lineage>.json so the held-out lineage travels
with the file.
"""
from __future__ import annotations

import argparse
import json
import re
import statistics
from pathlib import Path

TIERS = [(0.88, "usable"), (0.80, "moderate"), (0.0, "at chance")]
FOLD_RE = re.compile(r"fold_(?P<lineage>[A-Za-z0-9]+)\.json$")


def tier_for(auc: float) -> str:
    return next(name for floor, name in TIERS if auc >= floor)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--drug", required=True)
    ap.add_argument("--folds", nargs="+", type=Path, required=True)
    ap.add_argument("--feature-set", default="concordant")
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    per_lineage: dict[str, float] = {}
    n_resistant: dict[str, int] = {}
    randoms: list[float] = []
    optimistic: list[float] = []
    skipped: list[str] = []

    for path in sorted(a.folds):
        m = FOLD_RE.search(path.name)
        if not m:
            skipped.append(f"{path.name}: filename does not name a held-out lineage")
            continue
        lineage = m.group("lineage")

        fs = json.loads(path.read_text()).get("results", {}).get(a.feature_set)
        if not fs:
            skipped.append(f"{path.name}: no '{a.feature_set}' feature set")
            continue

        group = (fs.get("lineage_loo", {}).get("per_group") or {}).get(lineage)
        if not group or group.get("auc") is None:
            skipped.append(f"{path.name}: no held-out result for {lineage}")
            continue

        per_lineage[lineage] = group["auc"]
        n_resistant[lineage] = group.get("n_R")

        # context, not the headline
        if (r := fs.get("random", {}).get("overall_auc")) is not None:
            randoms.append(r)
        if (o := fs.get("lineage_loo", {}).get("mean_group_auc")) is not None:
            optimistic.append(o)

    aucs = list(per_lineage.values())
    mean = statistics.fmean(aucs) if aucs else None

    out = {
        "drug": a.drug,
        "feature_set": a.feature_set,
        "selection": "fold-internal",
        "primary_metric": "lineage-leave-one-out, concordance selection refitted within each fold",
        "auc": None if mean is None else round(mean, 4),
        "tier": None if mean is None else tier_for(mean),
        "per_lineage_auc": {k: round(v, 4) for k, v in sorted(per_lineage.items())},
        "n_resistant_held_out": n_resistant,
        "stability_std": round(statistics.pstdev(aucs), 4) if len(aucs) > 1 else None,
        "n_folds": len(per_lineage),
        "skipped": skipped,
    }
    if randoms:
        out["random_cv_auc"] = round(statistics.fmean(randoms), 4)
        if mean is not None:
            out["inflation_vs_random"] = round(out["random_cv_auc"] - mean, 4)
    if optimistic and mean is not None:
        # What the same folds report when every lineage in the manifest is
        # counted, including the ones the selection could see. The gap is the
        # cost of selection leakage, and it is the number this pipeline was
        # built to measure.
        out["auc_selection_on_all_lineages"] = round(statistics.fmean(optimistic), 4)
        out["selection_leakage"] = round(out["auc_selection_on_all_lineages"] - mean, 4)

    a.out.write_text(json.dumps(out, indent=2) + "\n")
    print(json.dumps(out, indent=2))


if __name__ == "__main__":
    main()
