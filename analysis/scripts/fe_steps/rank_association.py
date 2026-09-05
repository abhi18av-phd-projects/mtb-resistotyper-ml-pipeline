"""L4 — rank features by association with the phenotype, inside a fold.

This is where chi-square and mutual-information ranking belong. Both score a
feature against the outcome, so running them at mart-build time would choose the
column universe while looking at every lineage, including the one held out for
testing. The effect is milder than choosing the model's features that way -- it
only decides which columns exist -- but it is the same error, and a paper whose
argument is that unquantified leakage inflates AUC cannot carry one.

Passing held_out_lineage removes those rows before anything is counted, so the
ranking sees only the fold's training data. Omitting it is the deployment refit,
where using all data is correct; the note records which happened.
"""

from __future__ import annotations

import math

import pandas as pd

from analysis.scripts.fe_steps.registry import Step, register

LINEAGE_COL = "cov__lineage_raw"
LABEL_COL = "y__binary"


def _counts(x: pd.Series, y: pd.Series) -> tuple[int, int, int, int]:
    carried = x.astype(float) > 0
    a = int((carried & (y == 1)).sum())   # carried, resistant
    b = int((carried & (y == 0)).sum())   # carried, susceptible
    c = int((~carried & (y == 1)).sum())
    d = int((~carried & (y == 0)).sum())
    return a, b, c, d


def _chi2(a: int, b: int, c: int, d: int) -> float:
    n = a + b + c + d
    denom = (a + b) * (c + d) * (a + c) * (b + d)
    return 0.0 if denom == 0 else n * (a * d - b * c) ** 2 / denom


def _mutual_information(a: int, b: int, c: int, d: int) -> float:
    n = a + b + c + d
    if not n:
        return 0.0
    total = 0.0
    for joint, row, col in ((a, a + b, a + c), (b, a + b, b + d),
                            (c, c + d, a + c), (d, c + d, b + d)):
        if joint and row and col:
            total += (joint / n) * math.log((joint / n) / ((row / n) * (col / n)))
    return total


SCORERS = {"chi2": _chi2, "mi": _mutual_information}


def run(mart: pd.DataFrame, *, method: str = "chi2", top_n: int = 1000,
        held_out_lineage: str | None = None, **_) -> tuple[pd.DataFrame, dict]:
    if method not in SCORERS:
        raise ValueError(f"unknown method {method!r}; known: {sorted(SCORERS)}")
    if LABEL_COL not in mart.columns:
        raise ValueError(f"mart has no {LABEL_COL}; cannot rank by association")

    train = mart
    if held_out_lineage:
        if LINEAGE_COL not in mart.columns:
            raise ValueError(f"mart has no {LINEAGE_COL}; cannot hold a lineage out")
        train = mart.loc[mart[LINEAGE_COL] != held_out_lineage]
        if train.empty:
            raise ValueError(f"holding out {held_out_lineage} left no rows")

    y = train[LABEL_COL]
    scorer = SCORERS[method]
    features = [c for c in mart.columns if c.startswith("raw__")]
    scored = sorted(
        ((c, scorer(*_counts(train[c], y))) for c in features),
        key=lambda kv: -kv[1],
    )
    keep = {c for c, _ in scored[:top_n]}
    dropped = [c for c in features if c not in keep]
    kept = mart.drop(columns=dropped)

    note = {
        "step": "rank_association",
        "method": method,
        "top_n": top_n,
        # The arm is not a detail. A ranking fitted on all lineages and one fitted
        # without the held-out lineage support different claims, and nothing
        # downstream can tell them apart unless the artefact says so.
        "arm": "fold-internal" if held_out_lineage else "refit-on-all",
        "held_out_lineage": held_out_lineage,
        "rows_used": int(len(train)),
        "rows_total": int(len(mart)),
        "features_before": len(features),
        "features_after": len(features) - len(dropped),
        "top_features": [c for c, _ in scored[:10]],
    }
    return kept, note


register(Step(
    name="rank_association",
    layer="rank",
    label_free=False,
    summary="rank features by chi-square or mutual information with the phenotype",
    run=run,
))
