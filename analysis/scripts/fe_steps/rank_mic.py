"""L4 — rank features by association with the graded MIC, inside a fold.

Pei 2024 ranks mutation-phenotype association on the median MIC and an
FDR-corrected Wilcoxon rank-sum test, not on the binarised R/S label
(doi:10.1016/S2666-5247(24)00131-9). That uses strictly more information: the
binary label thresholds the MIC at a breakpoint and throws away how far past it
an isolate sits, so two mutations that shift the MIC by very different amounts
look identical once binarised.

The statistic here is the rank-biserial effect size from a Mann-Whitney U test
between carriers and non-carriers, with Benjamini-Hochberg correction across the
candidate set. Effect size rather than p-value alone, because with fourteen
thousand isolates almost everything is significant and the ordering would be
driven by carrier count.

Label-using: the MIC IS the phenotype. Runs inside the fold, or the ranking sees
the lineage it will be tested on.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from scipy import stats

from analysis.scripts.fe_steps.registry import Step, register

LINEAGE_COL = "cov__lineage_raw"
MIC_COL = "y__log2mic"


def _benjamini_hochberg(p: np.ndarray) -> np.ndarray:
    n = len(p)
    order = np.argsort(p)
    ranked = p[order] * n / (np.arange(n) + 1)
    # enforce monotonicity from the largest p downwards
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def run(mart: pd.DataFrame, *, top_n: int = 1000, min_carriers: int = 10,
        alpha: float = 0.05, held_out_lineage: str | None = None,
        **_) -> tuple[pd.DataFrame, dict]:
    if MIC_COL not in mart.columns:
        raise ValueError(
            f"mart has no {MIC_COL}; rebuild it with a build_mart that carries the "
            "graded phenotype. The binary label cannot substitute -- ranking on it "
            "is what rank_association already does.")

    train = mart
    if held_out_lineage:
        if LINEAGE_COL not in mart.columns:
            raise ValueError(f"mart has no {LINEAGE_COL}; cannot hold a lineage out")
        train = mart.loc[mart[LINEAGE_COL] != held_out_lineage]
        if train.empty:
            raise ValueError(f"holding out {held_out_lineage} left no rows")

    mic_all = pd.to_numeric(train[MIC_COL], errors="coerce")
    usable = mic_all.notna().to_numpy()
    mic = mic_all.to_numpy(dtype=float)[usable]
    if mic.size == 0:
        raise ValueError(f"{MIC_COL} is entirely null in the training rows")

    features = [c for c in mart.columns if c.startswith("raw__")]
    names, effects, pvals = [], [], []
    for col in features:
        # numpy throughout: mixing a boolean Series with a positional mask
        # silently realigns on the index and then fails on the negation
        col_vals = pd.to_numeric(train[col], errors="coerce").fillna(0).to_numpy(float)
        carried = col_vals[usable] > 0
        n1, n0 = int(carried.sum()), int(np.count_nonzero(~carried))
        if n1 < min_carriers or n0 < min_carriers:
            continue
        u, p = stats.mannwhitneyu(mic[carried], mic[~carried], alternative="two-sided")
        # rank-biserial correlation: a bounded, sign-carrying effect size
        names.append(col)
        effects.append(abs(2 * u / (n1 * n0) - 1))
        pvals.append(p)

    if not names:
        raise ValueError("no feature had enough carriers on both sides to test")

    q = _benjamini_hochberg(np.asarray(pvals))
    # significant first, then by descending effect. Not `-(bool)` — q is a
    # numpy float, so the comparison yields np.bool_, which has no unary minus.
    ranked = sorted(zip(names, effects, q),
                    key=lambda r: (0 if r[2] <= alpha else 1, -float(r[1])))
    keep = {n for n, _, _ in ranked[:top_n]}
    dropped = [c for c in features if c not in keep]
    kept = mart.drop(columns=dropped)

    note = {
        "step": "rank_mic",
        "statistic": "Mann-Whitney U on log2 MIC, rank-biserial effect, BH-corrected",
        "top_n": top_n,
        "alpha": alpha,
        "min_carriers": min_carriers,
        "arm": "fold-internal" if held_out_lineage else "refit-on-all",
        "held_out_lineage": held_out_lineage,
        "rows_used": int(np.count_nonzero(usable)),
        "rows_total": int(len(mart)),
        "features_before": len(features),
        "features_tested": len(names),
        "features_after": len(features) - len(dropped),
        "significant_after_fdr": int(sum(1 for _, _, qq in ranked if qq <= alpha)),
        "top_features": [n for n, _, _ in ranked[:10]],
    }
    return kept, note


register(Step(
    name="rank_mic",
    layer="rank",
    label_free=False,
    summary="rank features by Mann-Whitney effect on the graded MIC (Pei 2024)",
    run=run,
))
