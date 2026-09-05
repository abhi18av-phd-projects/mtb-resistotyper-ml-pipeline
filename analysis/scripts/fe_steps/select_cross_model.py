"""L5 — keep features two different model classes independently agree on.

Ghosh 2026 selects by taking the top-N features from a random forest and from a
multilayer perceptron, ranked by SHAP importance, then intersecting them
(doi:10.26508/lsa.202503539). On their data the union was 62 features and the
intersection 7 -- a wide precision/recall dial from one parameter.

This is a different KIND of voter from the concordance layer. Concordance asks
three causal questions of one candidate (glass-box importance, treatment effect,
lineage-residual signal). Cross-model agreement asks one question of two
inductive biases: does a feature matter to a model that splits, and to a model
that projects? A feature both rely on is unlikely to be an artefact of either
one's geometry.

Cheap, too: two fits and a set operation, no causal machinery.

Label-using. Runs inside the fold, or both models see the lineage they will be
tested on and the agreement is worthless.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from analysis.scripts.fe_steps.registry import Step, register

LINEAGE_COL = "cov__lineage_raw"
LABEL_COL = "y__binary"
MODES = ("intersection", "union")


def _importances(X: np.ndarray, y: np.ndarray, seed: int) -> dict[str, np.ndarray]:
    from sklearn.ensemble import RandomForestClassifier
    from sklearn.neural_network import MLPClassifier
    from sklearn.preprocessing import StandardScaler

    rf = RandomForestClassifier(n_estimators=300, n_jobs=-1, random_state=seed)
    rf.fit(X, y)

    # The MLP has no native importance, and SHAP on a net is slow enough to
    # dominate the step. Permutation importance on the fitted net answers the
    # same question -- how much does the model lean on this feature -- without
    # the extra dependency.
    from sklearn.inspection import permutation_importance
    Xs = StandardScaler().fit_transform(X)
    mlp = MLPClassifier(hidden_layer_sizes=(64, 32), alpha=1e-2, max_iter=300,
                        random_state=seed)
    mlp.fit(Xs, y)
    perm = permutation_importance(mlp, Xs, y, n_repeats=3, random_state=seed, n_jobs=-1)
    return {"rf": rf.feature_importances_, "mlp": perm.importances_mean}


def run(mart: pd.DataFrame, *, top_n: int = 50, mode: str = "intersection",
        seed: int = 42, held_out_lineage: str | None = None,
        **_) -> tuple[pd.DataFrame, dict]:
    if mode not in MODES:
        raise ValueError(f"unknown mode {mode!r}; known: {MODES}")
    if LABEL_COL not in mart.columns:
        raise ValueError(f"mart has no {LABEL_COL}")

    train = mart
    if held_out_lineage:
        if LINEAGE_COL not in mart.columns:
            raise ValueError(f"mart has no {LINEAGE_COL}; cannot hold a lineage out")
        train = mart.loc[mart[LINEAGE_COL] != held_out_lineage]
        if train.empty:
            raise ValueError(f"holding out {held_out_lineage} left no rows")

    features = [c for c in mart.columns if c.startswith("raw__")]
    X = train[features].apply(pd.to_numeric, errors="coerce").fillna(0).to_numpy(float)
    y = pd.to_numeric(train[LABEL_COL], errors="coerce").fillna(0).to_numpy(int)
    if len(np.unique(y)) < 2:
        raise ValueError("training rows carry only one class; cannot fit either model")

    imp = _importances(X, y, seed)
    picks = {k: {features[i] for i in np.argsort(v)[::-1][:top_n]} for k, v in imp.items()}
    inter, union = picks["rf"] & picks["mlp"], picks["rf"] | picks["mlp"]
    keep = inter if mode == "intersection" else union
    if not keep:
        raise ValueError(f"{mode} of the two top-{top_n} sets is empty")

    dropped = [c for c in features if c not in keep]
    kept = mart.drop(columns=dropped)

    note = {
        "step": "select_cross_model",
        "models": ["random_forest", "mlp_permutation"],
        "mode": mode,
        "top_n": top_n,
        "arm": "fold-internal" if held_out_lineage else "refit-on-all",
        "held_out_lineage": held_out_lineage,
        "rows_used": int(len(train)),
        "rows_total": int(len(mart)),
        "features_before": len(features),
        "n_rf": len(picks["rf"]),
        "n_mlp": len(picks["mlp"]),
        "n_intersection": len(inter),
        "n_union": len(union),
        "agreement": round(len(inter) / len(union), 4) if union else 0.0,
        "features_after": len(features) - len(dropped),
        "kept": sorted(keep)[:15],
    }
    return kept, note


register(Step(
    name="select_cross_model",
    layer="select",
    label_free=False,
    summary="keep features a random forest and an MLP independently rank top-N (Ghosh 2026)",
    run=run,
))
