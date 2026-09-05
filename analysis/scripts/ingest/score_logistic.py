"""Pure-numpy serving for the deployed L2-logistic per-drug models (Step 1).

Replaces the H2O MOJO scoring path for the usable/moderate tiers. Because the model
is StandardScaler + L2-logistic, scoring is a few lines of numpy — no H2O, no JVM, no
persistent cluster — and the per-feature **logit contributions** are an *exact*
attribution (coef_i · z_i), so prediction and explanation come from the same closed
form (no surrogate model, unlike the GLM-over-SE approach).

    p(R) = sigmoid( sum_i coef_i * (x_i - mean_i)/scale_i + intercept )

Model JSON is produced by `feature_mart/train_logistic.py`.

Usage:
    from analysis.scripts.ingest.score_logistic import load_model, score
    model = load_model("analysis/results/fulldata_ml/logistic_models/RIF_logistic.json")
    result = score(model, {"raw__rpoB_S450L": 1, "cov__median_coverage": 55.0, ...})
    # -> {"p_resistant": 0.97, "prediction": "R", "contributions": [...]}
"""

from __future__ import annotations

import json
import math
from pathlib import Path


def load_model(path: str | Path) -> dict:
    return json.loads(Path(path).read_text())


def score(model: dict, features: dict, threshold: float = 0.5) -> dict:
    """Score one isolate. `features` maps feature name -> value. Absent features
    default sensibly: a missing `raw__`/mutation feature is 0 (variant not carried),
    but a missing `cov__` covariate is mean-imputed (z=0, neutral) so a genotype-only
    call is not swamped by an out-of-distribution or unknown quality covariate."""
    names = model["features"]
    mean = model["scaler_mean"]
    scale = model["scaler_scale"]
    coef = model["coef"]
    logit = model["intercept"]
    contribs = []
    for i, name in enumerate(names):
        if name in features:
            x = float(features[name])
        elif name.startswith("cov__"):
            x = mean[i]           # neutral covariate (z = 0) when not supplied
        else:
            x = 0.0               # mutation not carried
        z = (x - mean[i]) / scale[i] if scale[i] else 0.0
        c = coef[i] * z
        logit += c
        if x != 0 or abs(coef[i]) > 0:  # report carried variants + all determinants
            contribs.append({"feature": name, "value": x, "coef": round(coef[i], 4),
                             "logit_contribution": round(c, 4),
                             "direction": "R+" if coef[i] > 0 else "S-"})
    p = 1.0 / (1.0 + math.exp(-logit))
    contribs.sort(key=lambda r: -abs(r["logit_contribution"]))
    return {"drug": model["drug"], "p_resistant": round(p, 4),
            "prediction": "R" if p >= threshold else "S",
            "tuned_C": model.get("C"),
            "contributions": contribs}


if __name__ == "__main__":  # tiny smoke test
    import argparse
    ap = argparse.ArgumentParser(description="Score an isolate with a deployed logistic model.")
    ap.add_argument("--model", type=Path, required=True)
    ap.add_argument("--features", type=Path, required=True, help="JSON dict feature->value")
    a = ap.parse_args()
    print(json.dumps(score(load_model(a.model), json.loads(a.features.read_text())), indent=2))
