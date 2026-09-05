"""Explainer: turn a prediction into per-mutation reasoning (design/16 Stage 2).

For each mutation that FIRED in the concordant model, assemble the evidence we
already produce at training time — no new modelling:

  catalogue_weight  the GLM signed coefficient (a learned-catalogue weight;
                    glm_catalogue.csv). + = resistance-associated, − = protective.
  effect_size       the causal CATE (average treatment effect + 95% CI) from
                    causal.py — is the effect large AND significant, or noise?
  ebm_importance    how much the glass-box model actually relies on it.
  on_target         is this gene a known determinant of THIS drug? (rpoB→RIF).
  role              derived, honest label:
       primary_driver              on-target + significant CATE (the real mechanism)
       co-selected passenger       off-target + significant CATE (MDR linkage —
                                   e.g. katG S315T riding along in a RIF model)
       weak/lineage-linked passenger  CATE CI crosses 0 / negligible EBM (noise
                                   the concordance filter let through)

This is the project's differentiator: the model predicts catalogue-free, and the
explanation says *why* AND flags which "reasons" are genuine drug mechanisms vs
co-travellers — something a catalogue lookup cannot do.

Consumed by predict.py; not usually run standalone.
"""

from __future__ import annotations

import pandas as pd

# Canonical TB resistance determinants (gene → the drugs it confers resistance
# to). Used only to label on- vs off-target; the PREDICTION never uses it.
DETERMINANTS: dict[str, set[str]] = {
    "RIF": {"rpoB"},
    "INH": {"katG", "inhA", "fabG1", "ahpC"},
    "EMB": {"embB", "embA", "embC", "embR"},
    "PZA": {"pncA", "rpsA", "panD"},
    "MXF": {"gyrA", "gyrB"}, "LEV": {"gyrA", "gyrB"},
    "AMI": {"rrs", "eis"}, "KAN": {"rrs", "eis"}, "STM": {"rpsL", "rrs", "gid"},
    "ETH": {"ethA", "ethR", "inhA", "fabG1"},
    "LZD": {"rrl", "rplC"},
    "BDQ": {"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"},
    "CFZ": {"Rv0678", "pepQ", "mmpL5", "mmpS5"},
    "DLM": {"ddn", "fbiA", "fbiB", "fbiC", "fgd1", "Rv0678"},
}

# Honest AUC tiers (design/03 §honest evaluation, design/15). Refresh from
# evaluate_cv.py's lineage-LOO once the causal venv is rebuilt.
CONFIDENCE_TIER = {
    "RIF": ("first-line", "usable"), "INH": ("first-line", "usable"),
    "EMB": ("first-line", "moderate"), "PZA": ("first-line", "moderate"),
    "MXF": ("second-line", "weak"), "LEV": ("second-line", "weak"),
    "AMI": ("second-line", "weak"), "KAN": ("second-line", "weak"),
    "ETH": ("second-line", "weak"),
    "BDQ": ("new-drug", "at-chance"), "DLM": ("new-drug", "at-chance"),
    "LZD": ("new-drug", "at-chance"), "CFZ": ("new-drug", "at-chance"),
}


def _feature_gene_garc(feature: str) -> tuple[str, str]:
    """raw__<gene>_<garc> → (gene, garc). GARC may contain '_', so split once."""
    body = feature[len("raw__"):]
    gene, _, garc = body.partition("_")
    return gene, garc


def load_evidence(catalogue_csv, causal_parquet) -> dict:
    """Index the GLM catalogue coefficients + causal votes by feature name."""
    cat = pd.read_csv(catalogue_csv)
    coef = dict(zip(cat["feature"], cat["coefficient"]))
    caus = pd.read_parquet(causal_parquet).set_index("feature_column")
    return {"coef": coef, "causal": caus}


def _role(on_target: bool, cate_sig: bool, cate_ate: float, ebm_imp: float) -> str:
    if on_target and cate_sig:
        return "primary_driver"
    if cate_sig and abs(cate_ate) >= 0.2:
        return "co-selected passenger (MDR linkage)"
    return "weak/lineage-linked passenger"


def explain_mutation(feature: str, drug: str, evidence: dict) -> dict:
    gene, garc = _feature_gene_garc(feature)
    coef = float(evidence["coef"].get(feature, 0.0))
    on_target = gene in DETERMINANTS.get(drug, set())
    c = evidence["causal"]
    cate_sig = cate_ate = ebm_imp = None
    lo = hi = None
    votes = None
    if feature in c.index:
        row = c.loc[feature]
        cate_sig = bool(row["cate_significant"])
        cate_ate = float(row["cate_ate"])
        lo, hi = float(row["cate_ci_lower"]), float(row["cate_ci_upper"])
        ebm_imp = float(row["ebm_importance"])
        votes = int(row["n_votes"])
    return {
        "mutation": f"{gene}@{garc}",
        "gene": gene,
        "catalogue_weight": round(coef, 3),
        "direction": "resistance" if coef > 0 else "protective",
        "on_target": on_target,
        "effect_size": None if cate_ate is None else {
            "cate_ate": round(cate_ate, 3), "ci": [round(lo, 3), round(hi, 3)],
            "significant": cate_sig,
        },
        "ebm_importance": None if ebm_imp is None else round(ebm_imp, 4),
        "concordance_votes": votes,
        "role": _role(on_target, bool(cate_sig), cate_ate or 0.0, ebm_imp or 0.0),
    }


def explain(fired_features: list[str], drug: str, evidence: dict, p_r: float) -> dict:
    reasons = [explain_mutation(f, drug, evidence) for f in fired_features]
    # rank: primary drivers first, then by |catalogue weight|
    order = {"primary_driver": 0, "co-selected passenger (MDR linkage)": 1,
             "weak/lineage-linked passenger": 2}
    reasons.sort(key=lambda r: (order.get(r["role"], 3), -abs(r["catalogue_weight"])))
    drivers = [r["mutation"] for r in reasons if r["role"] == "primary_driver"]
    passengers = [r["mutation"] for r in reasons if "passenger" in r["role"]]
    line, verdict = CONFIDENCE_TIER.get(drug, ("unknown", "unknown"))
    return {
        "primary_drivers": drivers,
        "passengers_flagged": passengers,
        "confidence": {"drug_class": line, "reliability": verdict,
                       "caveat": f"{drug} is {line}; model reliability: {verdict} "
                                 f"(honest lineage-aware CV tier, design/03)."},
        "reasons": reasons,
    }
