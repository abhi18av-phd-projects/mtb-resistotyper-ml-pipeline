"""Emit `feature_schema.json` — the training/serving contract (design/16 Stage 1).

The one artifact that prevents training/serving skew: it pins the EXACT feature
list the model consumes, in training order, with each feature's kind, dtype, and
fill value. The vectoriser (`vectorize.py`) is enforced against it; no adapter is
trusted until its output vectorises cleanly through this schema.

Scope: the **concordant** model (the deployment target, design/16) — its feature
set is the causal-concordant mutations (`qf__causal_concordant` from causal.py)
plus the lineage/coverage covariates, matching `train_h2o_stacked.py`'s
`concordant` feature columns exactly.

Feature kinds:
  mutation  raw__<gene>_<mutation> one-hot; fill 0 (absent variant = wild-type)
  lineage   cov__lineage_Lk one-hot; fill 0; carries the lineage string it tests
  coverage  cov__median_coverage / cov__tb_breadth; float32; fill = training
            median (used when a covariate is missing at serve time)

Run:
    .pixi/envs/default/bin/python -m analysis.scripts.ingest.feature_schema \
        --mart analysis/results/feature_mart/feature_mart_RIF_cryptic-slim-2026.05_v1.0.2.parquet \
        --concordance analysis/results/causal/RIF_mutation/causal_concordance_RIF_mutation_level.parquet \
        --out analysis/results/ingest/schemas
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

SCHEMA_VERSION = "1.0.0"
# training-time covariate order (must match train_h2o_stacked._covariates)
COVARIATE_ORDER = [
    "cov__lineage_L1", "cov__lineage_L2", "cov__lineage_L3", "cov__lineage_L4",
    "cov__median_coverage", "cov__tb_breadth",
]
LINEAGE_MAP = {
    "cov__lineage_L1": "lineage1", "cov__lineage_L2": "lineage2",
    "cov__lineage_L3": "lineage3", "cov__lineage_L4": "lineage4",
}
COVERAGE_SOURCE = {"cov__median_coverage": "median_coverage", "cov__tb_breadth": "breadth"}


def build_schema(mart_path: Path, concordance_path: Path, drug: str) -> dict:
    mart = pd.read_parquet(mart_path)
    res = pd.read_parquet(concordance_path)
    concordant = list(res.loc[res["qf__causal_concordant"], "feature_column"])

    missing = [c for c in concordant if c not in mart.columns]
    if missing:
        raise SystemExit(f"concordant features absent from mart: {missing}")

    features: list[dict] = []
    # 1) concordant mutations, in the concordance-parquet order (training order)
    for name in concordant:
        features.append({"name": name, "kind": "mutation", "dtype": "int8", "fill": 0})
    # 2) covariates, in the fixed training order
    for cov in COVARIATE_ORDER:
        if cov not in mart.columns:
            continue
        if cov.startswith("cov__lineage_L"):
            features.append({"name": cov, "kind": "lineage", "lineage": LINEAGE_MAP[cov],
                             "dtype": "int8", "fill": 0})
        else:
            features.append({"name": cov, "kind": "coverage", "source": COVERAGE_SOURCE[cov],
                             "dtype": "float32", "fill": round(float(mart[cov].median()), 6)})

    return {
        "schema_version": SCHEMA_VERSION,
        "drug": drug,
        "feature_set": "concordant",
        "n_features": len(features),
        "raw_key_rule": "raw__{gene}_{mutation}",
        "note": "mutation fill 0 = absent variant treated as wild-type; coverage fill = training median.",
        "provenance": {
            "mart": Path(mart_path).name,
            "concordance": Path(concordance_path).name,
        },
        "features": features,
    }


def _infer_drug(mart_path: Path) -> str:
    name = mart_path.name
    return name[len("feature_mart_"):].split("_")[0] if name.startswith("feature_mart_") else "UNKNOWN"


def main() -> None:
    p = argparse.ArgumentParser(description="Emit feature_schema.json for the concordant model (design/16 Stage 1).")
    p.add_argument("--mart", type=Path, required=True)
    p.add_argument("--concordance", type=Path, required=True)
    p.add_argument("--drug", default=None)
    p.add_argument("--out", type=Path, default=Path("analysis/results/ingest/schemas"))
    args = p.parse_args()

    drug = args.drug or _infer_drug(args.mart)
    schema = build_schema(args.mart, args.concordance, drug)
    args.out.mkdir(parents=True, exist_ok=True)
    path = args.out / f"feature_schema_{drug}_concordant.json"
    path.write_text(json.dumps(schema, indent=2))
    print(f"{drug}: {schema['n_features']} features "
          f"({sum(f['kind']=='mutation' for f in schema['features'])} mutations "
          f"+ {sum(f['kind']!='mutation' for f in schema['features'])} covariates) -> {path}")


if __name__ == "__main__":
    main()
