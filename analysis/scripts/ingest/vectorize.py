"""Vectoriser: ideal input (IR) -> the exact model feature row (design/16 Stage 1).

The bridge between the ingest contract (`input_spec.schema.json`) and the model
contract (`feature_schema.json`). It takes an ideal-input instance (GARC variants
+ covariates) and produces the one-row feature vector the concordant MOJO scores,
in the schema's order and dtypes. Getting this byte-identical to the training
mart row is the whole point — see `roundtrip.py` for the gate.

Construction rules (mirroring build_mart.py exactly):
  mutation feature raw__<gene>_<mutation>  = 1 iff the IR carries a variant whose
       f"raw__{gene}_{mutation}" equals the feature name, else 0. Key is built by
       VERBATIM concat (build_mart.py:172 `'raw__'||GENE||'_'||MUTATION`) — no
       sanitisation, so `X`-wildcard calls (rpoB S450X) simply don't match
       raw__rpoB_S450L and become 0 (documented Stage-0 default: failed call
       treated as wild-type for the concordant model).
  lineage  cov__lineage_Lk = 1 iff covariates.lineage == the k it tests.
  coverage cov__median_coverage/tb_breadth = the covariate value, CAST to the
       schema dtype (float32) so it matches the mart's DuckDB ::FLOAT storage;
       missing -> the schema's training-median fill.

Run:
    .pixi/envs/default/bin/python -m analysis.scripts.ingest.vectorize \
        --input analysis/results/ingest/examples/mdr.input.json \
        --schema analysis/results/ingest/schemas/feature_schema_RIF_concordant.json
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

_CAST = {"int8": np.int8, "float32": np.float32}


def vectorize(instance: dict, schema: dict) -> dict:
    """IR instance + feature schema -> ordered {feature_name: value} dict."""
    present = {f"raw__{v['gene']}_{v['mutation']}" for v in instance.get("variants", [])}
    covs = instance.get("covariates", {})
    lineage = covs.get("lineage")

    row: dict = {}
    for f in schema["features"]:
        kind, name = f["kind"], f["name"]
        if kind == "mutation":
            val = 1 if name in present else 0
        elif kind == "lineage":
            val = 1 if lineage == f["lineage"] else 0
        elif kind == "coverage":
            raw = covs.get(f["source"])
            val = f["fill"] if raw is None else raw
        else:
            raise ValueError(f"unknown feature kind {kind!r}")
        row[name] = _CAST[f["dtype"]](val)
    return row


def main() -> None:
    p = argparse.ArgumentParser(description="Vectorise an ideal-input JSON to the model feature row (design/16 Stage 1).")
    p.add_argument("--input", type=Path, required=True, help="ideal-input JSON (IR instance).")
    p.add_argument("--schema", type=Path, required=True, help="feature_schema.json.")
    p.add_argument("--out", type=Path, default=None, help="write the row as JSON; else print a summary.")
    args = p.parse_args()

    instance = json.loads(args.input.read_text())
    schema = json.loads(args.schema.read_text())
    row = vectorize(instance, schema)

    fired = [k for k, v in row.items() if k.startswith("raw__") and float(v) != 0.0]
    print(f"{instance['sample_id']}: {len(row)} features, {len(fired)} mutation(s) fired: {fired}")
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps({k: float(v) for k, v in row.items()}, indent=2))
        print(f"-> {args.out}")


if __name__ == "__main__":
    main()
