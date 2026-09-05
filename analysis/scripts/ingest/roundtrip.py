"""Stage-1 gate: the round-trip test (design/16).

Proves the serving path reproduces training with no skew:

    CRyPTIC DB ─▶ ideal input (IR) ─▶ vectorise ─▶ feature row
                                                        ║  must equal
    training mart row for the same isolate ═════════════╝  (byte-identical)

If the served feature row equals the row the model trained on, the MOJO score is
identical by construction (deterministic) — so byte-equality of the feature row
IS the correctness proof. `--with-mojo` additionally scores both through the MOJO
as belt-and-braces.

Because we standardised on gnomonicus (design/16 §Source tools) and the CRyPTIC
`mutations` table IS gnomonicus output, `reference_input.build_input` is the
gnomonicus-native adapter — so this round-trip is the reference every later
adapter (TB-Profiler, hAMRonization) is validated against.

Run (pixi env):
    .pixi/envs/default/bin/python -m analysis.scripts.ingest.roundtrip \
        --db analysis/databases/duckdb/cryptic-slim.duckdb \
        --mart analysis/results/feature_mart/feature_mart_RIF_cryptic-slim-2026.05_v1.0.2.parquet \
        --schema analysis/results/ingest/schemas/feature_schema_RIF_concordant.json \
        [--uniqueid <UID>] [--n 25] [--with-mojo <SE.zip>]
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd

from analysis.scripts.ingest.reference_input import build_input, build_inputs_batch, connect
from analysis.scripts.ingest.vectorize import vectorize


def _pick_isolates(mart: pd.DataFrame, schema: dict, uniqueid: str | None, n: int) -> list[str]:
    if uniqueid:
        return [uniqueid]
    ids = mart["id__UNIQUEID"].tolist()
    # Prefer some R-driving isolates (rpoB_S450L fired) so the gate exercises a
    # non-trivial feature row, then pad with arbitrary isolates.
    driver = next((f["name"] for f in schema["features"]
                   if f["kind"] == "mutation" and f["name"].endswith("rpoB_S450L")), None)
    picks: list[str] = []
    if driver and driver in mart.columns:
        picks += mart.loc[mart[driver] == 1, "id__UNIQUEID"].head(n // 2).tolist()
    for uid in ids:
        if uid not in picks:
            picks.append(uid)
        if len(picks) >= n:
            break
    return picks[:n]


def roundtrip(db_path: Path, mart_path: Path, schema_path: Path,
              uniqueid: str | None, n: int, mojo: Path | None) -> dict:
    schema = json.loads(schema_path.read_text())
    feats = [f["name"] for f in schema["features"]]
    mart = pd.read_parquet(mart_path, columns=["id__UNIQUEID", *feats])
    con = connect(db_path)

    isolates = _pick_isolates(mart, schema, uniqueid, n)
    instances = build_inputs_batch(con, isolates)  # one scan, not one-per-isolate
    mismatches: list[dict] = []
    checked = 0
    for uid in isolates:
        expected = mart.loc[mart["id__UNIQUEID"] == uid, feats]
        if expected.empty or uid not in instances:
            continue
        checked += 1
        exp = expected.iloc[0]
        got = vectorize(instances[uid], schema)
        # per-feature exact comparison, respecting dtype (float32 vs int8)
        diffs = {name: (float(exp[name]), float(got[name]))
                 for f in schema["features"]
                 for name in [f["name"]]
                 if not _equal(exp[name], got[name], f["dtype"])}
        if diffs:
            mismatches.append({"uniqueid": uid, "diffs": diffs})

    result = {
        "drug": schema["drug"], "feature_set": schema["feature_set"],
        "n_features": len(feats), "isolates_checked": checked,
        "isolates_matched": checked - len(mismatches),
        "mismatches": mismatches[:10],
        "passed": not mismatches,
    }

    if mojo and not mismatches:
        result["mojo_check"] = _mojo_agreement(mojo, mart, schema, con, isolates[: min(5, len(isolates))])
    return result


def _equal(a, b, dtype: str) -> bool:
    cast = np.float32 if dtype == "float32" else np.int8
    return cast(a) == cast(b)


def _mojo_agreement(mojo: Path, mart: pd.DataFrame, schema: dict,
                    con, isolates: list[str]) -> dict:
    """Belt-and-braces: score the mart row and the vectorised row through the
    actual MOJO; predictions must match to float precision."""
    import os
    os.environ.setdefault("JAVA_HOME", os.path.abspath(".pixi/envs/default/lib/jvm"))
    import h2o
    h2o.init(nthreads=2, max_mem_size="3G", name=f"rt-{os.getpid()}")
    model = h2o.import_mojo(str(mojo))
    feats = [f["name"] for f in schema["features"]]
    agree, maxdiff = 0, 0.0
    for uid in isolates:
        exp = mart.loc[mart["id__UNIQUEID"] == uid, feats]
        if exp.empty:
            continue
        got = vectorize(build_input(con, uid), schema)
        p_exp = _p1(model, h2o.H2OFrame(exp[feats]))
        p_got = _p1(model, h2o.H2OFrame(pd.DataFrame([{k: float(v) for k, v in got.items()}])[feats]))
        maxdiff = max(maxdiff, abs(p_exp - p_got))
        agree += int(abs(p_exp - p_got) < 1e-9)
    h2o.cluster().shutdown()
    return {"scored": len(isolates), "agree": agree, "max_abs_prob_diff": maxdiff}


def _p1(model, hf) -> float:
    pred = model.predict(hf).as_data_frame()
    return float(pred["p1"].iloc[0]) if "p1" in pred.columns else float(pred.iloc[0, -1])


def main() -> None:
    p = argparse.ArgumentParser(description="Stage-1 round-trip gate (design/16).")
    p.add_argument("--db", type=Path, default=Path("analysis/databases/duckdb/cryptic-slim.duckdb"))
    p.add_argument("--mart", type=Path, required=True)
    p.add_argument("--schema", type=Path, required=True)
    p.add_argument("--uniqueid", default=None)
    p.add_argument("--n", type=int, default=25)
    p.add_argument("--with-mojo", type=Path, default=None)
    args = p.parse_args()

    result = roundtrip(args.db, args.mart, args.schema, args.uniqueid, args.n, args.with_mojo)
    print(json.dumps(result, indent=2))
    if not result["passed"]:
        raise SystemExit(f"ROUND-TRIP FAILED: {len(result['mismatches'])} isolate(s) skewed")
    print(f"\n✓ round-trip PASSED — {result['isolates_matched']}/{result['isolates_checked']} "
          f"isolates byte-identical on {result['n_features']} features")


if __name__ == "__main__":
    main()
