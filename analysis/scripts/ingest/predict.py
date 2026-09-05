"""End-to-end predictor: ideal input → DR prediction + reasoning (design/16 Stage 2).

Closes the product loop (design/16 §Goal):

    input JSON ─▶ vectorise ─▶ concordant MOJO ─▶ P(R) ─▶ explain ─▶ prediction JSON
      (IR)        (Stage 1)      (train_h2o)              (Stage 2)

The prediction is **catalogue-free** — the H2O concordant model, not piezo,
calls R/S. The reasoning then annotates *why*, flagging genuine drug mechanisms
vs co-selected passengers (explain.py). This is the whole value proposition:
a prediction a catalogue can't make, with reasoning a catalogue can't give.

Run (pixi env: h2o + interpret):
    JAVA_HOME=$PWD/.pixi/envs/default/lib/jvm \
    .pixi/envs/default/bin/python -m analysis.scripts.ingest.predict \
        --input analysis/results/ingest/examples/mdr.input.json \
        --schema analysis/results/ingest/schemas/feature_schema_RIF_concordant.json \
        --mojo analysis/results/h2o_stacked/concordant/SE_top4_glm_71181.zip \
        --catalogue analysis/results/h2o_stacked/concordant/glm_catalogue.csv \
        --causal analysis/results/causal/RIF_mutation/causal_concordance_RIF_mutation_level.parquet
"""

from __future__ import annotations

import argparse
import json
import os
from pathlib import Path

import pandas as pd

from analysis.scripts.ingest.explain import explain, load_evidence
from analysis.scripts.ingest.vectorize import vectorize


def predict(input_path: Path, schema_path: Path, mojo_path: Path,
            catalogue_csv: Path, causal_parquet: Path) -> dict:
    instance = json.loads(input_path.read_text())
    schema = json.loads(schema_path.read_text())
    drug = schema["drug"]

    # Stage-1 vectorise
    row = vectorize(instance, schema)
    fired = [k for k, v in row.items() if k.startswith("raw__") and float(v) != 0.0]

    # Predict via the standalone MOJO (no persistent cluster needed for scoring,
    # but import_mojo is the simplest path in-process).
    os.environ.setdefault("JAVA_HOME", os.path.abspath(".pixi/envs/default/lib/jvm"))
    import h2o
    h2o.init(nthreads=2, max_mem_size="3G", name=f"predict-{os.getpid()}")
    model = h2o.import_mojo(str(mojo_path))
    frame = pd.DataFrame([{k: float(v) for k, v in row.items()}])[list(row.keys())]
    pred = model.predict(h2o.H2OFrame(frame)).as_data_frame()
    h2o.cluster().shutdown()

    p_r = float(pred["p1"].iloc[0]) if "p1" in pred.columns else float(pred.iloc[0, -1])
    call = "R" if p_r >= 0.5 else "S"

    # Stage-2 explain
    evidence = load_evidence(catalogue_csv, causal_parquet)
    reasoning = explain(fired, drug, evidence, p_r)

    return {
        "sample_id": instance["sample_id"],
        "drug": drug,
        "prediction": call,
        "probability_resistant": round(p_r, 4),
        "lineage": instance["covariates"].get("lineage"),
        "n_mutations_fired": len(fired),
        **reasoning,
        "provenance": {
            "model": Path(mojo_path).name, "feature_set": schema["feature_set"],
            "n_model_features": schema["n_features"], "catalogue_free_prediction": True,
        },
    }


def _print_human(r: dict) -> None:
    print(f"\n  {r['sample_id']}")
    print(f"  {r['drug']}: {r['prediction']}  (P(R)={r['probability_resistant']}, "
          f"lineage={r['lineage']})")
    print(f"  confidence: {r['confidence']['drug_class']} / {r['confidence']['reliability']}")
    print(f"  because ({r['n_mutations_fired']} mutation(s) fired):")
    for x in r["reasons"]:
        es = x["effect_size"]
        es_s = f"CATE={es['cate_ate']:+.2f} {es['ci']}" if es else "CATE=n/a"
        tgt = "on-target" if x["on_target"] else "OFF-target"
        print(f"    • {x['mutation']:<16} weight={x['catalogue_weight']:+.2f}  "
              f"{es_s:<26} {tgt:<10} → {x['role']}")


def main() -> None:
    p = argparse.ArgumentParser(description="File → DR prediction + reasoning (design/16 Stage 2).")
    p.add_argument("--input", type=Path, required=True)
    p.add_argument("--schema", type=Path, required=True)
    p.add_argument("--mojo", type=Path, required=True)
    p.add_argument("--catalogue", type=Path, required=True, help="glm_catalogue.csv")
    p.add_argument("--causal", type=Path, required=True, help="causal.py results parquet")
    p.add_argument("--out", type=Path, default=None)
    args = p.parse_args()

    result = predict(args.input, args.schema, args.mojo, args.catalogue, args.causal)
    _print_human(result)
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(json.dumps(result, indent=2))
        print(f"\n-> {args.out}")
    else:
        print("\n" + json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
