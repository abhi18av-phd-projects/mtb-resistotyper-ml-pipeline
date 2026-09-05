#!/usr/bin/env python3
"""Apply a chain of feature-engineering steps to a mart.

One uniform contract for every technique -- mart in, mart out, plus a note saying
what it did -- so adding a technique never means touching the pipeline.

The phase argument is checked against each step's label_free property before any
work happens, so a label-using step placed before the mart is a configuration
error rather than a mart that looks fine and is not.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from analysis.scripts.fe_steps import registry


def parse_steps(spec: str) -> list[tuple[str, dict]]:
    """'region_filter,rank_association:method=mi;top_n=500' -> [(name, params)]"""
    out: list[tuple[str, dict]] = []
    for chunk in (s.strip() for s in spec.split(",")):
        if not chunk:
            continue
        name, _, argstr = chunk.partition(":")
        params: dict[str, object] = {}
        for kv in (a for a in argstr.split(";") if a):
            k, _, v = kv.partition("=")
            params[k.strip()] = _coerce(v.strip())
        out.append((name.strip(), params))
    return out


def _coerce(v: str) -> object:
    if v.lower() in ("true", "false"):
        return v.lower() == "true"
    for cast in (int, float):
        try:
            return cast(v)
        except ValueError:
            pass
    return v


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--mart", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--steps", required=True,
                    help="comma-separated, each optionally name:k=v;k=v")
    ap.add_argument("--phase", choices=["pre-mart", "per-fold"], required=True)
    ap.add_argument("--annotations", type=Path, default=None)
    ap.add_argument("--held-out-lineage", default=None)
    ap.add_argument("--notes", type=Path, default=None)
    ap.add_argument("--list", action="store_true", help="print the registry and exit")
    a = ap.parse_args()

    registry.load_builtin_steps()

    if a.list:
        for s in registry.all_steps():
            print(f"{s.name:20} {s.layer:9} {s.phase:9} {s.summary}")
        return

    chain = parse_steps(a.steps)
    registry.assert_placement([n for n, _ in chain], phase=a.phase)

    mart = pd.read_parquet(a.mart)
    annotations = pd.read_parquet(a.annotations) if a.annotations else None

    notes = []
    for name, params in chain:
        step = registry.get(name)
        if not step.label_free and a.held_out_lineage:
            params.setdefault("held_out_lineage", a.held_out_lineage)
        before = mart.shape
        mart, note = step.run(mart, annotations=annotations, **params)
        note |= {"layer": step.layer, "label_free": step.label_free,
                 "shape_before": list(before), "shape_after": list(mart.shape)}
        notes.append(note)
        print(f"{name}: {before[1]} -> {mart.shape[1]} columns")

    a.out.parent.mkdir(parents=True, exist_ok=True)
    mart.to_parquet(a.out, index=False)

    summary = {
        "phase": a.phase,
        "held_out_lineage": a.held_out_lineage,
        # True only if every step in this chain was label-free. A mart that fails
        # this is fine for a deployment refit and invalid for evaluation.
        "label_safe": all(n["label_free"] for n in notes),
        "steps": notes,
    }
    (a.notes or a.out.with_suffix(".steps.json")).write_text(
        json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
