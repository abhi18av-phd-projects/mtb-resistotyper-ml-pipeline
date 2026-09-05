#!/usr/bin/env python3
"""Derive a reduced database from the full one, and say what was lost.

The reduction is a convenience for fast iteration, not a scientific artefact.
It drops the plate-level tables and the two wide `effects` columns, which is
most of the size and also exactly what three feature-engineering techniques
need. The manifest records the loss so nothing downstream has to infer it.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb

# Kept: everything the per-drug cohort and mart construction actually read.
KEEP_TABLES = [
    "_database_metadata", "countries", "drugs", "sites",
    "dst_samples", "wgs_samples", "genomes",
    "dst_measurements", "ukmyc_phenotypes",
    "effects", "predictions", "mutations",
]
# Dropped from `effects`: the two wide JSON/string columns.
DROP_EFFECTS_COLS = ["EVIDENCE", "PREDICTION_VALUES"]

CAPABILITY_LOSS = {
    "EVIDENCE": ["T3 EVIDENCE JSON — Observed_samples, FINAL CONFIDENCE GRADING, WHO_HGVS",
                 "T4 real WHO tier grading (cat__n_tier1/tier2/interim)"],
    "PREDICTION_VALUES": ["T2 RFUS / SOLO grading"],
}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--full", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--manifest", type=Path, required=True)
    a = ap.parse_args()

    if a.out.exists():
        a.out.unlink()
    con = duckdb.connect(str(a.out))
    con.execute(f"ATTACH '{a.full}' AS full_db (READ_ONLY)")

    present = {r[0] for r in con.execute(
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_catalog='full_db' AND table_schema='main'").fetchall()}

    copied, dropped_tables, dropped_cols = [], [], []
    for t in KEEP_TABLES:
        if t not in present:
            continue
        cols = [r[0] for r in con.execute(
            "SELECT column_name FROM information_schema.columns "
            f"WHERE table_catalog='full_db' AND table_name='{t}'").fetchall()]
        if t == "effects":
            keep = [c for c in cols if c not in DROP_EFFECTS_COLS]
            dropped_cols = [c for c in cols if c in DROP_EFFECTS_COLS]
            sel = ", ".join(f'"{c}"' for c in keep)
        else:
            sel = "*"
        con.execute(f'CREATE TABLE main."{t}" AS SELECT {sel} FROM full_db.main."{t}"')
        copied.append(t)
    dropped_tables = sorted(present - set(copied) - {"source_provenance"})

    con.execute("CREATE OR REPLACE TABLE slim_provenance (key VARCHAR, value VARCHAR)")
    con.executemany("INSERT INTO slim_provenance VALUES (?, ?)", [
        ("derived_from", str(a.full)),
        ("dropped_tables", json.dumps(dropped_tables)),
        ("dropped_effects_columns", json.dumps(dropped_cols)),
        ("blocked_techniques", json.dumps(
            sorted({t for c in dropped_cols for t in CAPABILITY_LOSS.get(c, [])}))),
    ])
    con.execute("DETACH full_db")
    con.execute("CHECKPOINT")
    con.close()

    manifest = {
        "derived_from": str(a.full),
        "tables_kept": copied,
        "tables_dropped": dropped_tables,
        "effects_columns_dropped": dropped_cols,
        "blocked_techniques": sorted(
            {t for c in dropped_cols for t in CAPABILITY_LOSS.get(c, [])}),
        "note": ("A reduced copy for fast iteration. Feature engineering that needs the "
                 "dropped columns must run against the full database."),
    }
    a.manifest.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"slim: {len(copied)} tables kept, {len(dropped_tables)} dropped, "
          f"effects columns dropped: {dropped_cols}")


if __name__ == "__main__":
    main()
