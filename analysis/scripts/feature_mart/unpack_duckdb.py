#!/usr/bin/env python3
"""Unpack a published DuckDB back into the files the pipeline already expects.

The inverse of create_duckdb, and the reason that artefact is worth depositing:
somebody who downloads it should be able to continue from where it was cut,
without the compendium, the cluster, or the twelve hours of feature engineering
and causal selection that produced it.

Deliberately dumb. It writes exactly the four files each training unit is fed
elsewhere in this pipeline:

    feature_mart_<DRUG>.parquet           the wide matrix
    feature_mart_<DRUG>.metadata.json     its provenance
    concordance.parquet                   the causal verdicts for one fold
    selection.json                        that fold's manifest

so TRAIN_H2O receives the same tuple whether the selection happened here an hour
ago or in a deposit a year old. Nothing downstream needs to know which.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--db", required=True, help="published DuckDB")
    ap.add_argument("--out", required=True, help="directory to write into")
    ap.add_argument("--drugs", default="", help="comma-separated subset")
    ap.add_argument("--arm", default="", help="restrict to one campaign arm")
    ap.add_argument("--index", default="units.json",
                    help="written to --out; lists the training units found")
    a = ap.parse_args()

    out = Path(a.out); out.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(a.db, read_only=True)

    tabs = {r[0] for r in con.execute(
        "SELECT table_name FROM duckdb_tables()").fetchall()}
    for need in ("marts_index", "selections"):
        if need not in tabs:
            raise SystemExit(f"{a.db} has no '{need}' table; it was not written by "
                             f"create_duckdb, or it is truncated.")

    want = [d.strip() for d in a.drugs.split(",") if d.strip()]
    q = "SELECT table_name, drug, arm, cryptic_version FROM marts_index WHERE 1=1"
    if want: q += " AND drug IN (" + ",".join(f"'{d}'" for d in want) + ")"
    if a.arm: q += f" AND arm = '{a.arm}'"
    marts = con.execute(q).fetchall()
    if not marts:
        raise SystemExit(f"no marts in {a.db} matching drugs={want or 'any'} arm={a.arm or 'any'}")

    has_man = "selection_manifests" in tabs
    units = []
    for tbl, drug, arm, rel in marts:
        d = out / f"{drug}"
        d.mkdir(parents=True, exist_ok=True)
        mart_p = d / f"feature_mart_{drug}.parquet"
        con.execute(f"COPY (SELECT * FROM \"{tbl}\") TO '{mart_p}' (FORMAT PARQUET)")

        prov = {k: v for k, v in con.execute(
            "SELECT key, value FROM provenance WHERE drug = ? AND arm = ?",
            [drug, arm]).fetchall()} if "provenance" in tabs else {}
        prov.setdefault("cryptic_version", rel)
        prov["unpacked_from"] = str(Path(a.db).name)
        (d / f"feature_mart_{drug}.metadata.json").write_text(
            json.dumps(prov, indent=2) + "\n")

        folds = [r[0] for r in con.execute(
            "SELECT DISTINCT held_out FROM selections WHERE drug = ? AND arm = ? "
            "AND held_out IS NOT NULL ORDER BY 1", [drug, arm]).fetchall()]
        if not folds:
            print(f"  {drug}/{arm}: mart written, but no selections — "
                  f"training can only run on the full feature set")
        for fold in folds:
            fd = d / f"fold-{fold}"
            fd.mkdir(parents=True, exist_ok=True)
            # Drop the keys we added on the way in, so the parquet is byte-shaped
            # like the one the selection process originally emitted.
            cols = [c[1] for c in con.execute("PRAGMA table_info(selections)").fetchall()]
            keep = [c for c in cols if c not in ("drug", "arm", "held_out", "cryptic_version")]
            sel_cols = ", ".join(f'"{c}"' for c in keep)
            con.execute(
                f"COPY (SELECT {sel_cols} FROM selections WHERE drug = ? AND arm = ? "
                f"AND held_out = ?) TO '{fd/'concordance.parquet'}' (FORMAT PARQUET)",
                [drug, arm, fold])
            man = {}
            if has_man:
                row = con.execute(
                    "SELECT manifest FROM selection_manifests WHERE drug = ? AND arm = ? "
                    "AND held_out = ? LIMIT 1", [drug, arm, fold]).fetchone()
                if row: man = json.loads(row[0])
            man.setdefault("drug", drug)
            man.setdefault("held_out_lineage", fold)
            man["unpacked_from"] = str(Path(a.db).name)
            (fd / "selection.json").write_text(json.dumps(man, indent=2) + "\n")
            units.append({"drug": drug, "arm": arm, "held_out": fold,
                          "cryptic_version": rel,
                          "mart": str(mart_p),
                          "meta": str(d / f"feature_mart_{drug}.metadata.json"),
                          "concordance": str(fd / "concordance.parquet"),
                          "selection": str(fd / "selection.json")})
        print(f"  {drug}/{arm} [cryptic {rel}]: mart + {len(folds)} fold(s)")

    (out / a.index).write_text(json.dumps(units, indent=2) + "\n")
    print(f"\n{len(units)} training unit(s) -> {out/a.index}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
