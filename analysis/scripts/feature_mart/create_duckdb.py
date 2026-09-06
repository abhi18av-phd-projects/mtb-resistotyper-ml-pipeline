#!/usr/bin/env python3
"""Assemble one publishable DuckDB from a campaign's FE and selection outputs.

This file is not read by the pipeline. It exists so the state of the data
BETWEEN feature engineering and model training can be deposited somewhere
citable, as a single artefact a reader can open without the pipeline, the
cluster, or the compendium it came from.

Self-contained on purpose. A deposit that references parquet paths outside
itself is not a deposit, so the marts are copied in rather than pointed at.
The feature matrix measured 56 % dense, so it is stored wide; a sparse
long-format triple would have been larger than the thing it replaced.

Shape:
    mart_<DRUG>_<ARM>   one wide table per drug and arm — the matrix that
                        trained that model: id__, y__, cov__, qf__, cat__,
                        raw__ and __fold__ columns, one row per isolate
    selections          every candidate feature, per drug, arm and fold, with
                        the three-way causal evidence and the concordance verdict
    marts_index         what tables exist, their grain and their size
    provenance          compendium release, catalogue, mart version, and the
                        leakage boundary each mart was built under
"""
from __future__ import annotations

import argparse
import json
import re
from pathlib import Path

import duckdb

SAFE = re.compile(r"[^A-Za-z0-9_]")


def _ident(*parts: str) -> str:
    return SAFE.sub("_", "_".join(p for p in parts if p)).strip("_")


def _q(v) -> str:
    """SQL literal, or NULL."""
    return "NULL" if v is None else "'" + str(v).replace("'", "''") + "'"


def _read_meta(p: Path) -> dict:
    try:
        return json.loads(p.read_text())
    except (OSError, json.JSONDecodeError):
        return {}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--marts", nargs="+", required=True,
                    help="feature_mart_<DRUG>.parquet files")
    ap.add_argument("--selections", nargs="*", default=[],
                    help="concordance parquet files")
    ap.add_argument("--manifests", nargs="*", default=[],
                    help="selection.json manifests, paired with --selections by "
                         "directory. They carry the drug and the held-out fold, "
                         "which the concordance parquet does not.")
    ap.add_argument("--out", required=True, help="DuckDB file to write")
    ap.add_argument("--arm", default="", help="campaign arm this run belongs to")
    a = ap.parse_args()

    out = Path(a.out)
    if out.exists():
        out.unlink()
    out.parent.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(out))
    con.execute("SET memory_limit='8GB'")

    # cryptic_version is a first-class column, not a provenance key/value. It is
    # the axis a reader filters on before anything else — two marts for the same
    # drug differ chiefly by which compendium release built them — and making it
    # a column means that question needs no join and no string parsing.
    con.execute("""CREATE TABLE marts_index(
        table_name VARCHAR, drug VARCHAR, arm VARCHAR, cryptic_version VARCHAR,
        mart_version VARCHAR, catalogue VARCHAR,
        n_rows BIGINT, n_columns BIGINT, n_features BIGINT,
        n_resistant BIGINT, n_susceptible BIGINT)""")
    con.execute("""CREATE TABLE provenance(
        drug VARCHAR, arm VARCHAR, key VARCHAR, value VARCHAR)""")

    for m in a.marts:
        mp = Path(m)
        if not mp.exists() or mp.stat().st_size == 0:
            print(f"  skip (missing): {mp}")
            continue
        drug = mp.stem.replace("feature_mart_", "")
        tbl = _ident("mart", drug, a.arm)
        con.execute(f'CREATE TABLE "{tbl}" AS SELECT * FROM read_parquet(?)', [str(mp)])
        cols = [r[1] for r in con.execute(f'PRAGMA table_info("{tbl}")').fetchall()]
        nrow = con.execute(f'SELECT count(*) FROM "{tbl}"').fetchone()[0]
        nfeat = sum(1 for c in cols if c.startswith("raw__"))
        nR = nS = None
        if "y__binary" in cols:
            nR, nS = con.execute(
                f'SELECT count(*) FILTER (WHERE y__binary = 1),'
                f'       count(*) FILTER (WHERE y__binary = 0) FROM "{tbl}"').fetchone()

        meta = _read_meta(mp.with_name(mp.stem + ".metadata.json"))
        cat = meta.get("catalogue")
        if isinstance(cat, dict):
            cat = cat.get("catalogue") or cat.get("catalogue_name") or json.dumps(cat)
        for k, v in meta.items():
            con.execute("INSERT INTO provenance VALUES (?,?,?,?)",
                        [drug, a.arm, k, json.dumps(v) if isinstance(v, (dict, list)) else str(v)])
        con.execute("INSERT INTO marts_index VALUES (?,?,?,?,?,?,?,?,?,?,?)",
                    [tbl, drug, a.arm,
                     str(meta.get("cryptic_version", "")) or None,
                     str(meta.get("mart_version", "")) or None,
                     str(cat) if cat else None,
                     nrow, len(cols), nfeat, nR, nS])
        print(f"  {tbl}: {nrow:,} rows x {len(cols):,} cols ({nfeat:,} features) "
              f"[cryptic {meta.get('cryptic_version', '?')}]")

    # A concordance parquet names neither its drug nor its fold; both live in the
    # selection.json beside it. Unioning the parquets alone produced one table in
    # which four folds were indistinguishable, which makes the artefact ambiguous
    # and the round-trip back to a training run impossible. Pair each parquet
    # with its manifest and carry those two keys as columns.
    by_dir = {Path(m).parent: Path(m) for m in a.manifests if Path(m).exists()}
    con.execute("""CREATE TABLE selection_manifests(
        drug VARCHAR, arm VARCHAR, held_out VARCHAR, manifest JSON)""")

    parts, n_sel = [], 0
    for sp in (Path(x) for x in a.selections):
        if not sp.exists() or not sp.stat().st_size:
            continue
        man = _read_meta(by_dir.get(sp.parent, Path("/nonexistent")))
        drug = str(man.get("drug", "")) or None
        # A selection with no held-out lineage is the refit-on-all arm, so "none"
        # is the right key rather than a guess. NULL would be unusable: the fold
        # is half the primary key a reader needs to find one training unit, and
        # SELECT_FULL's manifest is sometimes the empty {} the module falls back to.
        held = str(man.get("held_out_lineage", "") or "").strip() or "none"
        rel = con.execute("SELECT DISTINCT cryptic_version FROM marts_index "
                          "WHERE cryptic_version IS NOT NULL").fetchall()
        relv = rel[0][0] if len(rel) == 1 else None
        parts.append(
            f"SELECT *, {_q(drug)} AS drug, {_q(a.arm)} AS arm, "
            f"{_q(held)} AS held_out, {_q(relv)} AS cryptic_version "
            f"FROM read_parquet('{sp}')")
        if man:
            con.execute("INSERT INTO selection_manifests VALUES (?,?,?,?)",
                        [drug, a.arm, held, json.dumps(man)])
        n_sel += 1

    if parts:
        con.execute("CREATE TABLE selections AS " + " UNION ALL ".join(parts))
        n = con.execute("SELECT count(*) FROM selections").fetchone()[0]
        folds = con.execute(
            "SELECT count(DISTINCT held_out) FROM selections").fetchone()[0]
        print(f"  selections: {n:,} rows from {n_sel} file(s), {folds} distinct fold(s)")
    else:
        con.execute("""CREATE TABLE selections(
            feature_column VARCHAR, drug VARCHAR, arm VARCHAR, held_out VARCHAR,
            cryptic_version VARCHAR, qf__causal_concordant BOOLEAN)""")
        print("  selections: none supplied")

    con.execute("CHECKPOINT")
    tables = [r[0] for r in con.execute(
        "SELECT table_name FROM duckdb_tables() ORDER BY 1").fetchall()]
    con.close()
    print(f"\nwrote {out}  ({out.stat().st_size / 1e6:.1f} MB, {len(tables)} tables)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
