#!/usr/bin/env python3
"""Attach the per-stage shards and copy their tables into one database.

DuckDB allows only one writer per file, which is why the loads sharded. Here the
shards are attached read-only and a single writer copies them in, so the
concurrency stays where it is safe.

Indexes are NOT carried across by CREATE TABLE AS, so they are re-created from
each shard's own catalogue rather than from a list maintained by hand -- a list
would drift the moment a stage script added an index.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import duckdb


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--shards", nargs="+", type=Path, required=True)
    ap.add_argument("--out", type=Path, required=True)
    a = ap.parse_args()

    if a.out.exists():
        a.out.unlink()
    con = duckdb.connect(str(a.out))

    copied, indexes = [], []
    for i, shard in enumerate(sorted(a.shards)):
        alias = f"s{i}"
        con.execute(f"ATTACH '{shard}' AS {alias} (READ_ONLY)")
        tables = [r[0] for r in con.execute(
            "SELECT table_name FROM information_schema.tables "
            f"WHERE table_catalog='{alias}' AND table_schema='main' "
            "AND table_type='BASE TABLE'").fetchall()]
        for t in tables:
            con.execute(f'CREATE TABLE main."{t}" AS SELECT * FROM {alias}.main."{t}"')
            copied.append(t)
        try:
            indexes += [r[0] for r in con.execute(
                f"SELECT sql FROM duckdb_indexes() WHERE database_name='{alias}'").fetchall()
                if r[0]]
        except duckdb.Error:
            pass  # older DuckDB: indexes are rebuilt by the consumer instead
        con.execute(f"DETACH {alias}")

    for sql in indexes:
        try:
            con.execute(sql)
        except duckdb.Error as exc:
            print(f"  index skipped: {exc}")

    con.execute("CHECKPOINT")
    con.close()
    print(f"assembled {len(copied)} table(s) from {len(a.shards)} shard(s) -> {a.out}")


if __name__ == "__main__":
    main()
