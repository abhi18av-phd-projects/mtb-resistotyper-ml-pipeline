#!/usr/bin/env python3
"""Referential-integrity gate, scoped to the tables the analysis actually reads.

Two corrections to the first version of this gate.

It queried a table called `ri_summary`, which does not exist. Stage 11 creates
`v_ri_*` VIEWS and prints a summary to stdout; the table name was invented, so
the gate would have failed every run for the wrong reason. The checks are now
discovered from the catalogue.

And it blocked on any orphan at all. On the v3.4.0 build that means failing on
43,872 orphans in `ukmyc_growth_to_plates` and 1,190 in `bashthebug_to_plates`
-- plate-level readings and citizen-science classifications that no feature mart
reads. A gate that blocks on findings outside what it protects is a gate people
switch off, which is worse than a narrow one. So orphans in the analysis path
FAIL, and orphans elsewhere are REPORTED.

The scope is derived from the view name rather than listed by hand: each view is
`v_ri_<source>_to_<target>`, and a check blocks when <source> is a table the
marts read. A new check added to stage 11 is therefore classified automatically.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import duckdb

# Read empirically from cohort.py and build_mart.py, not from memory:
#   grep -ohE '\b(FROM|JOIN)\s+[a-zA-Z_]+' cohort.py build_mart.py
ANALYSIS_TABLES = {
    "effects", "genomes", "mutations", "predictions",
    "ukmyc_phenotypes", "wgs_samples",
    # lookups the mart joins through
    "drugs", "countries", "sites",
}

# Stage 11 abbreviates some source names in its view names.
SOURCE_ALIASES = {
    "ukmyc_phen": "ukmyc_phenotypes",
    "dstmeas": "dst_measurements",
}


def source_table(view: str) -> str:
    body = view[len("v_ri_"):]
    src = body.rsplit("_to_", 1)[0] if "_to_" in body else body
    return SOURCE_ALIASES.get(src, src)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--db", type=Path, required=True)
    ap.add_argument("--out", type=Path, default=Path("ri_summary.tsv"))
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--strict", action="store_true",
                    help="fail on advisory orphans too")
    a = ap.parse_args()

    con = duckdb.connect(str(a.db), read_only=True)
    views = [r[0] for r in con.execute(
        "SELECT table_name FROM information_schema.tables "
        "WHERE table_type='VIEW' AND table_name LIKE 'v_ri_%' ORDER BY 1").fetchall()]
    if not views:
        sys.exit("no v_ri_* views found — did stage 11_ri_checks run?")

    rows = []
    for v in views:
        src = source_table(v)
        n = con.execute(f'SELECT COUNT(*) FROM "{v}"').fetchone()[0]
        rows.append({"check": v[len("v_ri_"):], "view": v, "source": src,
                     "orphans": int(n),
                     "scope": "analysis" if src in ANALYSIS_TABLES else "peripheral"})
    con.close()

    blocking = [r for r in rows if r["scope"] == "analysis" and r["orphans"]]
    advisory = [r for r in rows if r["scope"] == "peripheral" and r["orphans"]]

    with a.out.open("w") as fh:
        fh.write("check\tsource\tscope\torphans\n")
        for r in sorted(rows, key=lambda r: (-r["orphans"], r["check"])):
            fh.write(f'{r["check"]}\t{r["source"]}\t{r["scope"]}\t{r["orphans"]}\n')
    if a.json:
        a.json.write_text(json.dumps(
            {"checks": rows, "n_blocking": len(blocking),
             "n_advisory": len(advisory)}, indent=2) + "\n")

    n_analysis = sum(1 for r in rows if r["scope"] == "analysis")
    print(f"{len(rows)} checks — {n_analysis} in the analysis path, "
          f"{len(rows) - n_analysis} peripheral")

    if advisory:
        print("\nADVISORY — orphans outside the analysis path, not blocking:")
        for r in sorted(advisory, key=lambda r: -r["orphans"]):
            print(f"  {r['check']:32} {r['orphans']:>9,}  ({r['source']})")

    if blocking or (a.strict and advisory):
        print("\nFAILED — orphans in the analysis path:", file=sys.stderr)
        for r in sorted(blocking or advisory, key=lambda r: -r["orphans"]):
            print(f"  {r['check']:32} {r['orphans']:>9,}  ({r['source']})",
                  file=sys.stderr)
        sys.exit(1)

    print("\nOK — every check on the analysis path is clean.")


if __name__ == "__main__":
    main()
