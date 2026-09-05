"""Generate the *ideal serving input* from the CRyPTIC DB (design/16, Stage 0).

The product goal (design/16 §Goal) is: a user hands us a file, we return DR
predictions + reasoning. This script defines and materialises the file we *want*
to receive, straight from real crypticdb data — so the contract is grounded in
data that exists, and each generated instance doubles as the Stage-1 round-trip
fixture (DB -> ideal input -> vectorise -> model == mart-fed training prediction).

What the ideal input IS (design/16):
  the CRyPTIC/gnomonicus GENOTYPE (variants, GARC) + sample QC covariates.

What it deliberately is NOT:
  any catalogue interpretation — no `effects`, no `predictions`, no R/S calls.
  Producing that interpretation is the model's job; shipping it in the input
  would leak the label. So this generator reads only `genomes` (covariates) and
  `mutations` (genotype); it never touches `effects`/`predictions`.

Absence = reference: only called variants are listed. Any panel mutation not
present is treated as wild-type (0) by the vectoriser.

Output: one self-describing JSON per sample, conforming to
`input_spec.schema.json` (the committed contract). Validated on write — with
`jsonschema` if importable, else a built-in structural check.

Run (pixi env; read-only DB):
    .pixi/envs/default/bin/python -m analysis.scripts.ingest.reference_input \
        --db analysis/databases/duckdb/cryptic-slim.duckdb \
        --auto --out analysis/results/ingest/examples
    # or a specific isolate:
    ... --uniqueid site.04.subj.01561.lab.729880.iso.1
"""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path

import duckdb

SCHEMA_VERSION = "1.0.0"
SCHEMA_PATH = Path(__file__).with_name("input_spec.schema.json")


def connect(db_path) -> duckdb.DuckDBPyConnection:
    """Open the CRyPTIC DB read-only.

    Historical note: cryptic-slim.duckdb once had corrupt row-group statistics
    that made the optimizer wrongly prune on single-`UNIQUEID` filters (returning
    a truncated genotype — 506 rows vs the true 1301). The `mutations` table was
    rebuilt to regenerate correct statistics (design/16 §Data-integrity), so the
    `PRAGMA disable_optimizer` workaround is no longer needed and pruning is both
    correct and fast (~40 ms per single-`UNIQUEID` scan of 60M rows).
    """
    return duckdb.connect(str(db_path), read_only=True)


def _clean_float(x) -> float | None:
    """DuckDB/pandas hand back NaN for missing FRS; JSON can't hold NaN."""
    if x is None:
        return None
    x = float(x)
    return None if math.isnan(x) else x


def _assemble(uniqueid, lineage, sublineage, species, tb_depth, tb_breadth, variant_rows) -> dict:
    variants = [
        {"gene": g, "mutation": m, "frs": _clean_float(frs),
         "coverage": _clean_float(cov_), "is_minor": bool(is_minor)}
        for (g, m, frs, cov_, is_minor) in variant_rows
    ]
    return {
        "schema_version": SCHEMA_VERSION,
        "sample_id": uniqueid,
        "nomenclature": "GARC",
        "source": {"producer": "cryptic-slim", "reference": "NC_000962.3"},
        "covariates": {
            "lineage": lineage,
            "median_coverage": _clean_float(tb_depth),
            "breadth": _clean_float(tb_breadth),
            "sublineage": sublineage,
            "species": species,
        },
        "variants": variants,
    }


def build_input(con: duckdb.DuckDBPyConnection, uniqueid: str) -> dict:
    """Assemble the ideal-input dict for one isolate from genomes + mutations."""
    cov = con.execute(
        "SELECT LINEAGE, SUBLINEAGE, SPECIES, TB_DEPTH, TB_COVERAGE "
        "FROM genomes WHERE UNIQUEID = ?", [uniqueid]
    ).fetchone()
    if cov is None:
        raise SystemExit(f"no genome row for {uniqueid!r}")
    # Deterministic variant order so the fixture is reproducible byte-for-byte.
    rows = con.execute(
        "SELECT GENE, MUTATION, FRS, COVERAGE, IS_MINOR "
        "FROM mutations WHERE UNIQUEID = ? ORDER BY GENE, MUTATION", [uniqueid]
    ).fetchall()
    return _assemble(uniqueid, *cov, rows)


def build_inputs_batch(con: duckdb.DuckDBPyConnection, uniqueids: list[str]) -> dict[str, dict]:
    """Assemble many isolates in ONE scan each of genomes/mutations.

    With the optimizer disabled (see `connect`), a single `UNIQUEID IN (...)` scan
    returns all rows correctly — vastly cheaper than one full scan per isolate,
    which is what makes the round-trip over dozens of isolates feasible.
    """
    from collections import defaultdict
    ph = ",".join(["?"] * len(uniqueids))
    covs = {r[0]: r[1:] for r in con.execute(
        f"SELECT UNIQUEID, LINEAGE, SUBLINEAGE, SPECIES, TB_DEPTH, TB_COVERAGE "
        f"FROM genomes WHERE UNIQUEID IN ({ph})", uniqueids).fetchall()}
    mut_by: dict[str, list] = defaultdict(list)
    for uid, g, m, frs, cov_, minor in con.execute(
        f"SELECT UNIQUEID, GENE, MUTATION, FRS, COVERAGE, IS_MINOR "
        f"FROM mutations WHERE UNIQUEID IN ({ph}) ORDER BY UNIQUEID, GENE, MUTATION",
        uniqueids).fetchall():
        mut_by[uid].append((g, m, frs, cov_, minor))
    return {uid: _assemble(uid, *covs[uid], mut_by[uid]) for uid in uniqueids if uid in covs}


def _pick_examples(con: duckdb.DuckDBPyConnection) -> list[tuple[str, str]]:
    """Two grounded, deterministic exemplars: an MDR case (rpoB S450L + katG
    S315T -> RIF-R/INH-R) and a pan-susceptible case (no R prediction)."""
    mdr = con.execute("""
        SELECT g.UNIQUEID FROM genomes g
        WHERE g.LINEAGE IS NOT NULL AND g.TB_DEPTH IS NOT NULL AND g.TB_COVERAGE IS NOT NULL
          AND EXISTS (SELECT 1 FROM mutations m WHERE m.UNIQUEID=g.UNIQUEID AND m.GENE='rpoB' AND m.MUTATION='S450L')
          AND EXISTS (SELECT 1 FROM mutations m WHERE m.UNIQUEID=g.UNIQUEID AND m.GENE='katG' AND m.MUTATION='S315T')
        ORDER BY g.UNIQUEID LIMIT 1""").fetchone()
    pans = con.execute("""
        SELECT g.UNIQUEID FROM genomes g
        WHERE g.LINEAGE IS NOT NULL AND g.TB_DEPTH IS NOT NULL AND g.TB_COVERAGE IS NOT NULL
          AND NOT EXISTS (SELECT 1 FROM predictions p WHERE p.UNIQUEID=g.UNIQUEID AND p.PREDICTION='R')
          AND (SELECT count(*) FROM mutations m WHERE m.UNIQUEID=g.UNIQUEID) BETWEEN 2 AND 6
        ORDER BY g.UNIQUEID LIMIT 1""").fetchone()
    out = []
    if mdr:
        out.append(("mdr", mdr[0]))
    if pans:
        out.append(("pan_susceptible", pans[0]))
    return out


# --------------------------------------------------------------------------- #
# validation against the committed schema
# --------------------------------------------------------------------------- #
def validate(instance: dict) -> None:
    schema = json.loads(SCHEMA_PATH.read_text())
    try:
        import jsonschema  # optional — real Draft 2020-12 validation if present
        jsonschema.validate(instance, schema)
        return
    except ImportError:
        _structural_check(instance, schema)


def _structural_check(instance: dict, schema: dict) -> None:
    """Dependency-free fallback: required keys + top-level types + per-variant
    required keys. Not a full JSON-Schema validator, but catches contract drift."""
    for key in schema["required"]:
        if key not in instance:
            raise ValueError(f"missing required top-level key: {key!r}")
    if instance["schema_version"] != SCHEMA_VERSION:
        raise ValueError(f"schema_version {instance['schema_version']!r} != {SCHEMA_VERSION!r}")
    if instance["nomenclature"] != "GARC":
        raise ValueError("nomenclature must be 'GARC'")
    for key in ("lineage", "median_coverage", "breadth"):
        if key not in instance["covariates"]:
            raise ValueError(f"missing required covariate: {key!r}")
    for i, v in enumerate(instance["variants"]):
        for key in ("gene", "mutation"):
            if key not in v:
                raise ValueError(f"variant[{i}] missing required key: {key!r}")


def main() -> None:
    p = argparse.ArgumentParser(description="Generate the ideal serving input from the CRyPTIC DB (design/16 Stage 0).")
    p.add_argument("--db", type=Path, default=Path("analysis/databases/duckdb/cryptic-slim.duckdb"))
    p.add_argument("--uniqueid", action="append", default=[], help="isolate(s) to emit; repeatable.")
    p.add_argument("--auto", action="store_true", help="auto-pick an MDR + a pan-susceptible exemplar.")
    p.add_argument("--out", type=Path, default=Path("analysis/results/ingest/examples"))
    args = p.parse_args()

    args.out.mkdir(parents=True, exist_ok=True)
    con = connect(args.db)

    targets: list[tuple[str, str]] = [("sample", u) for u in args.uniqueid]
    if args.auto:
        targets += _pick_examples(con)
    if not targets:
        raise SystemExit("nothing to do — pass --auto or --uniqueid")

    for label, uid in targets:
        inst = build_input(con, uid)
        validate(inst)
        fname = f"{label}.input.json" if label != "sample" else f"{uid}.input.json"
        (args.out / fname).write_text(json.dumps(inst, indent=2))
        print(f"[{label}] {uid}: {len(inst['variants'])} variants, "
              f"lineage={inst['covariates']['lineage']}, "
              f"depth={inst['covariates']['median_coverage']:.1f} -> {fname}")


if __name__ == "__main__":
    main()
