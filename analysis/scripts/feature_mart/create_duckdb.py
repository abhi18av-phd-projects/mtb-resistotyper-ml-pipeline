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

try:
    from analysis.scripts.feature_mart import fe_identity
except ImportError:  # run as a plain script from this directory
    import fe_identity

SAFE = re.compile(r"[^A-Za-z0-9_]")


def _ident(*parts: str) -> str:
    return SAFE.sub("_", "_".join(p for p in parts if p)).strip("_")


def _q(v) -> str:
    """SQL literal, or NULL."""
    return "NULL" if v is None else "'" + str(v).replace("'", "''") + "'"


def _unit_of(p) -> tuple[str, str]:
    """(drug, fold) from `concordance_<DRUG>_<fold>.parquet` / `selection_<DRUG>_<fold>.json`.

    Both come from the filename. Naming by fold alone was unique for one drug and
    collided for two: CREATE_DUCKDB collects every drug's selections into one
    task, so RIF and INH each staged a `concordance_lineage1.parquet`. Neither a
    drug code (RIF, BDQ) nor a fold (lineage1..lineage4, none) contains an
    underscore, so the LAST underscore separates them.

    A legacy fold-only name (`concordance_lineage1`) yields drug "" so the
    caller can fall back to the manifest; an unrecognised name yields ("", ""),
    so an unknown fold is never mistaken for the refit-on-all arm.
    """
    stem = Path(p).stem
    for prefix in ("concordance_", "selection_"):
        if stem.startswith(prefix):
            rest = stem[len(prefix):]
            if "_" in rest:
                drug, fold = rest.rsplit("_", 1)
                return drug, fold
            return "", rest
    return "", ""


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
                    help="selection_<DRUG>_<fold>.json manifests, paired with "
                         "--selections by (drug, fold) parsed from the filename.")
    ap.add_argument("--out", required=True, help="DuckDB file to write")
    ap.add_argument("--arm", default="", help="campaign arm this run belongs to")
    ap.add_argument("--fe-pre-steps", default="", help="label-free FE steps that ran")
    ap.add_argument("--fe-fold-steps", default="", help="label-using FE steps that ran")
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
        n_resistant BIGINT, n_susceptible BIGINT,
        fe_version VARCHAR, fe_id VARCHAR, mart_id VARCHAR)""")
    # The FE technique behind each mart, once per fe_id. fe_id excludes the
    # database, so marts of one technique on two releases share a row here and
    # differ only in mart_id -- the join a cross-release comparison needs.
    con.execute("""CREATE TABLE fe_identity(
        fe_id VARCHAR, fe_version VARCHAR, fe_config JSON)""")
    fe_seen: set[str] = set()
    # Selection manifests by (drug, fold), read up front: a mart's FE identity
    # includes the selection settings, so it is needed while indexing marts.
    by_unit = {_unit_of(m): Path(m) for m in a.manifests if Path(m).exists()}
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
        # Prefer the refit-on-all selection; any fold's settings are the same.
        sel = next((by_unit[k] for k in ((drug, "none"),) if k in by_unit), None) \
            or next((v for (d, _f), v in by_unit.items() if d == drug), None)
        try:
            ident = fe_identity.identity(meta, _read_meta(sel) if sel else {},
                                         a.fe_pre_steps, a.fe_fold_steps)
            fe_v, fe_i = ident["fe_version"], ident["fe_id"]
            m_i = fe_identity.mart_id(fe_i, (meta.get("source_db") or {}).get("checksum_sha256"))
            if fe_i not in fe_seen:
                con.execute("INSERT INTO fe_identity VALUES (?,?,?)",
                            [fe_i, fe_v, json.dumps(ident["fe_config"])])
                fe_seen.add(fe_i)
        except SystemExit as exc:
            # A mart built before the FE block existed has no recorded identity,
            # and this process cannot vouch for the code that built it. NULL is
            # the honest value; the deposit says so rather than guessing.
            print(f"  {drug}: no FE identity ({exc})")
            fe_v = fe_i = m_i = None
        con.execute("INSERT INTO marts_index VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    [tbl, drug, a.arm,
                     str(meta.get("cryptic_version", "")) or None,
                     str(meta.get("mart_version", "")) or None,
                     str(cat) if cat else None,
                     nrow, len(cols), nfeat, nR, nS, fe_v, fe_i, m_i])
        print(f"  {tbl}: {nrow:,} rows x {len(cols):,} cols ({nfeat:,} features) "
              f"[cryptic {meta.get('cryptic_version', '?')}]")

    # A concordance parquet names neither its drug nor its fold; both live in the
    # selection.json beside it. Unioning the parquets alone produced one table in
    # which four folds were indistinguishable, which makes the artefact ambiguous
    # and the round-trip back to a training run impossible. Pair each parquet
    # with its manifest and carry those two keys as columns.
    # Pair by (DRUG, FOLD), taken from the filename, not by directory and not by
    # fold alone. Directory pairing broke the moment the two inputs were staged
    # into separate per-file dirs: every lookup missed, every selection silently
    # fell back to drug=NULL and held_out='none', and five folds collapsed into
    # one indistinguishable table. Fold-only pairing then held for exactly one
    # drug -- with two, RIF's lineage1 manifest and INH's share a key and one
    # silently overwrites the other.
    by_unit = {_unit_of(m): Path(m) for m in a.manifests if Path(m).exists()}
    con.execute("""CREATE TABLE selection_manifests(
        drug VARCHAR, arm VARCHAR, held_out VARCHAR, manifest JSON)""")

    # One campaign is one drug in practice, and marts_index already knows it. Use
    # it when a manifest is absent, so the selections table is never keyed on NULL.
    drugs = [r[0] for r in con.execute(
        "SELECT DISTINCT drug FROM marts_index WHERE drug IS NOT NULL").fetchall()]
    only_drug = drugs[0] if len(drugs) == 1 else None

    parts, n_sel = [], 0
    for sp in (Path(x) for x in a.selections):
        if not sp.exists() or not sp.stat().st_size:
            continue
        name_drug, fold = _unit_of(sp)
        man = _read_meta(by_unit.get((name_drug, fold), Path("/nonexistent")))
        # The filename is authoritative for the drug, as it is for the fold. The
        # manifest and the single-drug fallback only serve legacy fold-only names.
        drug = name_drug or str(man.get("drug", "")) or only_drug
        # A selection with no held-out lineage is the refit-on-all arm, so "none"
        # is the right key rather than a guess. NULL would be unusable: the fold
        # is half the primary key a reader needs to find one training unit, and
        # SELECT_FULL's manifest is sometimes the empty {} the module falls back to.
        # The filename is authoritative for the fold; the manifest only enriches it.
        # A manifest that failed to arrive must not silently relabel a fold "none".
        held = fold or str(man.get("held_out_lineage", "") or "").strip() or "none"
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
