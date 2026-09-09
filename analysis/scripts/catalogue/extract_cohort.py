"""Emit catomatic's two input tables for one drug, from the compendium.

catomatic asks for very little: a samples table of `UNIQUEID, PHENOTYPE` and a
mutations table of `UNIQUEID, MUTATION`. Both vocabularies are already ours --
the projects descend from the same CRyPTIC ecosystem -- so this is two filtered
SELECTs and a provenance record, not a conversion.

The provenance record is the part worth care. A catalogue is only interpretable
against the cohort that produced it, and the cohort here is defined by a
*release tag*, `wgs_samples.dataset`, rather than by a sample list. That tag is
what makes the reproduction checkable: `CRyPTIC-v1.0` recovers 39,545 samples,
39,402 of them phenotyped, against the 39,358 the catomatic paper reports.

Read-support filtering happens HERE and not in catomatic, which has no FRS flag.
The published parameter sweep varies FRS from 0.1 to 0.9 at a fixed background
and p, so the threshold is an input to this step and is recorded as one.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import duckdb


class CohortUnavailable(SystemExit):
    """The database cannot express the requested cohort. Never approximated."""


def _columns(con, table: str) -> set[str]:
    rows = con.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE lower(table_name) = ?", [table.lower()]).fetchall()
    return {str(c).upper() for (c,) in rows}


def _release(con) -> str:
    """The compendium release, from the database's own metadata."""
    try:
        rows = con.execute("SELECT key, value FROM _database_metadata").fetchall()
    except Exception:
        return "unknown"
    meta = {k: v for k, v in rows}
    import re
    v = (meta.get("source_version") or "").strip()
    m = re.search(r"cryptic-tables-(v[0-9]+(?:\.[0-9]+)*)", v)
    if m:
        v = m.group(1)
    elif not v.startswith("v"):
        m = re.search(r"cryptic-tables-(v[0-9]+(?:\.[0-9]+)*)",
                      meta.get("source_dir") or "")
        v = m.group(1) if m else "unknown"
    return f"{v}-slim" if any(k.startswith("slim_") for k in meta) else v


PHENOTYPE_SOURCES = {
    # The DST phenotype: what the catomatic paper's cohort is built on.
    "dst_measurements": {"column": "PHENOTYPE"},
    # The MIC-plate reading the feature marts use. Available deliberately, so a
    # comparison between the two is a run rather than a rewrite.
    "ukmyc_phenotypes": {"column": "BINARY_PHENOTYPE"},
}


def extract(db: Path, drug: str, dataset_tag: str, frs: float,
            out_dir: Path, phenotype_source: str = "dst_measurements",
            genes: list[str] | None = None,
            allow_sparse_frs: bool = False) -> dict:
    con = duckdb.connect(str(db), read_only=True)

    # The cohort tag is a hard requirement, not a preference. CRyPTIC v2.1.2
    # publishes no `dataset` column at all, and a build that quietly fell back
    # to the whole compendium would produce a catalogue that looks like the
    # paper's and answers a different question. Fail, and say which database.
    if phenotype_source not in PHENOTYPE_SOURCES:
        raise CohortUnavailable(f"unknown phenotype source {phenotype_source!r}")
    if not _columns(con, phenotype_source):
        raise CohortUnavailable(
            f"{db.name} has no {phenotype_source} table; the slim build drops "
            "several phenotype tables. Build the catalogue against a full "
            "compendium.")
    if "DATASET" not in _columns(con, "wgs_samples"):
        raise CohortUnavailable(
            f"{db.name} has no wgs_samples.dataset column, so the cohort "
            f"{dataset_tag!r} cannot be selected. This is a property of the "
            "CRyPTIC release, not of the build: v2.1.2 does not publish it. "
            "Build the catalogue against a v3.4.0 compendium, or supply a "
            "cohort by another means and change this script deliberately."
        )

    tags = [t for (t,) in con.execute(
        "SELECT DISTINCT dataset FROM wgs_samples WHERE dataset IS NOT NULL"
    ).fetchall()]
    if dataset_tag not in tags:
        raise CohortUnavailable(
            f"cohort tag {dataset_tag!r} is not in {db.name}. Available: "
            f"{sorted(tags)}")

    out_dir.mkdir(parents=True, exist_ok=True)
    samples_path = out_dir / f"samples_{drug}.csv"
    mutations_path = out_dir / f"mutations_{drug}.csv"

    con.execute("""
        CREATE OR REPLACE TEMP TABLE cohort AS
        SELECT ws.UNIQUEID
        FROM   wgs_samples ws
        WHERE  ws.dataset = ?
    """, [dataset_tag])

    # One phenotype per sample per drug, R or S only. catomatic's binary builder
    # has no vocabulary for an intermediate call, so an ambiguous phenotype is
    # dropped here rather than coerced into one of the two it understands.
    #
    # The source is dst_measurements, NOT ukmyc_phenotypes. Both tables carry a
    # binary call and either will produce a catalogue, which is what makes the
    # choice dangerous: on this cohort ukmyc_phenotypes yields 12,324 phenotyped
    # samples and dst_measurements yields 39,402 -- the latter being the figure
    # the catomatic paper's cohort matches. ukmyc_phenotypes is the MIC-plate
    # reading used by the feature marts; dst_measurements is the drug
    # susceptibility test. Picking the wrong one produces a catalogue built from
    # a third of the intended cohort, and nothing downstream would say so.
    src = PHENOTYPE_SOURCES[phenotype_source]

    # A sample can carry MORE THAN ONE measurement for a drug, and they do not
    # always agree: 271 samples in this cohort have both an R and an S row for
    # rifampicin. Grouping by (sample, phenotype) emits such a sample twice,
    # once as resistant and once as susceptible, and catomatic would count it on
    # both sides of every test it runs. Contradictory samples are DROPPED rather
    # than resolved by a rule nobody stated, and the count is reported.
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE pheno AS
        SELECT ph.UNIQUEID,
               count(DISTINCT ph.{src['column']}) AS n_distinct,
               min(ph.{src['column']})            AS PHENOTYPE
        FROM   {phenotype_source} ph
        JOIN   cohort USING (UNIQUEID)
        WHERE  ph.DRUG = ?
          AND  ph.{src['column']} IN ('R','S')
        GROUP BY ph.UNIQUEID
    """, [drug])
    n_conflicting = con.execute(
        "SELECT count(*) FROM pheno WHERE n_distinct > 1").fetchone()[0]
    con.execute(f"""
        COPY (
            SELECT UNIQUEID, PHENOTYPE FROM pheno WHERE n_distinct = 1
        ) TO '{samples_path}' (FORMAT CSV, HEADER)
    """)

    # Gene restriction is not an optimisation. The published per-drug catalogues
    # are built over the drug's own genes -- rifampicin over rpoB, bedaquiline
    # over Rv0678, atpE and pepQ -- and an unrestricted build grades the whole
    # genome, including the PE/PPE families no catalogue entry mentions. On this
    # cohort that is 1,801 genes and ~994k mutation rows against a few thousand,
    # and the resulting catalogue is not comparable with the reference at all.
    # FRS coverage is checked BEFORE it is used as a filter. In the v3.4.0
    # build the column is populated for 10.19% of mutation rows (v2.1.2: 0.33%),
    # so `FRS >= 0.9` silently discards ~96% of rpoB's rows -- 78,542 down to
    # 60 -- and still returns a catalogue. Five such catalogues, one per sweep
    # setting, would each be built from a tenth of the cohort and none would
    # say so. A sparse column is refused rather than applied.
    gene_filter = ""
    cov_params: list = []
    if genes:
        gene_filter = "WHERE m.GENE IN (" + ", ".join("?" for _ in genes) + ")"
        cov_params = list(genes)
    total, with_frs = con.execute(
        f"SELECT count(*), count(m.FRS) FROM mutations m {gene_filter}",
        cov_params).fetchone()
    frs_coverage = round(100.0 * with_frs / total, 2) if total else 0.0
    if frs is not None and frs_coverage < 50.0 and not allow_sparse_frs:
        raise CohortUnavailable(
            f"FRS is populated for only {frs_coverage}% of the {total} mutation "
            f"rows in scope, so filtering at FRS >= {frs} would drop most of the "
            "cohort and still produce a catalogue. Build without --frs, or pass "
            "--allow-sparse-frs if the sparsity is itself the subject."
        )

    frs_clause = "AND m.FRS >= ?" if frs is not None else ""
    params = [frs] if frs is not None else []
    gene_clause = ""
    if genes:
        placeholders = ", ".join("?" for _ in genes)
        gene_clause = f"AND m.GENE IN ({placeholders})"
        params = params + list(genes)
    con.execute(f"""
        COPY (
            SELECT m.UNIQUEID, m.GENE || '@' || m.MUTATION AS MUTATION
            FROM   mutations m
            JOIN   cohort USING (UNIQUEID)
            WHERE  m.UNIQUEID IN (SELECT UNIQUEID FROM read_csv_auto('{samples_path}'))
              {frs_clause}
              {gene_clause}
        ) TO '{mutations_path}' (FORMAT CSV, HEADER)
    """, params)

    n_cohort = con.execute("SELECT count(*) FROM cohort").fetchone()[0]
    # The paper's 39,358 is a COHORT figure -- samples carrying an R/S
    # phenotype for any drug -- not a per-drug one. Comparing this drug's
    # sample count against it produced a -68.9% "delta" on the first run, which
    # is not a discrepancy but a category error. Both counts are reported, and
    # only the cohort-level one is ever compared.
    n_cohort_phenotyped = con.execute(f"""
        SELECT count(DISTINCT ph.UNIQUEID)
        FROM   {phenotype_source} ph
        JOIN   cohort USING (UNIQUEID)
        WHERE  ph.{src['column']} IN ('R','S')
    """).fetchone()[0]
    n_samples = con.execute(
        f"SELECT count(*) FROM read_csv_auto('{samples_path}')").fetchone()[0]
    n_mutations = con.execute(
        f"SELECT count(*) FROM read_csv_auto('{mutations_path}')").fetchone()[0]
    n_r = con.execute(
        f"SELECT count(*) FROM read_csv_auto('{samples_path}') "
        "WHERE PHENOTYPE = 'R'").fetchone()[0]

    report = {
        "drug": drug,
        "cryptic_version": _release(con),
        "dataset_tag": dataset_tag,
        "phenotype_source": phenotype_source,
        "genes_requested": genes,
        "frs_threshold": frs,
        "frs_coverage_pct": frs_coverage,
        "n_cohort_samples": n_cohort,
        "n_cohort_phenotyped_any_drug": n_cohort_phenotyped,
        "n_phenotyped_this_drug": n_samples,
        "n_dropped_conflicting_phenotype": n_conflicting,
        "n_resistant": n_r,
        "n_susceptible": n_samples - n_r,
        "n_mutation_rows": n_mutations,
        "samples_file": samples_path.name,
        "mutations_file": mutations_path.name,
        # The paper's own count for the training cohort, carried so a reader of
        # the report can see the delta without going to the publication. The
        # paper notes its figure is itself seven fewer than the list used to
        # build WHOv1, "the reason for this is unknown", so agreement here is
        # stated as a delta and never as identity.
        "reference_n_phenotyped": 39358 if dataset_tag == "CRyPTIC-v1.0" else None,
    }
    ref = report["reference_n_phenotyped"]
    if ref:
        report["delta_vs_reference"] = n_cohort_phenotyped - ref
        report["delta_pct"] = round(100.0 * (n_cohort_phenotyped - ref) / ref, 3)
        report["delta_is_over"] = "cohort phenotyped for any drug"

    # ---- wildcards -------------------------------------------------------
    # catomatic's piezo export requires a wildcards file: the GARC1 rules a
    # catalogue carries beyond its explicit entries. catomatic ships a TEMPLATE
    # whose key is the literal string "gene", not a gene list, so it cannot be
    # used as-is against a real cohort.
    #
    # The rules below are catomatic's own defaults, verbatim. What is derived
    # here is only the GENE LIST they are expanded over, taken from the genes
    # actually present in this cohort rather than from a panel written
    # elsewhere. The distinction is recorded in the report: the rules are
    # upstream's, the scope is this cohort's. Supply --wildcards to override
    # both, which is the right move if the grading rules are ever a finding
    # rather than a default.
    genes = [g for (g,) in con.execute(
        f"SELECT DISTINCT split_part(MUTATION, '@', 1) "
        f"FROM read_csv_auto('{mutations_path}') ORDER BY 1").fetchall() if g]
    rules = {
        "*=": {"pred": "S", "evid": {"default_rule": "True"}},
        "-*_indel": {"pred": "U", "evid": {"default_rule": "True"}},
        "*_indel": {"pred": "U", "evid": {"default_rule": "True"}},
        "-*?": {"pred": "U", "evid": {"default_rule": "True"}},
        "*?": {"pred": "U", "evid": {"default_rule": "True"}},
        "del_0.0": {"pred": "U", "evid": {"default_rule": "True"}},
    }
    wildcards = {f"{g}@{suffix}": rule
                 for g in genes for suffix, rule in rules.items()}
    (out_dir / f"wildcards_{drug}.json").write_text(
        json.dumps(wildcards, indent=2) + "\n")
    report["n_genes"] = len(genes)
    report["wildcards_file"] = f"wildcards_{drug}.json"
    report["wildcards_origin"] = (
        "catomatic default rules, expanded over the genes present in this "
        "cohort; not a curated rule set")

    (out_dir / f"cohort_{drug}.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--db", type=Path, required=True)
    p.add_argument("--drug", required=True)
    p.add_argument("--dataset-tag", default="CRyPTIC-v1.0",
                   help="wgs_samples.dataset value defining the cohort")
    p.add_argument("--phenotype-source", default="dst_measurements",
                   choices=sorted(PHENOTYPE_SOURCES),
                   help="dst_measurements is the DST phenotype the catomatic "
                        "cohort uses; ukmyc_phenotypes is the MIC-plate reading "
                        "the feature marts use, and yields a third of the cohort")
    p.add_argument("--genes", default=None,
                   help="comma-separated gene restriction, e.g. rpoB. Omit to "
                        "grade every gene, which is not what the published "
                        "per-drug catalogues do")
    p.add_argument("--frs", type=float, default=None,
                   help="minimum read support; omit for no FRS filter")
    p.add_argument("--allow-sparse-frs", action="store_true",
                   help="apply an FRS filter even where the column is mostly "
                        "null; the default refuses, because the result looks "
                        "like a catalogue either way")
    p.add_argument("--out", type=Path, default=Path("."))
    a = p.parse_args()
    genes = [g.strip() for g in a.genes.split(",")] if a.genes else None
    r = extract(a.db, a.drug, a.dataset_tag, a.frs, a.out, a.phenotype_source, genes,
                a.allow_sparse_frs)
    print(f"{r['drug']:<5} {r['cryptic_version']:<12} cohort={r['n_cohort_samples']} "
          f"phenotyped(any)={r['n_cohort_phenotyped_any_drug']} "
          f"phenotyped({r['drug']})={r['n_phenotyped_this_drug']} (R={r['n_resistant']}) "
          f"mutations={r['n_mutation_rows']} frs={r['frs_threshold']}"
          + (f" delta_vs_paper={r['delta_vs_reference']:+}"
             if r.get("delta_vs_reference") is not None else ""))


if __name__ == "__main__":
    main()
