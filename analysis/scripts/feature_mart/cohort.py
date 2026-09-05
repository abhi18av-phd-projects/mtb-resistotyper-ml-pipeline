"""Path-1 cohort builder for the slim DB model-vs-WHO-catalogue comparison.

Realises §9 step 1 of `abc-universe/brainstorms/mtb-resistotyper-ml/
2026-06-30-slim-db-catalogue-comparison-handoff.md`: the labelled
substrate that the per-drug feature marts (`design/09`) and the modelling
pipeline (`design/03`) consume.

Output: one row per (UNIQUEID, DRUG) with
    truth (DST R/S), catalogue (WHO v2 R/S), lineage, sublineage,
    phenotype_quality

Joins (per handoff §3):
    ukmyc_phenotypes ⋈ predictions ⋈ genomes
    on UNIQUEID and (UNIQUEID, DRUG) — DST quality keeps R/S only;
    catalogue scope keeps R/S only (U/F excluded); catalogue version
    pinned to WHO-UCN-GTB-PCI-2023.5 (the slim DB's bundled v2).

Run:
    python -m analysis.scripts.feature_mart.cohort \
        --db analysis/databases/duckdb/cryptic-slim.duckdb \
        --out analysis/results/cohort/

Outputs:
    cohort.parquet                     full labelled cohort
    per_drug_summary.csv               n_total, n_R/S, %R, catalogue calls
    per_drug_lineage_summary.csv       fairness power table (L1–L4 + minor)
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import duckdb


CATALOGUE_NAME = "WHO-UCN-GTB-PCI-2023.5"


COHORT_SQL = f"""
SELECT ph.UNIQUEID,
       ph.DRUG,
       ph.BINARY_PHENOTYPE   AS truth,
       pr.PREDICTION         AS catalogue,
       g.LINEAGE,
       g.SUBLINEAGE,
       ph.PHENOTYPE_QUALITY
FROM   ukmyc_phenotypes ph
JOIN   predictions       pr USING (UNIQUEID, DRUG)
JOIN   genomes           g  USING (UNIQUEID)
WHERE  ph.BINARY_PHENOTYPE IN ('R','S')
  AND  pr.PREDICTION       IN ('R','S')
  AND  pr.CATALOGUE_NAME = '{CATALOGUE_NAME}'
"""


PER_DRUG_SQL = """
SELECT DRUG,
       COUNT(*)                                           AS n_total,
       SUM(CASE WHEN truth = 'R' THEN 1 ELSE 0 END)       AS n_R,
       SUM(CASE WHEN truth = 'S' THEN 1 ELSE 0 END)       AS n_S,
       ROUND(100.0 * SUM(CASE WHEN truth='R' THEN 1 ELSE 0 END)
                       / COUNT(*), 2)                     AS pct_R,
       SUM(CASE WHEN catalogue = 'R' THEN 1 ELSE 0 END)   AS cat_called_R,
       SUM(CASE WHEN catalogue = 'R' AND truth = 'S' THEN 1 ELSE 0 END)
                                                          AS cat_FP,
       SUM(CASE WHEN catalogue = 'S' AND truth = 'R' THEN 1 ELSE 0 END)
                                                          AS cat_FN
FROM cohort
GROUP BY DRUG
ORDER BY n_total DESC
"""

PER_DRUG_LINEAGE_SQL = """
SELECT DRUG, LINEAGE,
       COUNT(*)                                           AS n,
       SUM(CASE WHEN truth = 'R' THEN 1 ELSE 0 END)       AS n_R
FROM cohort
WHERE LINEAGE IS NOT NULL
GROUP BY DRUG, LINEAGE
ORDER BY DRUG, LINEAGE
"""


def build_cohort(db_path: Path, out_dir: Path) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    con = duckdb.connect(str(db_path), read_only=True)

    con.execute(f"CREATE OR REPLACE TEMP VIEW cohort AS {COHORT_SQL}")

    cohort_parquet = out_dir / "cohort.parquet"
    con.execute(
        f"COPY (SELECT * FROM cohort) "
        f"TO '{cohort_parquet}' (FORMAT 'parquet', COMPRESSION 'zstd')"
    )

    per_drug_csv = out_dir / "per_drug_summary.csv"
    rows = con.execute(PER_DRUG_SQL).fetchall()
    cols = [d[0] for d in con.description]
    with per_drug_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)

    per_lineage_csv = out_dir / "per_drug_lineage_summary.csv"
    rows = con.execute(PER_DRUG_LINEAGE_SQL).fetchall()
    cols = [d[0] for d in con.description]
    with per_lineage_csv.open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(cols)
        w.writerows(rows)

    counts = con.execute(
        "SELECT COUNT(*) AS rows, "
        "       COUNT(DISTINCT UNIQUEID) AS isolates, "
        "       COUNT(DISTINCT DRUG) AS drugs "
        "FROM cohort"
    ).fetchone()

    manifest = {
        "source_db": str(db_path),
        "catalogue_name": CATALOGUE_NAME,
        "cohort_rows": counts[0],
        "distinct_isolates": counts[1],
        "distinct_drugs": counts[2],
        "outputs": {
            "cohort_parquet": str(cohort_parquet.relative_to(out_dir.parent.parent)),
            "per_drug_summary": str(per_drug_csv.relative_to(out_dir.parent.parent)),
            "per_drug_lineage_summary": str(per_lineage_csv.relative_to(out_dir.parent.parent)),
        },
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    print(json.dumps(manifest, indent=2))


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Path-1 labelled cohort for the slim DB comparison."
    )
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("analysis/databases/duckdb/cryptic-slim.duckdb"),
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=Path("analysis/results/cohort"),
    )
    args = parser.parse_args()
    build_cohort(args.db, args.out)


if __name__ == "__main__":
    main()
