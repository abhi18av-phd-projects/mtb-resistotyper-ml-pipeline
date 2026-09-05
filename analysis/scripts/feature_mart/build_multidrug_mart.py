"""Long-format multi-drug feature mart (design/15 Part B).

The build path for the multi-drug single stacked ensemble: ONE row per
(isolate, drug), `id__DRUG` a categorical feature, so a single H2O model can
score any (isolate, drug) and — the research bet — exploit shared genetic
background to lift low-prevalence drugs (design/03 "single-task" revisit
trigger; the transfer experiment in design/15's staged strategy).

Structure (per design/15 Part B):
  id__UNIQUEID, id__DRUG          identifiers; id__DRUG is ALSO a feature
  y__binary                       R/S for that (isolate, drug)
  cov__lineage_* / coverage /     per-isolate covariates (broadcast across
    country                        the isolate's drug-rows)
  raw__<gene>_<mutation>          the shared genomic feature universe =
                                   UNION of the per-drug top-K-by-correlation
                                   candidate pools, so every drug's own
                                   determinants (rpoB, katG, gyrA, rrs, …)
                                   are present; id__DRUG gates their relevance
  __group__                       = id__UNIQUEID (the GroupKFold key)
  __fold__                        StratifiedGroupKFold-by-isolate: all of an
                                   isolate's drug-rows share a fold, so random
                                   leakage across drug-rows is impossible
                                   (NON-NEGOTIABLE for an honest multi-drug AUC)

Why GroupKFold-by-isolate: the same isolate appears in up to 12 drug-rows;
plain k-fold would put its rows in both train and test. Grouping by isolate
(then layering the lineage-aware regime in evaluation) is the multi-drug
analogue of the per-drug shared-fold protocol.

Feature universe: union of per-drug top-K-by-|correlation| mutation columns
(the G3 prefilter, no full concordance needed — a mutation relevant for ANY
drug is present). This is the de-novo multi-drug feature set; a concordant
variant (union of per-drug concordant sets) follows once all 12 drugs have
concordance runs.

Run (causal venv: duckdb + pandas + sklearn):
    python -m analysis.scripts.feature_mart.build_multidrug_mart \
        --db analysis/databases/duckdb/cryptic-slim.duckdb \
        --marts analysis/results/feature_mart \
        --out analysis/results/feature_mart_multidrug \
        --per-drug-top 300
"""

from __future__ import annotations

import argparse
import glob
import json
import time
from pathlib import Path

import duckdb
import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedGroupKFold

CATALOGUE_NAME = "WHO-UCN-GTB-PCI-2023.5"
RANDOM_STATE = 42
N_FOLDS = 5
MIN_CARRIER = 10
DRUGS = ["RIF", "INH", "EMB", "MXF", "LEV", "ETH", "KAN", "AMI", "CFZ", "LZD", "BDQ", "DLM"]

# Mechanism class per drug (T14). The transfer bet: a single SE with id__DRUG +
# cov__drug_class lets same-class drugs share signal — FQ (MXF/LEV) and
# aminoglycoside (KAN/AMI) most obviously, and the efflux cross-resistance link
# BDQ↔CFZ↔(CLZ) is carried by the shared Rv0678 feature in the universe.
DRUG_CLASS = {
    "RIF": "rifamycin", "INH": "isoniazid", "EMB": "ethambutol",
    "MXF": "fluoroquinolone", "LEV": "fluoroquinolone",
    "ETH": "thioamide", "KAN": "aminoglycoside", "AMI": "aminoglycoside",
    "CFZ": "riminophenazine", "BDQ": "diarylquinoline",
    "LZD": "oxazolidinone", "DLM": "nitroimidazole",
}


def _feature_universe(marts_dir: Path, per_drug_top: int, mart_version: str = "v1.0.2") -> list[str]:
    """Union of per-drug top-K mutation columns by |corr with y__binary|.

    Uses the per-drug marts' already-carrier-windowed raw__ columns; ranks
    each drug's by label correlation (G3) and unions the top-K. Every drug's
    own resistance determinants land in the universe, so the multi-drug model
    is not handicapped on any single drug.
    """
    universe: set[str] = set()
    for drug in DRUGS:
        matches = glob.glob(str(marts_dir / f"feature_mart_{drug}_*{mart_version}.parquet"))
        if not matches:
            continue
        df = pd.read_parquet(matches[0])
        raw = [
            c for c in df.columns
            if c.startswith("raw__") and not c.startswith("raw__gene_")
            and c not in ("raw__n_mutations", "raw__mean_frs", "raw__mean_coverage", "raw__n_minor")
        ]
        carriers = df[raw].sum()
        eligible = carriers[carriers >= MIN_CARRIER].index
        corr = df[eligible].corrwith(df["y__binary"]).abs().sort_values(ascending=False)
        universe |= set(corr.head(per_drug_top).index)
    return sorted(universe)


def _parse_gene_mut(col: str) -> tuple[str, str]:
    body = col.removeprefix("raw__")
    gene, _, mut = body.partition("_")
    return gene, mut


def build(db_path: Path, marts_dir: Path, out_dir: Path, per_drug_top: int,
          mart_version: str = "v1.0.2") -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    con = duckdb.connect(str(db_path), read_only=True)
    drug_class_case = ("CASE DRUG "
                       + " ".join(f"WHEN '{d}' THEN '{c}'" for d, c in DRUG_CLASS.items())
                       + " ELSE 'other' END")

    # ---- long-format cohort: one row per (isolate, drug) --------------------
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE cohort AS
        SELECT ph.UNIQUEID, ph.DRUG,
               CASE WHEN ph.BINARY_PHENOTYPE = 'R' THEN 1 ELSE 0 END AS y__binary,
               g.LINEAGE, g.SUBLINEAGE, g.TB_DEPTH, g.TB_COVERAGE,
               ws.country, ws.dataset
        FROM   ukmyc_phenotypes ph
        JOIN   predictions  pr USING (UNIQUEID, DRUG)
        JOIN   genomes      g  USING (UNIQUEID)
        LEFT JOIN wgs_samples ws USING (UNIQUEID)
        WHERE  ph.BINARY_PHENOTYPE IN ('R','S')
          AND  pr.PREDICTION       IN ('R','S')
          AND  pr.CATALOGUE_NAME = '{CATALOGUE_NAME}'
    """)

    # ---- feature universe: union of per-drug top-K candidate pools ----------
    universe = _feature_universe(marts_dir, per_drug_top, mart_version)
    gene_muts = pd.DataFrame([_parse_gene_mut(c) for c in universe], columns=["GENE", "MUTATION"])
    gene_muts["feature_name"] = universe
    con.register("universe", gene_muts)

    # ---- pivot the universe mutations per isolate (one row per isolate) -----
    con.execute("""
        CREATE OR REPLACE TEMP TABLE mut_long AS
        SELECT DISTINCT m.UNIQUEID, u.feature_name, 1 AS present
        FROM   mutations m
        JOIN   universe u ON u.GENE = m.GENE AND u.MUTATION = m.MUTATION
        WHERE  m.UNIQUEID IN (SELECT DISTINCT UNIQUEID FROM cohort)
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE feat_wide AS
        PIVOT mut_long ON feature_name USING MAX(present) GROUP BY UNIQUEID
    """)

    # ---- assemble long-format: cohort core + covariates + broadcast features -
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE core AS
        SELECT
            UNIQUEID                                        AS id__UNIQUEID,
            DRUG                                            AS id__DRUG,
            {drug_class_case}                               AS cov__drug_class,
            y__binary,
            CASE WHEN LINEAGE='lineage1' THEN 1 ELSE 0 END  AS cov__lineage_L1,
            CASE WHEN LINEAGE='lineage2' THEN 1 ELSE 0 END  AS cov__lineage_L2,
            CASE WHEN LINEAGE='lineage3' THEN 1 ELSE 0 END  AS cov__lineage_L3,
            CASE WHEN LINEAGE='lineage4' THEN 1 ELSE 0 END  AS cov__lineage_L4,
            COALESCE(LINEAGE,'unknown')                     AS cov__lineage_raw,
            COALESCE(country,'')                            AS cov__country,
            TB_DEPTH::FLOAT                                 AS cov__median_coverage,
            TB_COVERAGE::FLOAT                              AS cov__tb_breadth
        FROM cohort
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE mart AS
        SELECT core.*, feat_wide.* EXCLUDE (UNIQUEID)
        FROM core
        LEFT JOIN feat_wide ON feat_wide.UNIQUEID = core.id__UNIQUEID
    """)

    df = con.execute("SELECT * FROM mart").fetchdf()
    # NULL pivot cells (isolate carries none of that mutation) -> 0
    feat_cols = [c for c in df.columns if c.startswith("raw__")]
    df[feat_cols] = df[feat_cols].fillna(0).astype("int8")

    # ---- StratifiedGroupKFold by isolate ------------------------------------
    df["__group__"] = df["id__UNIQUEID"]
    sgkf = StratifiedGroupKFold(n_splits=N_FOLDS, shuffle=True, random_state=RANDOM_STATE)
    fold = np.zeros(len(df), dtype="int8")
    for f, (_, val) in enumerate(sgkf.split(df, df["y__binary"], groups=df["__group__"])):
        fold[val] = f
    df["__fold__"] = fold

    # sanity: no isolate split across folds
    leak = df.groupby("id__UNIQUEID")["__fold__"].nunique()
    assert (leak == 1).all(), f"{(leak > 1).sum()} isolates leaked across folds"

    mart_path = out_dir / "feature_mart_MULTIDRUG_cryptic-slim-2026.05_v1.0.0.parquet"
    df.to_parquet(mart_path, index=False, compression="zstd")

    manifest = {
        "mart": "multidrug-longformat",
        "mart_version": "v1.0.0",
        "cryptic_version": "slim-2026.05",
        "catalogue": CATALOGUE_NAME,
        "n_rows": len(df),
        "n_isolates": int(df["id__UNIQUEID"].nunique()),
        "n_drugs": int(df["id__DRUG"].nunique()),
        "drugs": sorted(df["id__DRUG"].unique().tolist()),
        "feature_universe_size": len(feat_cols),
        "feature_universe_rule": f"union of per-drug top-{per_drug_top} by |corr with label| (G3)",
        "rows_per_drug": df["id__DRUG"].value_counts().to_dict(),
        "pct_R_per_drug": (df.groupby("id__DRUG")["y__binary"].mean() * 100).round(2).to_dict(),
        "fold_assignment": {
            "method": "StratifiedGroupKFold",
            "n_splits": N_FOLDS,
            "group": "id__UNIQUEID",
            "shuffle": True,
            "random_state": RANDOM_STATE,
            "isolate_leakage_across_folds": int((leak > 1).sum()),
        },
        "wall_time_seconds": round(time.time() - t0, 2),
        "outputs": {"parquet": mart_path.name},
    }
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def main() -> None:
    p = argparse.ArgumentParser(description="Long-format multi-drug feature mart (design/15 Part B).")
    p.add_argument("--db", type=Path, default=Path("analysis/databases/duckdb/cryptic-slim.duckdb"))
    p.add_argument("--marts", type=Path, default=Path("analysis/results/feature_mart"))
    p.add_argument("--out", type=Path, default=Path("analysis/results/feature_mart_multidrug"))
    p.add_argument("--per-drug-top", type=int, default=300)
    p.add_argument("--mart-version", default="v1.0.2", help="per-drug mart version glob (e.g. vFULL for the full dataset)")
    args = p.parse_args()

    m = build(args.db, args.marts, args.out, args.per_drug_top, args.mart_version)
    print(json.dumps({k: m[k] for k in
                      ["n_rows", "n_isolates", "n_drugs", "feature_universe_size", "wall_time_seconds"]}, indent=2))


if __name__ == "__main__":
    main()
