"""Per-drug Path-1 feature mart builder (pass 1: no mdl__/causal__/cf__).

Realises §9 step 3 of the 2026-06-30 slim-DB catalogue-comparison
handoff, per the design/09 column-naming convention.

Sections built in pass 1:
    id__, y__, raw__, cat__, cov__, qf__, __fold__

Sections deferred to pass 2 (per user decision):
    mdl__ (EBM/LR OOF stacking columns)
    causal__ (lineage residuals, CATE, refutation, propensity, ICP, mediation)
    cf__ (DiCE counterfactuals)

Leakage boundary (per handoff §4 + locked user decision):
    cat__who_predict, cat__who_p, cat__catalogue_disagrees_with_label
        -> NEVER features in either arm
    cat__n_per_mut_R / _S / _U / _F
        -> CATALOGUE-INFORMED arm only (structural per-mutation counts)
        -> NOTE: design/09 spec called these `cat__n_tier1/tier2/interim`
           but the slim DB's effects table is trimmed (no EVIDENCE /
           PREDICTION_VALUES dicts), so tier grading isn't recoverable.
           Per-mutation R/S/U/F counts are the closest honest substitute.
    No cat__ in the de-novo arm (downstream filter:
        `[c for c in cols if not c.startswith('cat__')]`)

Run for one drug:
    python -m analysis.scripts.feature_mart.build_mart \
        --drug RIF \
        --db analysis/databases/duckdb/cryptic-slim.duckdb \
        --out analysis/results/feature_mart \
        --version v1.0.0

Run for all 12 drugs (sequential, ~1-3 min each on a laptop):
    for drug in RIF INH EMB MXF LEV ETH KAN AMI CFZ LZD BDQ DLM; do
      python -m analysis.scripts.feature_mart.build_mart --drug $drug
    done
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import random
import time
from pathlib import Path

import duckdb

try:
    from analysis.scripts.feature_mart import fe_identity
except ImportError:  # run as a plain script from this directory
    import fe_identity


CATALOGUE_NAME = "WHO-UCN-GTB-PCI-2023.5"
TOP_N_MUTATIONS = 1000

# Canonical TB determinant genes per drug — their carrier-windowed mutations are
# UNIONed into the top-N pool so the true driver is never crowded out of the
# prevalence ranking. Critical on the FULL dataset: with PE/PPE + indels present,
# the naive top-1000 drops rpoB_S450L for RIF (prevalence-rank swamped by extra
# lineage/indel mutations) — see design/16 full-vs-slim experiment. Mirrors
# causal.py's DETERMINANT_GENES (kept local to avoid importing its heavy deps).
DETERMINANT_GENES: dict[str, set[str]] = {
    "RIF": {"rpoB"},
    "INH": {"katG", "inhA", "fabG1", "ahpC"},
    "EMB": {"embB", "embA", "embC", "embR"},
    "PZA": {"pncA", "rpsA", "panD"},
    "MXF": {"gyrA", "gyrB"}, "LEV": {"gyrA", "gyrB"},
    "AMI": {"rrs", "eis"}, "KAN": {"rrs", "eis"}, "STM": {"rpsL", "rrs", "gid"},
    "ETH": {"ethA", "ethR", "inhA", "fabG1"},
    "LZD": {"rrl", "rplC"},
    "BDQ": {"Rv0678", "atpE", "pepQ", "mmpL5", "mmpS5"},
    "CFZ": {"Rv0678", "pepQ", "mmpL5", "mmpS5"},
    # Nitroimidazole activation pathway (DLM/PMD) — 6 genes per in-vitro selection
    # spectrum (Int J Antimicrob Agents 2026, 10.1016/j.ijantimicag.2026.107789).
    # Rv0678 dropped: efflux (BDQ/CFZ) confounder, appears in DLM-R only via MDR
    # co-occurrence. See design/19 §8.
    "DLM": {"ddn", "fgd1", "fbiA", "fbiB", "fbiC", "fbiD"},
}
N_FOLDS = 5
RANDOM_STATE = 42
MART_VERSION_DEFAULT = "v1.0.2"

# Co-resistance confounder panel (added v1.0.2, design/15 finding). MDR
# co-selection means an isolate resistant to one drug is disproportionately
# resistant to its MDR partner (RIF<->INH especially) — so a drug's genetic
# marker (katG/S315T for INH) spuriously "predicts" another drug (RIF) via
# strain background, not causation. Exposing each isolate's DST status for
# the other first-line drugs as cov__ confounders lets the causal layer
# condition on the MDR background. The current drug is always excluded from
# its own mart's co-resistance columns (else it would leak the label).
CORESIST_PANEL = ("RIF", "INH", "EMB")

# Mutation-selection filter (deviation from design/09 default).
#
# `design/09` defaults to "top-N by sample-prevalence". Run against the
# slim DB this picks up only near-universal lineage markers (variants in
# >99% of isolates carry zero discriminative signal — they're handled by
# cov__lineage_*) and EXCLUDES the canonical resistance mutations like
# rpoB/S450L (which sits at carrier-rate ~25% per the RIF cohort).
#
# Fix: apply a sensible carrier-rate window BEFORE the top-N cut:
#   - drop carrier-rate > MAX_CARRIER_FRAC (lineage / near-universal noise)
#   - drop carrier-rate < MIN_CARRIER_COUNT (singletons / rare-call noise)
# Empirically this lifts the canonical resistance mutations into the
# top-N feature set without growing N. Documented in the mart's metadata
# sidecar so consumers know the deviation from design/09.
MAX_CARRIER_FRAC = 0.95
MIN_CARRIER_COUNT = 10



def _read_db_metadata(con) -> dict:
    """Everything the database records about its own build.

    Empty when the table is absent, which is itself worth recording: a mart built
    from a database that carries no provenance should say so rather than inherit
    a plausible-looking default.
    """
    try:
        rows = con.execute("SELECT key, value FROM _database_metadata").fetchall()
    except Exception:
        return {}
    return {k: v for k, v in rows}


_RELEASE_RE = re.compile(r"cryptic-tables-(v[0-9]+(?:\.[0-9]+)*)")


def _genome_columns(con) -> set[str]:
    """Columns the genomes table actually has, upper-cased.

    Not every CRyPTIC build carries the same optional covariates: v3.4.0 has
    SUBLINEAGE, v2.1.2 does not. Selecting it unconditionally is a binder error
    that stops the mart before any drug is built, and hard-coding it out would
    silently drop a real covariate from the release that does have it. So the
    column is selected when present and NULL-substituted when absent, and the
    mart metadata records which happened -- two releases whose feature spaces
    differ must say so, or a comparison between their models is not a
    comparison of the data.
    """
    return _table_columns(con, "genomes")


def _table_columns(con, table: str) -> set[str]:
    """Columns a table actually has, upper-cased."""
    rows = con.execute(
        "SELECT column_name FROM information_schema.columns "
        "WHERE lower(table_name) = ?", [table.lower()]).fetchall()
    return {str(c).upper() for (c,) in rows}


def _read_release(con) -> str:
    """The CRyPTIC release this database was built from, plus its scope.

    Refuses to guess. A wrong release label silently mis-attributes every number
    downstream, and the failure is invisible because the mart still builds.

    The label follows CRyPTIC's own convention -- ``v3.4.0``, ``v2.1.2`` -- and
    never a build date of ours. ``slim-2026.05`` was such a date, and it names
    nothing a reader of the compendium could look up: the slim database is built
    from cryptic-tables-v3.4.0, the same release as the full one. What differs
    is scope, so scope is what the suffix carries: a database whose metadata
    records the slim pruning is ``v3.4.0-slim``, and the unpruned build of the
    same release is ``v3.4.0``. Two sets that differ only in scope must still be
    distinguishable, or a comparison between them cannot be stated.

    The older slim build predates ``source_version`` and records the release
    only inside ``source_dir``; that path is parsed rather than rejected,
    because refusing it would make the one shipped bundle set unnameable.
    """
    meta = _read_db_metadata(con)
    # source_version is not always the version. Several builds recorded the
    # source DIRECTORY BASENAME there -- "src" in the full builds, and
    # "cryptic-tables-v3.4.0" in the slim one -- so the field is normalised
    # rather than trusted: whatever carries a cryptic-tables-<version> is
    # reduced to the version, and a bare basename falls through to source_dir.
    v = (meta.get("source_version") or "").strip()
    m = _RELEASE_RE.search(v)
    if m:
        v = m.group(1)
    elif not v or v in {"src", "source", "."} or not v.startswith("v"):
        m = _RELEASE_RE.search(meta.get("source_dir") or "")
        v = m.group(1) if m else ""
    if not v:
        raise SystemExit(
            "the database records no usable source_version "
            "and no release in its source_dir. Rebuild with MTB_CRYPTIC_VERSION "
            "set, or backfill _database_metadata, before building a mart: a mart "
            "that cannot name its release cannot be compared with one from "
            "another release."
        )
    if any(k.startswith("slim_") for k in meta):
        return f"{v}-slim"
    return v


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def _stratified_fold(y_binary: list[int], n_splits: int, random_state: int) -> list[int]:
    """Deterministic StratifiedKFold equivalent using only stdlib.

    Within each class, shuffles indices with the seeded RNG, then deals
    round-robin into folds. Matches sklearn's StratifiedKFold(shuffle=True)
    behaviour to within RNG implementation.
    """
    rng = random.Random(random_state)
    by_class: dict[int, list[int]] = {0: [], 1: []}
    for i, y in enumerate(y_binary):
        by_class[int(y)].append(i)
    fold = [0] * len(y_binary)
    for cls, idxs in by_class.items():
        rng.shuffle(idxs)
        for k, i in enumerate(idxs):
            fold[i] = k % n_splits
    return fold


def build_mart(
    drug: str,
    db_path: Path,
    out_dir: Path,
    *,
    mart_version: str = MART_VERSION_DEFAULT,
    top_n: int = TOP_N_MUTATIONS,
    n_folds: int = N_FOLDS,
    random_state: int = RANDOM_STATE,
    # --- feature-engineering knobs -------------------------------------------
    # Previously module constants, so changing an FE choice meant editing this
    # file. They are arguments now because the mart is the reproducibility unit:
    # every one of these lands in the provenance sidecar, so a mart states how it
    # was built and two marts are diffable.
    mutation_selector: str = "top_n",
    max_carrier_frac: float = MAX_CARRIER_FRAC,
    min_carrier_count: int = MIN_CARRIER_COUNT,
    gene_level_flags: bool = True,
    fold_strategy: str = "stratified",
    fe_name: str = "default",
) -> dict:
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    con = duckdb.connect(str(db_path), read_only=True)
    # The optional genome covariates, and the type each is read back as. The
    # two CRyPTIC builds share almost nothing here: v3.4.0's genomes table
    # carries SUBLINEAGE, TB_DEPTH and TB_COVERAGE; v2.1.2's carries none of
    # them. Missing ones become NULL of the right type rather than a binder
    # error, and every one that was missing is named in the mart metadata, so a
    # model trained on a release without coverage covariates is identifiable as
    # such instead of looking like one whose coverage happened to be unknown.
    _present = _genome_columns(con)
    _optional = {"SUBLINEAGE": "VARCHAR", "TB_DEPTH": "DOUBLE",
                 "TB_COVERAGE": "DOUBLE"}
    genome_optional_select = ", ".join(
        f"g.{col}" if col in _present else f"CAST(NULL AS {typ}) AS {col}"
        for col, typ in _optional.items())
    missing_covariates = sorted(c for c in _optional if c not in _present)

    # wgs_samples diverges the same way: v3.4.0 publishes COUNTRY and DATASET,
    # v2.1.2 publishes neither. They feed cov__country and cov__dataset_version,
    # which are cohort descriptors rather than genotype, so a NULL is an honest
    # "not published by this release" and not a silently wrong value.
    _wgs = _table_columns(con, "wgs_samples")
    _wgs_optional = {"COUNTRY": "VARCHAR", "DATASET": "VARCHAR"}
    wgs_optional_select = ", ".join(
        f"ws.{col}" if col in _wgs else f"CAST(NULL AS {typ}) AS {col}"
        for col, typ in _wgs_optional.items())
    missing_covariates += sorted(c for c in _wgs_optional if c not in _wgs)

    # The two NUMERIC covariates are OMITTED when their source column is absent,
    # not emitted as NULL. A column that is NULL for every row carries no
    # information, and it is not inert downstream: causal.py hands the covariate
    # block to econml, which rejects it outright with "Input contains NaN" and
    # takes the whole run with it. The string covariates above COALESCE to '' and
    # are harmless either way. Consumers already select covariates by
    # intersection (`[c for c in COV if c in df.columns]`), so an absent column
    # is the case they are written for; an all-NaN one is not.
    _quality_cov = [("TB_DEPTH", "cov__median_coverage"),
                    ("TB_COVERAGE", "cov__tb_breadth")]
    quality_covariate_select = "".join(
        f"{src}::FLOAT AS {name},\n            "
        for src, name in _quality_cov if src in _present)

    # ---- cohort ------------------------------------------------------------
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE cohort AS
        SELECT ph.UNIQUEID, ph.DRUG,
               ph.BINARY_PHENOTYPE AS truth,
               pr.PREDICTION       AS catalogue,
               g.LINEAGE, {genome_optional_select},
               {wgs_optional_select},
               ph.PHENOTYPE_QUALITY
        FROM   ukmyc_phenotypes ph
        JOIN   predictions  pr USING (UNIQUEID, DRUG)
        JOIN   genomes      g  USING (UNIQUEID)
        LEFT JOIN wgs_samples ws USING (UNIQUEID)
        WHERE  ph.DRUG = ?
          AND  ph.BINARY_PHENOTYPE IN ('R','S')
          AND  pr.PREDICTION       IN ('R','S')
          AND  pr.CATALOGUE_NAME = ?
    """, [drug, CATALOGUE_NAME])

    cohort_size = con.execute("SELECT COUNT(*) FROM cohort").fetchone()[0]
    if cohort_size == 0:
        raise SystemExit(f"empty cohort for drug={drug} — check the slim DB")

    # ---- raw__ section: top-N mutations per drug, pivoted ------------------
    # Pre-filter the candidate pool by carrier-rate window (see
    # MAX_CARRIER_FRAC / MIN_CARRIER_COUNT comments) so near-universal
    # lineage markers don't crowd out resistance variants.
    cohort_n = con.execute("SELECT COUNT(*) FROM cohort").fetchone()[0]
    max_carriers = int(cohort_n * max_carrier_frac)
    # Determinant-gene seed clause: carrier-windowed mutations in the drug's
    # determinant genes are UNIONed in regardless of prevalence rank.
    det_genes = DETERMINANT_GENES.get(drug, set())
    seed_filter = (
        "GENE IN (" + ", ".join(f"'{g}'" for g in sorted(det_genes)) + ")"
        if det_genes else "FALSE"
    )
    # The ranking statistic is the FE choice. Prevalence (top_n) is cheap and
    # label-blind; chi2 and mutual information rank by association with the
    # phenotype instead, which surfaces rarer variants that actually discriminate
    # but risks selecting on the label, so downstream selection stays fold-internal.
    #
    # Counts for the 2x2 table (carried x resistant), computed once and reused:
    #   a = carried & R    b = carried & S
    #   c = absent  & R    d = absent  & S
    # LABEL BOUNDARY. Everything build_mart does is label-free by construction, so
    # the mart it produces is valid for every fold. Ranking by association with the
    # phenotype (chi2, mutual information) is NOT label-free: it would choose the
    # column universe while looking at the held-out lineage, which is the same class
    # of selection leakage the fold-internal concordance step exists to remove --
    # milder, because it only decides which columns exist, but the same error.
    # Those rankers live downstream of the boundary, in fe_steps/rank_association.
    LABEL_USING_RANKERS = {"chi2", "mi"}
    if mutation_selector in LABEL_USING_RANKERS:
        raise ValueError(
            f"mutation_selector={mutation_selector!r} ranks by association with the "
            "phenotype, so it cannot run at mart-build time: the mart is built once "
            "and shared by every fold, and the ranking would see the held-out lineage. "
            "Build the mart with --mutation-selector all, then apply the ranking per "
            "fold with the rank_association FE step."
        )

    RANKERS = {
        "top_n": "n_isolates DESC",
        "all": "n_isolates DESC",
    }
    if mutation_selector not in RANKERS:
        raise ValueError(f"unknown mutation_selector {mutation_selector!r}")
    order_by = RANKERS[mutation_selector]
    limit_clause = "" if mutation_selector == "all" else f"LIMIT {top_n}"

    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE top_muts AS
        WITH carriers AS (
            SELECT m.GENE, m.MUTATION, m.UNIQUEID
            FROM   mutations m
            WHERE  m.UNIQUEID IN (SELECT UNIQUEID FROM cohort)
            GROUP BY m.GENE, m.MUTATION, m.UNIQUEID
        ),
        eligible AS (
            SELECT c.GENE, c.MUTATION, COUNT(*) AS n_isolates
            FROM   carriers c JOIN cohort co USING (UNIQUEID)
            GROUP BY c.GENE, c.MUTATION
            HAVING COUNT(*) BETWEEN {min_carrier_count} AND {max_carriers}
        )
        SELECT GENE, MUTATION, n_isolates FROM (
            SELECT GENE, MUTATION, n_isolates FROM eligible
            ORDER BY {order_by}, GENE, MUTATION
            {limit_clause}
        )
        UNION
        SELECT GENE, MUTATION, n_isolates FROM eligible WHERE {seed_filter}
    """)


    # Per-(UID, gene, mutation) presence — narrow EAV table, then PIVOT
    con.execute("""
        CREATE OR REPLACE TEMP TABLE mut_long AS
        SELECT DISTINCT
               m.UNIQUEID,
               'raw__' || m.GENE || '_' || m.MUTATION AS feature_name,
               1 AS present
        FROM   mutations m
        JOIN   top_muts  t USING (GENE, MUTATION)
        WHERE  m.UNIQUEID IN (SELECT UNIQUEID FROM cohort)
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE raw_wide AS
        PIVOT mut_long
        ON feature_name
        USING MAX(present)
        GROUP BY UNIQUEID
    """)

    # Gene-level any-mut flags (one per gene that appears in top_muts)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE gene_long AS
        SELECT DISTINCT
               m.UNIQUEID,
               'raw__gene_' || m.GENE AS feature_name,
               1 AS present
        FROM   mutations m
        WHERE  m.UNIQUEID IN (SELECT UNIQUEID FROM cohort)
          AND  m.GENE IN (SELECT DISTINCT GENE FROM top_muts)
    """)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE gene_wide AS
        PIVOT gene_long
        ON feature_name
        USING MAX(present)
        GROUP BY UNIQUEID
    """)

    # The minor-allele flag is called IS_MINOR in v3.4.0 and IS_MINOR_ALLELE in
    # v2.1.2 -- the same measurement under two names, so it is mapped rather than
    # lost. COVERAGE has no v2.1.2 equivalent: that release never published it,
    # which is a fact about CRyPTIC and not about this build, so raw__mean_coverage
    # becomes NULL and is named in missing_mutation_columns. Both are recorded, so
    # a cross-release comparison can state which features it could not hold equal.
    _mut = _table_columns(con, "mutations")
    _minor = next((c for c in ("IS_MINOR", "IS_MINOR_ALLELE") if c in _mut), None)
    n_minor_select = (f"SUM(CASE WHEN {_minor} THEN 1 ELSE 0 END)::INTEGER"
                      if _minor else "CAST(NULL AS INTEGER)")
    mean_coverage_select = ("AVG(COVERAGE)::FLOAT" if "COVERAGE" in _mut
                            else "CAST(NULL AS FLOAT)")
    missing_mutation_columns = ([] if "COVERAGE" in _mut else ["COVERAGE"]) + \
                               ([] if _minor else ["IS_MINOR"])

    # Summary columns (n_mutations, mean_frs, mean_coverage) — over ALL
    # mutations carried by the sample, not just the top-N. Honest signal.
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE raw_summary AS
        SELECT UNIQUEID,
               COUNT(*)::INTEGER AS raw__n_mutations,
               AVG(FRS)::FLOAT   AS raw__mean_frs,
               {mean_coverage_select} AS raw__mean_coverage,
               {n_minor_select} AS raw__n_minor
        FROM   mutations
        WHERE  UNIQUEID IN (SELECT UNIQUEID FROM cohort)
        GROUP BY UNIQUEID
    """)

    # ---- cat__ section: per-mutation catalogue R/S/U/F counts --------------
    # (Catalogue-informed arm only — downstream filter drops `cat__` for de-novo.
    #  cat__who_predict / cat__who_p / cat__catalogue_disagrees_with_label
    #  are deliberately NOT computed: handoff §4 says the catalogue
    #  prediction itself is never a feature.)
    con.execute("""
        CREATE OR REPLACE TEMP TABLE cat_counts AS
        SELECT e.UNIQUEID,
               SUM(CASE WHEN e.PREDICTION='R' THEN 1 ELSE 0 END)::INTEGER AS cat__n_per_mut_R,
               SUM(CASE WHEN e.PREDICTION='S' THEN 1 ELSE 0 END)::INTEGER AS cat__n_per_mut_S,
               SUM(CASE WHEN e.PREDICTION='U' THEN 1 ELSE 0 END)::INTEGER AS cat__n_per_mut_U,
               SUM(CASE WHEN e.PREDICTION='F' THEN 1 ELSE 0 END)::INTEGER AS cat__n_per_mut_F
        FROM   effects e
        JOIN   cohort  c USING (UNIQUEID, DRUG)
        WHERE  e.CATALOGUE_NAME = ?
        GROUP BY e.UNIQUEID
    """, [CATALOGUE_NAME])

    # ---- co-resistance confounders (cov__coresist_*) -----------------------
    # Per-isolate DST R/S for the other first-line drugs (current drug
    # excluded to avoid leaking the label). NULL when the isolate wasn't
    # tested for that drug — filled to 0 downstream (treat untested as the
    # non-resistant baseline; a coverage flag could refine this later).
    coresist_drugs = [d for d in CORESIST_PANEL if d != drug]
    coresist_cols_sql = ",\n               ".join(
        f"MAX(CASE WHEN ph.DRUG = '{d}' AND ph.BINARY_PHENOTYPE = 'R' THEN 1 "
        f"WHEN ph.DRUG = '{d}' AND ph.BINARY_PHENOTYPE = 'S' THEN 0 END)::TINYINT "
        f'AS "cov__coresist_{d}"'
        for d in coresist_drugs
    )
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE coresist AS
        SELECT ph.UNIQUEID,
               {coresist_cols_sql},
               SUM(CASE WHEN ph.DRUG <> '{drug}' AND ph.BINARY_PHENOTYPE = 'R'
                        THEN 1 ELSE 0 END)::SMALLINT AS "cov__coresist_n_other_R"
        FROM   ukmyc_phenotypes ph
        WHERE  ph.UNIQUEID IN (SELECT UNIQUEID FROM cohort)
          AND  ph.BINARY_PHENOTYPE IN ('R','S')
        GROUP BY ph.UNIQUEID
    """)

    # ---- cov__ + qf__ + y__ + id__ section ---------------------------------
    # f-string: the quality covariates are interpolated, because which of them
    # exist depends on the release (see quality_covariate_select above).
    con.execute(f"""
        CREATE OR REPLACE TEMP TABLE core AS
        SELECT
            UNIQUEID                                       AS id__UNIQUEID,
            DRUG                                           AS id__DRUG,
            ?  AS id__catalogue_name,
            -- y__
            CASE WHEN truth = 'R' THEN 1 ELSE 0 END        AS y__binary,
            truth                                          AS y__truth_raw,
            PHENOTYPE_QUALITY                              AS y__quality,
            mic.LOG2MIC                                    AS y__log2mic,
            -- cov__
            CASE WHEN LINEAGE = 'lineage1' THEN 1 ELSE 0 END AS cov__lineage_L1,
            CASE WHEN LINEAGE = 'lineage2' THEN 1 ELSE 0 END AS cov__lineage_L2,
            CASE WHEN LINEAGE = 'lineage3' THEN 1 ELSE 0 END AS cov__lineage_L3,
            CASE WHEN LINEAGE = 'lineage4' THEN 1 ELSE 0 END AS cov__lineage_L4,
            COALESCE(LINEAGE, 'unknown')                   AS cov__lineage_raw,
            COALESCE(SUBLINEAGE, '')                       AS cov__sublineage,
            COALESCE(country, '')                          AS cov__country,
            COALESCE(dataset, '')                          AS cov__dataset_version,
            {quality_covariate_select}
            -- qf__
            (PHENOTYPE_QUALITY = 'HIGH')                   AS qf__phenotype_quality_high,
            (PHENOTYPE_QUALITY IN ('HIGH','MEDIUM'))       AS qf__phenotype_quality_med_or_high,
            (LINEAGE IN ('lineage1','lineage2','lineage3','lineage4')
              AND TB_DEPTH > 30
              AND PHENOTYPE_QUALITY IN ('HIGH','MEDIUM'))  AS qf__train_eligible
        FROM cohort
        -- One isolate can be plated more than once; take the median LOG2MIC,
        -- matching evaluate_mic_regression so the graded phenotype means the
        -- same thing wherever it is read.
        LEFT JOIN (
            SELECT UNIQUEID, DRUG, median(LOG2MIC) AS LOG2MIC
            FROM   ukmyc_phenotypes
            WHERE  LOG2MIC IS NOT NULL
            GROUP BY UNIQUEID, DRUG
        ) mic USING (UNIQUEID, DRUG)
    """, [CATALOGUE_NAME])

    # ---- final wide assembly ----------------------------------------------
    # Left-join everything onto core. Missing raw__ / gene_wide / cat / summary
    # values get filled with 0 / NULL by COALESCE downstream.
    con.execute("""
        CREATE OR REPLACE TEMP TABLE mart AS
        SELECT
            core.*,
            COALESCE(raw_summary.raw__n_mutations,   0) AS raw__n_mutations,
            COALESCE(raw_summary.raw__mean_frs,    NULL) AS raw__mean_frs,
            COALESCE(raw_summary.raw__mean_coverage, NULL) AS raw__mean_coverage,
            COALESCE(raw_summary.raw__n_minor,       0) AS raw__n_minor,
            COALESCE(cat_counts.cat__n_per_mut_R, 0) AS cat__n_per_mut_R,
            COALESCE(cat_counts.cat__n_per_mut_S, 0) AS cat__n_per_mut_S,
            COALESCE(cat_counts.cat__n_per_mut_U, 0) AS cat__n_per_mut_U,
            COALESCE(cat_counts.cat__n_per_mut_F, 0) AS cat__n_per_mut_F,
            coresist.* EXCLUDE (UNIQUEID),
            raw_wide.*  EXCLUDE (UNIQUEID),
            gene_wide.* EXCLUDE (UNIQUEID)
        FROM core
        LEFT JOIN raw_summary ON raw_summary.UNIQUEID = core.id__UNIQUEID
        LEFT JOIN cat_counts  ON cat_counts.UNIQUEID  = core.id__UNIQUEID
        LEFT JOIN coresist    ON coresist.UNIQUEID    = core.id__UNIQUEID
        LEFT JOIN raw_wide    ON raw_wide.UNIQUEID    = core.id__UNIQUEID
        LEFT JOIN gene_wide   ON gene_wide.UNIQUEID   = core.id__UNIQUEID
    """)

    # NULLs -> 0 in (a) pivot columns (no rows in mut_long for that
    # UID×feature) and (b) co-resistance columns (isolate untested for that
    # drug, or absent from the coresist join): treat untested as the
    # non-resistant baseline, consistent with the CORESIST_PANEL comment.
    raw_cols = [
        r[0] for r in con.execute(
            "SELECT column_name FROM (DESCRIBE mart) "
            "WHERE column_name LIKE 'raw\\_\\_%' ESCAPE '\\' "
            "  AND column_name NOT IN ('raw__n_mutations','raw__mean_frs','raw__mean_coverage','raw__n_minor')"
        ).fetchall()
    ]
    coresist_cols = [
        r[0] for r in con.execute(
            "SELECT column_name FROM (DESCRIBE mart) "
            "WHERE column_name LIKE 'cov\\_\\_coresist\\_%' ESCAPE '\\'"
        ).fetchall()
    ]
    fill_cols = raw_cols + coresist_cols
    if fill_cols:
        coalesce_clauses = ",\n        ".join(
            f'COALESCE("{c}", 0)::TINYINT AS "{c}"' for c in raw_cols
        )
        if coresist_cols:
            coalesce_clauses += ",\n        " + ",\n        ".join(
                f'COALESCE("{c}", 0)::SMALLINT AS "{c}"' for c in coresist_cols
            )
        other_cols = [
            r[0] for r in con.execute("SELECT column_name FROM (DESCRIBE mart)").fetchall()
            if r[0] not in fill_cols
        ]
        other_clause = ", ".join(f'"{c}"' for c in other_cols)
        con.execute(f"""
            CREATE OR REPLACE TEMP TABLE mart2 AS
            SELECT {other_clause},
                   {coalesce_clauses}
            FROM mart
        """)
        con.execute("DROP TABLE mart")
        con.execute("ALTER TABLE mart2 RENAME TO mart")

    # ---- __fold__ assignment (StratifiedKFold-equivalent) ------------------
    group_col = {"lineage": "cov__lineage_raw", "country": "cov__country"}.get(fold_strategy)
    select_cols = "id__UNIQUEID, y__binary" + (f', "{group_col}"' if group_col else "")
    y_rows = con.execute(
        f"SELECT {select_cols} FROM mart ORDER BY id__UNIQUEID"
    ).fetchall()
    uids = [r[0] for r in y_rows]
    y_bin = [r[1] for r in y_rows]
    group_key = [r[2] for r in y_rows] if group_col else None
    if fold_strategy == "stratified":
        folds = _stratified_fold(y_bin, n_splits=n_folds, random_state=random_state)
    elif fold_strategy in ("lineage", "country"):
        # Group folds. Grouping on lineage is what makes a downstream lineage-LOO
        # comparison like-for-like: without it a "cross-validated" number is the
        # random regime, which this project's own results show is inflated.
        groups = sorted({g for g in group_key if g and g != "unknown"})
        if not groups:
            raise ValueError(f"fold_strategy={fold_strategy!r} but {group_col} is empty")
        index = {g: i for i, g in enumerate(groups)}
        # One fold per group. An isolate whose group is missing goes to the last
        # fold rather than being dropped, so the mart's row count is unchanged by
        # the fold strategy and two strategies stay comparable.
        folds = [index.get(g, len(groups) - 1) for g in group_key]
    else:
        raise ValueError(f"unknown fold_strategy {fold_strategy!r}")
    con.execute("CREATE OR REPLACE TEMP TABLE fold_map(UNIQUEID VARCHAR, fold TINYINT)")
    con.executemany(
        "INSERT INTO fold_map VALUES (?, ?)",
        list(zip(uids, folds)),
    )
    con.execute("""
        CREATE OR REPLACE TEMP TABLE mart_final AS
        SELECT mart.*, fold_map.fold AS "__fold__"
        FROM   mart
        LEFT JOIN fold_map ON fold_map.UNIQUEID = mart.id__UNIQUEID
    """)

    # ---- write parquet -----------------------------------------------------
    # Read the release from the database rather than asserting it. This was
    # hardcoded to "slim-2026.05", so a mart built from the FULL v3.4.0 release
    # labelled itself as the slim one. Every mart in a database-comparison
    # campaign then carried the same stamp, which would have made the arms
    # indistinguishable in their own metadata: a provenance record that is
    # constant across the thing it is supposed to distinguish is worse than none,
    # because it looks like provenance.
    cryptic_version = _read_release(con)
    mart_filename = f"feature_mart_{drug}_cryptic-{cryptic_version}_{mart_version}.parquet"
    mart_path = out_dir / mart_filename
    con.execute(
        f"COPY mart_final TO '{mart_path}' (FORMAT 'parquet', COMPRESSION 'zstd')"
    )

    # ---- metadata sidecar --------------------------------------------------
    all_cols = [r[0] for r in con.execute("SELECT column_name FROM (DESCRIBE mart_final)").fetchall()]
    sections = {"id": [], "y": [], "raw": [], "cat": [], "cov": [], "qf": [], "internal": []}
    for c in all_cols:
        if c.startswith("id__"):       sections["id"].append(c)
        elif c.startswith("y__"):      sections["y"].append(c)
        elif c.startswith("raw__"):    sections["raw"].append(c)
        elif c.startswith("cat__"):    sections["cat"].append(c)
        elif c.startswith("cov__"):    sections["cov"].append(c)
        elif c.startswith("qf__"):     sections["qf"].append(c)
        elif c.startswith("__"):       sections["internal"].append(c)
    cat_FP = con.execute(
        "SELECT SUM(CASE WHEN core.y__binary=0 AND coh.catalogue='R' THEN 1 ELSE 0 END) "
        "FROM core LEFT JOIN cohort coh ON coh.UNIQUEID = core.id__UNIQUEID"
    ).fetchone()[0]
    cat_FN = con.execute(
        "SELECT SUM(CASE WHEN core.y__binary=1 AND coh.catalogue='S' THEN 1 ELSE 0 END) "
        "FROM core LEFT JOIN cohort coh ON coh.UNIQUEID = core.id__UNIQUEID"
    ).fetchone()[0]
    n_R = sum(1 for y in y_bin if y == 1)
    n_S = len(y_bin) - n_R

    meta = {
        # The FE version and settings, recorded by the stage that applies them.
        # `mart_version` below is NOT the FE version: it is the scope label the
        # filename carries (every real run passes "vFULL" so downstream globs
        # find the file), and it overwrote the only FE version this project had.
        # fe_identity.identity() assembles fe_id from this block.
        "fe_version": fe_identity.FE_VERSION,
        "fe": fe_identity.producer_block("mart", {
            "mutation_selector": mutation_selector,
            "top_n_mutations": top_n,
            "max_carrier_frac": max_carrier_frac,
            "min_carrier_count": min_carrier_count,
            "gene_level_flags": gene_level_flags,
            "fold_strategy": fold_strategy,
            "n_folds": n_folds,
            "random_state": random_state,
        }, {"determinant_genes": {k: sorted(v) for k, v in sorted(DETERMINANT_GENES.items())}}),
        "mart_version": mart_version,
        "drug": drug,
        "cryptic_version": cryptic_version,
        # False means the release has no SUBLINEAGE column and
        # cov__sublineage is empty for every row here -- a feature-space
        # difference between releases, recorded so a cross-release
        # comparison can account for it instead of discovering it.
        # Empty means the release carried every optional genome covariate.
        # A non-empty list is a feature-space difference between releases,
        # recorded so a cross-release comparison can account for it
        # instead of discovering it.
        "missing_genome_covariates": missing_covariates,
        "missing_mutation_columns": missing_mutation_columns,
        "minor_allele_column": _minor,
        "cryptic_provenance": _read_db_metadata(con),
        "n_samples": len(y_bin),
        "n_R": n_R,
        "n_S": n_S,
        "pct_R": round(100.0 * n_R / len(y_bin), 2),
        "catalogue": {
            "name": CATALOGUE_NAME,
            "FP_count": cat_FP,
            "FN_count": cat_FN,
        },
        "n_columns_total": len(all_cols),
        "n_columns_per_section": {k: len(v) for k, v in sections.items()},
        "source_db": {
            "path": str(db_path),
            "checksum_sha256": _sha256(db_path),
        },
        "fold_assignment": {
            "method": (
                "stratified_round_robin_stdlib_random" if fold_strategy == "stratified"
                else f"grouped_by_{group_col}"
            ),
            "strategy": fold_strategy,
            "group_column": group_col,
            "n_splits": n_folds,
            "random_state": random_state,
            "stratify_on": "y__binary" if fold_strategy == "stratified" else None,
            "column": "__fold__",
            "fold_class_counts": {
                str(f): {
                    "R": sum(1 for i, fold_i in enumerate(folds) if fold_i == f and y_bin[i] == 1),
                    "S": sum(1 for i, fold_i in enumerate(folds) if fold_i == f and y_bin[i] == 0),
                }
                for f in range(n_folds)
            },
        },
        "leakage_boundary": {
            "de_novo_arm_filter": "drop columns starting with cat__",
            "catalogue_informed_arm_filter": "include cat__ tier-count columns; drop cat__who_*",
            "never_features": [
                "cat__who_predict (not computed; would be the catalogue prediction itself)",
                "cat__who_p (not computed)",
                "cat__catalogue_disagrees_with_label (not computed; uses label)",
            ],
            "tier_grading_note": (
                "design/09 spec called for cat__n_tier1/tier2/interim but the slim DB's "
                "effects table is trimmed (no EVIDENCE / PREDICTION_VALUES dicts) so tier "
                "grading isn't recoverable. cat__n_per_mut_{R,S,U,F} are the honest "
                "structural substitutes."
            ),
        },
        "build": {
            "wall_time_seconds": round(time.time() - t0, 2),
            "top_n_mutations": top_n,
            "carrier_rate_window": {
                "max_carrier_frac": max_carrier_frac,
                "min_carrier_count": min_carrier_count,
                "mutation_selector": mutation_selector,
                "gene_level_flags": gene_level_flags,
                "fold_strategy": fold_strategy,
                "fe_name": fe_name,
                "note": (
                    "Deviation from design/09: pre-filter mutations by carrier "
                    "rate before top-N. Pure prevalence ranking picks up only "
                    "lineage markers on the slim DB. Window keeps the canonical "
                    "resistance variants (rpoB/S450L, katG/S315T, embB/M306V, "
                    "gyrA/A90V, etc.) in the feature set."
                ),
            },
            "deferred_sections": ["mdl__", "causal__", "cf__"],
        },
        "outputs": {
            "parquet": mart_filename,
            "parquet_sha256": _sha256(mart_path),
        },
    }
    meta_path = mart_path.with_suffix(".metadata.json")
    meta_path.write_text(json.dumps(meta, indent=2))
    return meta


def main() -> None:
    parser = argparse.ArgumentParser(description="Per-drug Path-1 feature mart (pass 1).")
    parser.add_argument("--drug", required=True, help="DRUG code (e.g. RIF)")
    parser.add_argument("--db", type=Path, default=Path("analysis/databases/duckdb/cryptic-slim.duckdb"))
    parser.add_argument("--out", type=Path, default=Path("analysis/results/feature_mart"))
    parser.add_argument("--version", default=MART_VERSION_DEFAULT)
    parser.add_argument("--top-n", type=int, default=TOP_N_MUTATIONS)
    parser.add_argument("--mutation-selector", default="top_n",
                        choices=["top_n", "all"],
                        help="label-free ranking only; chi2/mi are FE steps that run "
                             "per fold, downstream of the label boundary")
    parser.add_argument("--max-carrier-frac", type=float, default=MAX_CARRIER_FRAC)
    parser.add_argument("--min-carrier-count", type=int, default=MIN_CARRIER_COUNT)
    parser.add_argument("--fold-strategy", default="stratified",
                        choices=["stratified", "lineage", "country"])
    parser.add_argument("--n-folds", type=int, default=N_FOLDS)
    parser.add_argument("--gene-level-flags", dest="gene_level_flags",
                        action="store_true", default=True)
    parser.add_argument("--no-gene-level-flags", dest="gene_level_flags",
                        action="store_false")
    parser.add_argument("--fe-name", default="default",
                        help="label for this FE configuration; lands in the sidecar")
    args = parser.parse_args()

    meta = build_mart(
        drug=args.drug,
        db_path=args.db,
        out_dir=args.out,
        mart_version=args.version,
        top_n=args.top_n,
        n_folds=args.n_folds,
        mutation_selector=args.mutation_selector,
        max_carrier_frac=args.max_carrier_frac,
        min_carrier_count=args.min_carrier_count,
        gene_level_flags=args.gene_level_flags,
        fold_strategy=args.fold_strategy,
        fe_name=args.fe_name,
    )
    print(json.dumps(meta, indent=2))


if __name__ == "__main__":
    main()
