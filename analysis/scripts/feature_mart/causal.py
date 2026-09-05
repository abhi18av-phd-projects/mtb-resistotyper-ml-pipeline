"""Causal-concordance feature selection — R1 experiment (design/15).

Implements the milestone design/15 §Status & resume marker prescribes: one
drug (RIF), T1 tier, no cluster. Family A (causal identification) votes via
EconML CausalForestDML CATE + DoWhy refutation; Family B (glass-box) votes
via EBM main-effect importance. A feature is `causal_concordant` when >=2 of
the 3 votes agree it is a driver (voting mode, per the operator decision —
design/15 leaves strict-AND vs voting open and asks for the sensitivity to
be reported; this is the voting-mode arm).

Feature granularity — CANDIDATE_LEVEL — two supported modes:
  "gene"     raw__gene_<gene> any-mutation-in-gene flags (default; scaled-up
             pool per user direction — coarser per candidate, covers more
             biological ground, 742 gene-level columns available in the
             RIF mart vs the original 60-mutation pilot)
  "mutation" raw__<gene>_<mutation> per-variant flags (the original R1
             pilot scope — kept available, not deleted, since it's already
             validated: see the git history / commit for the G3 fix that
             made it surface rpoB/S450L and katG/S315T correctly)

Candidate pool (this run only — not the deployed mart's feature set):
top-CANDIDATE_POOL_SIZE columns (at the chosen level) by |correlation| with
y__binary — G3 selection (design/02), NOT raw prevalence. Raw prevalence
within the mart's already-windowed columns still favours near-universal
lineage-linked passengers over resistance drivers (verified: rpoB/S450L
ranks 938/1000 by prevalence but has 94.7% R-rate among carriers vs 38.5%
base rate) — G3 label-correlation is what actually surfaces known drivers.

Votes, each independently computed per candidate:
  ebm_important    — candidate's EBM main-effect importance ranks in the
                      top EBM_TOP_K of the candidate pool (kept at the same
                      50% ratio as the original top-30-of-60 pilot)
  cate_significant — EconML CausalForestDML's 95% CI on the ATE excludes 0
  refutation_pass  — DoWhy placebo-treatment refuter does not invalidate
                      the backdoor-linear-regression effect estimate
                      (p > 0.05 on the refuter's own significance test)

concordant = (ebm_important + cate_significant + refutation_pass) >= 2

Every stochastic voter is seeded per design/15 R1:
  EconML CausalForestDML   random_state=42
  DoWhy placebo refuter    random_seed=42
  EBM                      random_state=42

WHAT THIS SCRIPT IS NOT (yet — see design/15 §Reproducibility, R2-R5):
  - NOT fold-internal. This fits on the full mart (no train/test split,
    no nested CV). Design/15's "no selection leakage" discipline requires
    fold-internal selection before any AUC claim is honest — that's the
    R2 hardening step, once this R1 mechanism proof is validated.
  - NOT recording realised-selection artifacts to MinIO, NOT
    content-addressing the CRyPTIC input, NOT pinning a container digest.
  - Family C (bacterial GWAS / phylogenetics) is out of scope per design/15
    (deferred until whole-genome + phylogeny inputs are available).

Environment note (2026-07): econml + dowhy could not be installed locally
this session — pixi's solver crashes on any new dependency (reproduced
independently with duckdb-cli and again with econml/dowhy; looks like a
rattler-solve 1.4.4 bug, not package-specific), and a pip/uv-based
workaround hit a broken native C toolchain building dowhy's cvxpy->scs
dependency. `pixi.toml` was NOT changed to avoid leaving a broken lock.
This module is written against design/15's spec and is intended to be run
once econml/dowhy are available (e.g. on the abc-cluster / a workbench
with a working causal-inference stack) — see the module-level TODO below
for the two pixi.toml lines to add first.

TODO before running: add to pixi.toml under [dependencies]:
    econml = "*"
    dowhy = "*"
then `pixi lock` (or, if the solver crash recurs, install into the pixi
env's Python directly with a working pip/uv on that machine).

Run (gene-level, pool 300 — the current default level):
    python -m analysis.scripts.feature_mart.causal \
        --mart analysis/results/feature_mart/feature_mart_RIF_cryptic-slim-2026.05_v1.0.1.parquet \
        --out analysis/results/causal/RIF/

Run (mutation-level, pool 300 — the sharper per-variant signal):
    python -m analysis.scripts.feature_mart.causal --level mutation \
        --mart analysis/results/feature_mart/feature_mart_RIF_cryptic-slim-2026.05_v1.0.1.parquet \
        --out analysis/results/causal/RIF_mutation_level/

Run (original per-mutation pilot scope — pool 60):
    python -m analysis.scripts.feature_mart.causal --level mutation --pool-size 60 \
        --mart analysis/results/feature_mart/feature_mart_RIF_cryptic-slim-2026.05_v1.0.1.parquet \
        --out analysis/results/causal/RIF_mutation_pilot/

Both levels default to a 300-candidate pool. Gene-level (raw__gene_<gene>)
is coarser — rpoB dilutes to ~0.49 R-rate because it counts any mutation
in the gene. Mutation-level (raw__<gene>_<mutation>) isolates the exact
resistance variant — rpoB/S450L is 0.95 R-rate, katG/S315T 0.76 — so the
concordance test has a much sharper true-driver-vs-hitchhiker contrast to
resolve (a steep drop to a ~0.46 lineage-marker cluster below the top 2).
"""

from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats

try:
    from econml.dml import CausalForestDML
except ImportError:  # pragma: no cover — surfaced at runtime, see module docstring
    CausalForestDML = None

try:
    import dowhy
    from dowhy import CausalModel
except ImportError:  # pragma: no cover
    dowhy = None
    CausalModel = None

try:
    from interpret.glassbox import ExplainableBoostingClassifier
except ImportError:  # pragma: no cover
    ExplainableBoostingClassifier = None


RANDOM_STATE = 42  # design/15 R1: every stochastic voter takes this seed
CANDIDATE_LEVEL = "gene"  # "gene" (raw__gene_<gene>) or "mutation" (raw__<gene>_<mutation>)
CANDIDATE_POOL_SIZE = 300  # widened per user direction — gene-level scale-up (was 60, mutation-level pilot)
MIN_CARRIER_COUNT_FOR_CANDIDACY = 10  # exclude singleton/noise candidates

# Canonical TB resistance determinant genes per drug. Their eligible candidate
# columns are ALWAYS seeded into the pool (union with the top-N-by-correlation),
# so a small/fast pool never silently drops the drug's true driver — e.g. EMB's
# embB, which a pool=60 correlation cut missed (design/16 fan-out finding).
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
    # Nitroimidazole activation pathway (DLM/PMD), per in-vitro selection spectrum
    # (Int J Antimicrob Agents 2026, 10.1016/j.ijantimicag.2026.107789): six genes,
    # SNP/LoF/indel all causal, extensive DLM<->PMD cross-resistance. Rv0678 dropped
    # — it is a BDQ/CFZ efflux regulator that appears in DLM-R isolates only as MDR
    # co-occurrence (an efflux confounder), not a nitroimidazole determinant.
    "DLM": {"ddn", "fgd1", "fbiA", "fbiB", "fbiC", "fbiD"},
}
EBM_TOP_K = 150  # 50% of pool — same ratio as the original top-30-of-60 mutation-level design
CONCORDANCE_K = 2  # of 3 votes — voting mode per the operator decision
CATE_N_ESTIMATORS = 200  # CausalForestDML trees; lower (e.g. 50) for a faster multi-drug sweep


def _require_libs() -> None:
    missing = [
        name
        for name, mod in [("econml", CausalForestDML), ("dowhy", CausalModel), ("interpret", ExplainableBoostingClassifier)]
        if mod is None
    ]
    if missing:
        raise SystemExit(
            f"Missing required libraries: {', '.join(missing)}. "
            "See this module's docstring — add econml/dowhy to pixi.toml "
            "and run on a machine with a working pixi/pip install first."
        )


def _candidate_columns(df: pd.DataFrame, level: str) -> list[str]:
    if level == "gene":
        return [c for c in df.columns if c.startswith("raw__gene_")]
    if level == "mutation":
        return [
            c
            for c in df.columns
            if c.startswith("raw__")
            and c not in {"raw__n_mutations", "raw__mean_frs", "raw__mean_coverage", "raw__n_minor"}
            and not c.startswith("raw__gene_")
        ]
    raise ValueError(f"level must be 'gene' or 'mutation', got {level!r}")


def _determinant_seed(df: pd.DataFrame, eligible, drug: str, level: str) -> list[str]:
    """Eligible candidate columns whose gene is a known determinant of `drug`.

    Guarantees the drug's true drivers are evaluated regardless of pool size.
    Mutation level: raw__<gene>_<mut>; gene level: raw__gene_<gene>.
    """
    genes = DETERMINANT_GENES.get(drug, set())
    if not genes:
        return []
    seed = []
    for c in eligible:
        if level == "gene":
            gene = c[len("raw__gene_"):]
        else:
            gene = c[len("raw__"):].split("_", 1)[0]
        if gene in genes:
            seed.append(c)
    return seed


def _select_candidate_pool(
    df: pd.DataFrame,
    pool_size: int = CANDIDATE_POOL_SIZE,
    level: str = CANDIDATE_LEVEL,
    drug: str | None = None,
) -> list[str]:
    """Top-N candidate columns by |correlation| with y__binary (G3, design/02).

    The mart's own top-1000 raw__ columns are already G2-filtered (frequency
    + carrier-rate window, build_mart.py) but the resulting pool still spans
    a wide prevalence range (10th to 94th percentile carrier rate) — mostly
    lineage-linked passengers, not resistance drivers. Ranking by raw
    prevalence within that pool reselects the same problem one level up
    (the most prevalent survivors are still disproportionately lineage
    markers). G3 (point-biserial correlation against the label) is the
    correct next filter in design/15's documented pipeline ("G2 frequency
    top-N -> G3 chi-square/MI prefilter -> EBM proposes...") and is what
    actually surfaces known drivers: at the mutation level, rpoB/S450L sits
    at prevalence-rank 938/1000 in the RIF mart but has among the strongest
    label associations (94.7% R-rate among carriers vs 38.5% base rate).
    """
    candidate_cols = _candidate_columns(df, level)
    carrier_counts = df[candidate_cols].sum(axis=0)
    eligible = carrier_counts[carrier_counts >= MIN_CARRIER_COUNT_FOR_CANDIDACY].index

    y = df["y__binary"]
    association = df[eligible].corrwith(y).abs().sort_values(ascending=False)
    pool = list(association.head(pool_size).index)
    # Seed determinant-gene candidates that the top-N-by-correlation cut missed.
    for c in _determinant_seed(df, eligible, drug or "", level):
        if c not in pool:
            pool.append(c)
    return pool


def _confounder_matrix(df: pd.DataFrame) -> pd.DataFrame:
    """Causal confounder set: lineage one-hots + coverage. Co-resistance is
    DELIBERATELY EXCLUDED.

    We tried adding cov__coresist_* (DST for the other first-line drugs) to
    deconfound the MDR co-selection finding (katG/S315T "predicts" RIF only
    because RIF-R strains are ~92% INH-R). It backfired: co-resistance is a
    *downstream consequence* of the outcome, not an upstream confounder —
    rpoB/S450L -> RIF-R -> (co-selection) -> INH-R. Conditioning on INH-R is
    therefore collider / over-control bias: it partially conditions on
    RIF-resistance itself and collapses the TRUE driver's effect (rpoB ATE
    0.75 -> 0.20, CATE significance lost), not just the hitchhiker's. Since
    rpoB (own-drug) and katG (cross-drug) have opposite causal relationships
    to co-resistance, no single "condition on all co-resistance" rule is
    correct for both.

    So co-resistance stays in the mart as a *covariate + diagnostic*
    (the 92%-vs-22% MDR signal is real and useful) but is NOT a causal
    confounder here. The katG-style cross-drug hitchhiker is a documented
    known limitation for the R2 fold-internal stage — deconfounding it
    correctly needs per-mutation causal structure (own-drug vs cross-drug),
    e.g. the attenuation-ratio diagnostic, not blanket conditioning.
    """
    cols = ["cov__lineage_L1", "cov__lineage_L2", "cov__lineage_L3", "cov__lineage_L4", "cov__median_coverage"]
    return df[cols].fillna(df[cols].median(numeric_only=True))


def _fit_ebm_importance(df: pd.DataFrame, candidates: list[str]) -> dict[str, float]:
    """Fit one EBM on the candidate pool + confounders; return per-candidate main-effect importance."""
    X = pd.concat([df[candidates], _confounder_matrix(df)], axis=1)
    y = df["y__binary"].to_numpy()

    ebm = ExplainableBoostingClassifier(interactions=0, random_state=RANDOM_STATE)
    ebm.fit(X, y)

    importances = dict(zip(ebm.term_names_, ebm.term_importances()))
    # Keep only the candidate-column importances (drop confounder terms).
    return {c: importances[c] for c in candidates if c in importances}


def _ebm_votes(importances: dict[str, float], top_k: int = EBM_TOP_K) -> dict[str, bool]:
    ranked = sorted(importances, key=importances.get, reverse=True)
    top_set = set(ranked[:top_k])
    return {c: (c in top_set) for c in importances}


def _cate_significant(df: pd.DataFrame, candidate_col: str,
                      n_estimators: int = CATE_N_ESTIMATORS) -> tuple[bool, float, float, float]:
    """EconML CausalForestDML: treatment=candidate indicator, outcome=y__binary.

    Returns (significant, ate, ci_lower, ci_upper). Significant when the
    95% CI on the average treatment effect excludes 0.

    CausalForestDML is a CATE (heterogeneous-effect) estimator and requires
    X (effect-modifier features) — it rejects X=None. We pass the confounder
    matrix as X, which lets the mutation's effect vary across lineage /
    coverage (design/15's "effect invariant across lineages" intent) and
    then averages those per-sample CATEs into the ATE + CI. The same matrix
    is adjusted for as a confounder either way, so this does not weaken the
    deconfounding — it only permits (and measures) heterogeneity.

    `n_estimators` trades runtime for CI precision. The RIF experiment used
    200 (44 min for 300 candidates); significance (CI excludes 0) is robust
    for strong effects, so a smaller forest (e.g. 50) is a defensible ~4x
    speedup for multi-drug sweeps — the borderline calls shift, not the
    clear drivers.
    """
    X = _confounder_matrix(df).to_numpy()
    T = df[candidate_col].to_numpy()
    Y = df["y__binary"].to_numpy().astype(float)

    # econml requires n_estimators divisible by subforest_size (default 4);
    # round to the nearest multiple of 4 so any --cate-estimators works.
    n_est = max(4, round(n_estimators / 4) * 4)
    est = CausalForestDML(
        n_estimators=n_est,
        random_state=RANDOM_STATE,
        discrete_treatment=True,
    )
    est.fit(Y, T, X=X)
    ate = float(est.ate(X))
    lower, upper = est.ate_interval(X, alpha=0.05)
    lower, upper = float(lower), float(upper)
    significant = not (lower <= 0.0 <= upper)
    return significant, ate, lower, upper


def _residual_signal_survives(df: pd.DataFrame, candidate_col: str) -> tuple[bool, float, float]:
    """Family-A lineage/confounder residualisation vote (design/15 spec'd).

    Frisch-Waugh-Lovell partial correlation: residualise both the candidate
    indicator and y__binary on the confounder matrix (lineage + coverage +
    co-resistance) via OLS, then test whether the residuals still correlate.
    If they do (significant partial correlation), the candidate's signal
    SURVIVES partialling out the confounders → a driver vote. A pure lineage
    marker's correlation with y IS the lineage, so its residual correlation
    collapses to ~0 → correctly not a driver. With co-resistance in the
    confounders (mart v1.0.2+), an MDR-hitchhiker like katG/S315T should
    likewise collapse once INH-co-resistance is partialled out.

    Returns (survives, partial_r, p_value). Survives when the two-sided
    partial-correlation p-value < 0.05 AND |partial_r| exceeds a small floor
    (0.05) so a statistically-significant-but-negligible correlation on
    14k rows doesn't trivially pass.

    This REPLACES the DoWhy placebo refuter as the third vote: refutation
    tests the estimation METHOD's robustness (it returns an identical,
    non-discriminating verdict for every candidate — verified empirically),
    whereas residualisation tests THIS candidate's signal against the
    confounders, which is the causal-vs-hitchhiker question the vote needs.
    `_refutation_pass` is retained below as an optional GLOBAL method-validity
    gate (run once), not a per-candidate vote.
    """
    W = _confounder_matrix(df).to_numpy()
    W = np.column_stack([np.ones(len(W)), W])  # intercept
    T = df[candidate_col].to_numpy().astype(float)
    Y = df["y__binary"].to_numpy().astype(float)

    # OLS residuals of T and Y on the confounders (least-squares projection).
    beta_t, *_ = np.linalg.lstsq(W, T, rcond=None)
    beta_y, *_ = np.linalg.lstsq(W, Y, rcond=None)
    res_t = T - W @ beta_t
    res_y = Y - W @ beta_y

    if np.std(res_t) < 1e-12 or np.std(res_y) < 1e-12:
        return False, 0.0, 1.0

    partial_r = float(np.corrcoef(res_t, res_y)[0, 1])
    n, k = len(T), W.shape[1]
    dof = n - k - 1
    # t-statistic for a partial correlation with dof residual degrees of freedom.
    t_stat = partial_r * np.sqrt(dof / max(1e-12, 1.0 - partial_r**2))
    p_value = float(2.0 * stats.t.sf(abs(t_stat), dof))
    survives = (p_value < 0.05) and (abs(partial_r) >= 0.05)
    return survives, partial_r, p_value


def _refutation_pass(df: pd.DataFrame, candidate_col: str) -> tuple[bool, float]:
    """DoWhy placebo-treatment refuter — GLOBAL method-validity gate, NOT a vote.

    Retained for the once-per-run method sanity check design/15 expects
    (does the estimation pipeline hallucinate an effect from a randomised
    placebo treatment?). NOT used as a per-candidate concordance vote: the
    placebo refuter returns a byte-identical, treatment-independent verdict
    for every candidate (verified — p_value=0.4016853 for effects 0.75 and
    0.25 alike), because it probes the METHOD, not the feature. See
    `_residual_signal_survives` for the per-candidate Family-A vote that
    replaced it.

    Returns (passes, refuter_p_value). Passes when the placebo effect is not
    significantly different from zero (p > 0.05) — the method behaves as
    expected under a known-null substitution.
    """
    confounder_cols = list(_confounder_matrix(df).columns)
    data = df[[candidate_col, "y__binary"] + confounder_cols].fillna(
        df[confounder_cols].median(numeric_only=True)
    )

    model = CausalModel(
        data=data,
        treatment=candidate_col,
        outcome="y__binary",
        common_causes=confounder_cols,
    )
    identified = model.identify_effect(proceed_when_unidentifiable=True)
    estimate = model.estimate_effect(
        identified,
        method_name="backdoor.linear_regression",
    )
    refutation = model.refute_estimate(
        identified,
        estimate,
        method_name="placebo_treatment_refuter",
        random_seed=RANDOM_STATE,
        num_simulations=20,
    )
    p_value = float(refutation.refutation_result.get("p_value", 1.0))
    passes = p_value > 0.05
    return passes, p_value


def concordance_sensitivity(results: pd.DataFrame) -> dict:
    """Report how the concordant-set size depends on the rule (design/15 asks
    for this sensitivity explicitly; the operator picks strict-AND vs voting).

    The R1 finding: CATE-significance is the sole discriminating vote — EBM
    (top-K of the pool) and residual (partial-r floor) are permissive and
    add lineage-marker noise, so any CATE-requiring rule collapses the
    concordant set to the true drivers. Reported per-run so the choice of
    CONCORDANCE_K is grounded in data, not asserted.
    """
    votes = results["n_votes"]
    return {
        "voting_ge_1": int((votes >= 1).sum()),
        "voting_ge_2": int((votes >= 2).sum()),
        "strict_and_3": int((votes >= 3).sum()),
        "cate_required": int(results["cate_significant"].sum()),
        "cate_and_ebm": int((results["cate_significant"] & results["ebm_important"]).sum()),
        "cate_and_residual": int((results["cate_significant"] & results["residual_survives"]).sum()),
    }


def _build_manifest(results: pd.DataFrame, drug: str, level: str, pool_size: int, ebm_top_k: int,
                    wall_seconds: float, results_path: Path, out_dir: Path,
                    cate_estimators: int = CATE_N_ESTIMATORS) -> dict:
    return {
        "drug": drug,
        "stage": "R1 experiment (design/15)",
        "candidate_level": level,
        "candidate_pool_size": pool_size,
        "ebm_top_k": ebm_top_k,
        "cate_n_estimators": cate_estimators,
        "concordance_rule": f"voting mode, >={CONCORDANCE_K}-of-3",
        "n_concordant": int(results["qf__causal_concordant"].sum()),
        "concordance_sensitivity": concordance_sensitivity(results),
        "seeds": {"econml": RANDOM_STATE, "dowhy_refuter": RANDOM_STATE, "ebm": RANDOM_STATE},
        "not_yet_implemented": [
            "fold-internal selection (no selection leakage) — R2",
            "realised-selection artifacts to MinIO — R2",
            "content-addressed CRyPTIC input — R3",
            "container digest pinning — R4",
            "fold-assignment artifact — R5",
        ],
        "wall_time_seconds": round(wall_seconds, 2),
        "outputs": {"results_parquet": str(results_path.relative_to(out_dir.parent.parent))},
    }


def run_r1_experiment(
    mart_path: Path,
    out_dir: Path,
    drug: str = "RIF",
    level: str = CANDIDATE_LEVEL,
    pool_size: int = CANDIDATE_POOL_SIZE,
    ebm_top_k: int = EBM_TOP_K,
    cate_estimators: int = CATE_N_ESTIMATORS,
    held_out_lineage: str | None = None,
) -> dict:
    _require_libs()
    out_dir.mkdir(parents=True, exist_ok=True)
    t0 = time.time()

    df = pd.read_parquet(mart_path)

    # Fold-internal selection. When a lineage is held out, the voters never see
    # its rows, so the resulting feature set can be evaluated on that lineage
    # without the selection having already looked at it. Selecting on all data
    # and then holding a lineage out for scoring is selection leakage, and it
    # inflates exactly the number this project reports as honest.
    n_all = len(df)
    if held_out_lineage:
        if "cov__lineage_raw" not in df.columns:
            raise ValueError("mart has no cov__lineage_raw; cannot hold a lineage out")
        df = df.loc[df["cov__lineage_raw"] != held_out_lineage].reset_index(drop=True)
        if df.empty:
            raise ValueError(f"holding out {held_out_lineage} left no rows")
    candidates = _select_candidate_pool(df, pool_size=pool_size, level=level, drug=drug)

    importances = _fit_ebm_importance(df, candidates)
    ebm_votes = _ebm_votes(importances, top_k=ebm_top_k)

    rows = []
    for candidate in candidates:
        ebm_vote = ebm_votes.get(candidate, False)
        cate_sig, ate, ci_lo, ci_hi = _cate_significant(df, candidate, n_estimators=cate_estimators)
        resid_survives, partial_r, resid_p = _residual_signal_survives(df, candidate)

        # Three discriminating per-candidate votes across Families A + B:
        #   Family B: ebm_important (glass-box predictive importance)
        #   Family A: cate_significant (CausalForestDML ATE CI excludes 0)
        #   Family A: residual_survives (FWL partial corr after confounders)
        votes = int(ebm_vote) + int(cate_sig) + int(resid_survives)
        concordant = votes >= CONCORDANCE_K

        rows.append(
            {
                "feature_column": candidate,
                "ebm_important": ebm_vote,
                "ebm_importance": importances.get(candidate),
                "cate_significant": cate_sig,
                "cate_ate": ate,
                "cate_ci_lower": ci_lo,
                "cate_ci_upper": ci_hi,
                "residual_survives": resid_survives,
                "residual_partial_r": partial_r,
                "residual_p_value": resid_p,
                "n_votes": votes,
                "qf__causal_concordant": concordant,
            }
        )

    results = pd.DataFrame(rows).sort_values("n_votes", ascending=False)
    results_path = out_dir / f"causal_concordance_{drug}_{level}_level.parquet"
    results.to_parquet(results_path, index=False)

    manifest = _build_manifest(
        results, drug, level, len(candidates), ebm_top_k, time.time() - t0, results_path, out_dir,
        cate_estimators=cate_estimators,
    )
    # Which arm produced this selection is not a detail. A feature set chosen on
    # all data and one chosen without the held-out lineage support different
    # claims, and nothing downstream can tell them apart unless the artefact says.
    manifest["selection_arm"] = "fold-internal" if held_out_lineage else "refit-on-all"
    manifest["held_out_lineage"] = held_out_lineage
    manifest["n_rows_used"] = int(len(df))
    manifest["n_rows_total"] = int(n_all)
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def reprocess_manifest(results_path: Path, out_dir: Path) -> dict:
    """Add the concordance sensitivity to an existing run's manifest.json
    from its results parquet — no 44-minute vote recompute. Preserves the
    original wall_time, catalogue_benchmark, etc.; only adds/refreshes the
    sensitivity block. Used to enrich prior runs after adding the field."""
    results = pd.read_parquet(results_path)
    manifest_path = out_dir / "manifest.json"
    manifest = json.loads(manifest_path.read_text()) if manifest_path.exists() else {}
    manifest["concordance_sensitivity"] = concordance_sensitivity(results)
    manifest["reprocessed"] = True
    manifest_path.write_text(json.dumps(manifest, indent=2))
    return manifest


def _parse_candidate(column: str, level: str = CANDIDATE_LEVEL) -> tuple[str, str | None]:
    """Gene-level: `raw__gene_rpoB` -> (`rpoB`, None).
    Mutation-level: `raw__rpoB_S450L` -> (`rpoB`, `S450L`)."""
    if level == "gene":
        return column.removeprefix("raw__gene_"), None
    body = column.removeprefix("raw__")
    gene, _, mutation = body.partition("_")
    return gene, mutation


def benchmark_vs_catalogue(
    results: pd.DataFrame,
    db_path: Path,
    drug: str = "RIF",
    catalogue_name: str = "WHO-UCN-GTB-PCI-2023.5",
    level: str = CANDIDATE_LEVEL,
) -> pd.DataFrame:
    """Join the concordant set against the WHO catalogue's verdict.

    Gene level: does ANY mutation in this gene carry a catalogue verdict for
    this drug? Reports counts of catalogue R/total mutations per gene, since
    a single (gene, mutation) match doesn't make sense at this granularity.
    Mutation level: exact (gene, mutation) catalogue verdict, as in the
    original R1 pilot.

    Uses only DuckDB — no econml/dowhy/interpret dependency — so this half
    of the R1 deliverable ("benchmarked vs the WHO-catalogue baseline") is
    runnable and testable independent of the causal-library environment
    blocker documented in this module's docstring.
    """
    import duckdb

    parsed = [_parse_candidate(c, level) for c in results["feature_column"]]
    lookup = pd.DataFrame(parsed, columns=["GENE", "MUTATION"])
    lookup["feature_column"] = results["feature_column"].values

    con = duckdb.connect(str(db_path), read_only=True)
    con.register("candidates", lookup)

    if level == "gene":
        catalogue = con.execute(
            """
            SELECT c.feature_column,
                   COUNT(DISTINCT e.MUTATION)                                             AS n_catalogue_mutations,
                   COUNT(DISTINCT CASE WHEN e.PREDICTION = 'R' THEN e.MUTATION END)        AS n_catalogue_R_mutations
            FROM   candidates c
            LEFT JOIN effects e
              ON  e.GENE = c.GENE
              AND e.DRUG = ? AND e.CATALOGUE_NAME = ?
            GROUP BY c.feature_column
            """,
            [drug, catalogue_name],
        ).fetchdf()
        catalogue["catalogue_verdict"] = np.where(
            catalogue["n_catalogue_R_mutations"].fillna(0) > 0,
            "has_catalogue_R_mutation",
            "no_catalogue_R_mutation",
        )
        merge_cols = ["feature_column", "catalogue_verdict", "n_catalogue_mutations", "n_catalogue_R_mutations"]
    else:
        catalogue = con.execute(
            """
            SELECT c.feature_column,
                   e.PREDICTION AS catalogue_verdict
            FROM   candidates c
            LEFT JOIN effects e
              ON  e.GENE = c.GENE AND e.MUTATION = c.MUTATION
              AND e.DRUG = ? AND e.CATALOGUE_NAME = ?
            """,
            [drug, catalogue_name],
        ).fetchdf()
        # A mutation can appear in multiple isolates' effect rows with the
        # same verdict; collapse to one row per mutation (verdict is
        # deterministic per (gene, mutation, catalogue_name)).
        catalogue = catalogue.drop_duplicates(subset="feature_column")
        catalogue["catalogue_verdict"] = catalogue["catalogue_verdict"].fillna("not_in_catalogue")
        merge_cols = ["feature_column", "catalogue_verdict"]

    merged = results.merge(catalogue[merge_cols], on="feature_column", how="left")
    merged["catalogue_verdict"] = merged["catalogue_verdict"].fillna(
        "no_catalogue_R_mutation" if level == "gene" else "not_in_catalogue"
    )
    return merged


def _infer_drug(mart_path: Path) -> str:
    """feature_mart_<DRUG>_cryptic-... -> <DRUG>."""
    name = mart_path.name
    if name.startswith("feature_mart_"):
        return name[len("feature_mart_"):].split("_")[0]
    return "RIF"


def main() -> None:
    parser = argparse.ArgumentParser(description="Causal-concordance experiment (design/15), any drug.")
    parser.add_argument("--mart", type=Path, required=True, help="Path to the per-drug feature mart Parquet.")
    parser.add_argument("--out", type=Path, default=None, help="Output dir (default: analysis/results/causal/<DRUG>_<level>).")
    parser.add_argument("--drug", default=None, help="Drug label (inferred from the mart filename if omitted).")
    parser.add_argument("--level", choices=["gene", "mutation"], default=CANDIDATE_LEVEL)
    parser.add_argument("--pool-size", type=int, default=None, help="Override the candidate pool size.")
    arm = parser.add_mutually_exclusive_group()
    arm.add_argument("--held-out-lineage", default=None,
                     help="Withhold this lineage from selection: the EVALUATION arm. "
                          "Without it the voters see every row, which is only correct "
                          "for the deployment refit.")
    arm.add_argument("--refit-on-all", dest="held_out_lineage",
                     action="store_const", const=None,
                     help="Explicit DEPLOYMENT arm: select on all data (the default).")
    parser.add_argument("--cate-estimators", type=int, default=CATE_N_ESTIMATORS,
                        help="CausalForestDML trees (lower = faster multi-drug sweep).")
    parser.add_argument(
        "--db",
        type=Path,
        default=Path("analysis/databases/duckdb/cryptic-slim.duckdb"),
        help="Slim DB path, for the WHO-catalogue benchmark join.",
    )
    parser.add_argument(
        "--skip-catalogue-benchmark",
        action="store_true",
        help="Skip the WHO-catalogue join (e.g. if the DB isn't available locally).",
    )
    args = parser.parse_args()

    drug = args.drug or _infer_drug(args.mart)
    out_dir = args.out or Path(f"analysis/results/causal/{drug}_{args.level}")
    pool_size = args.pool_size or CANDIDATE_POOL_SIZE  # 300 at both levels; --pool-size 60 -> original mutation pilot
    ebm_top_k = pool_size // 2  # preserve the 50% ratio regardless of level/pool size

    manifest = run_r1_experiment(args.mart, out_dir, drug=drug, level=args.level,
                                 held_out_lineage=args.held_out_lineage,
                                 pool_size=pool_size, ebm_top_k=ebm_top_k,
                                 cate_estimators=args.cate_estimators)

    if not args.skip_catalogue_benchmark:
        results = pd.read_parquet(out_dir / f"causal_concordance_{drug}_{args.level}_level.parquet")
        benchmarked = benchmark_vs_catalogue(results, args.db, drug=drug, level=args.level)
        benchmarked.to_parquet(out_dir / f"causal_concordance_{drug}_{args.level}_level_vs_catalogue.parquet", index=False)
        concordant = benchmarked[benchmarked["qf__causal_concordant"]]
        manifest["catalogue_benchmark"] = {
            "n_concordant": len(concordant),
            "n_concordant_with_catalogue_signal": int(
                concordant["catalogue_verdict"].isin(["R", "has_catalogue_R_mutation"]).sum()
            ),
            "n_concordant_novel": int(
                concordant["catalogue_verdict"].isin(["not_in_catalogue", "no_catalogue_R_mutation"]).sum()
            ),
        }
        (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))

    print(json.dumps(manifest, indent=2))


if __name__ == "__main__":
    main()
