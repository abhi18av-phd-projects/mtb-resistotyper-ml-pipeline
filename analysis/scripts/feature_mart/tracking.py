"""MLflow experiment tracking for the FE -> metric chain (garden template).

The garden template (`template-analysis-and-writeup`) chose MLflow as its
`experiment_tracking` tool and already ships an `mlflow ui` task
(`analysis/analysis.just:37`). This module is the bridge: it reads the
committed run manifests (feature-mart builds + causal-concordance runs) and
logs their feature-engineering *params* and resulting *metrics* to MLflow,
so `mlflow ui` and the Quarto dashboard
(`analysis/data/08_reporting/083_dashboards/`) can answer
"which feature-engineering choice moved which metric".

Design decisions:
  - ONE experiment ("mtb-fe-to-metric") spanning the full FE-decision chain;
    a `stage` tag ("feature_mart" | "causal_concordance") separates the
    upstream mart-build decisions from the downstream causal-selection
    decisions (the user's "full FE-decision chain" granularity).
  - MLflow is a regenerable INDEX, not the source of truth. The committed
    manifest.json files are the source of truth; `mlruns/` is gitignored
    and rebuilt by `backfill()`. Same treatment as the result parquets
    (gitignored, manifests committed) — a collaborator reproduces the
    dashboard with `backfill()` then `mlflow ui` / `quarto render`.
  - Tracking backend is sqlite (`sqlite:///mlflow.db`), NOT the legacy file
    store. MLflow 3.x put the `./mlruns` file store in maintenance mode and
    refuses to write to it; the template's original `mlflow ui` task (from
    the MLflow 1.x/2.x era) is updated in analysis/analysis.just to point at
    this sqlite backend. Both `mlflow.db` and `mlartifacts/` are gitignored
    and regenerable via `backfill()`.

Run (backfill every manifest found under analysis/results/):
    python -m analysis.scripts.feature_mart.tracking

Then browse:
    mlflow ui --backend-store-uri sqlite:///mlflow.db
    # or: just -f analysis/analysis.just mlflow-ui
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import mlflow

EXPERIMENT_NAME = "mtb-fe-to-metric"
# sqlite backend — MLflow 3.x rejects the legacy ./mlruns file store.
DEFAULT_TRACKING_URI = "sqlite:///mlflow.db"


def _run_name(manifest: dict, stage: str) -> str:
    if stage == "feature_mart":
        return f"mart:{manifest.get('drug')}:{manifest.get('mart_version')}"
    return (
        f"causal:{manifest.get('drug')}:{manifest.get('candidate_level')}"
        f":pool{manifest.get('candidate_pool_size')}"
    )


def _log_feature_mart(manifest: dict, manifest_path: Path) -> str:
    """Log a feature-mart build: FE params (top-N, carrier window, folds) ->
    metrics (sample/column counts, class balance, catalogue error counts)."""
    build = manifest.get("build", {})
    window = build.get("carrier_rate_window", {})
    sections = manifest.get("n_columns_per_section", {})

    with mlflow.start_run(run_name=_run_name(manifest, "feature_mart")) as run:
        mlflow.set_tag("stage", "feature_mart")
        mlflow.set_tag("drug", manifest.get("drug"))
        mlflow.set_tag("mart_version", manifest.get("mart_version"))
        mlflow.set_tag("manifest_path", str(manifest_path))

        # --- FE params (the choices) ---
        mlflow.log_params(
            {
                "drug": manifest.get("drug"),
                "mart_version": manifest.get("mart_version"),
                "cryptic_version": manifest.get("cryptic_version"),
                "fe_top_n_mutations": build.get("top_n_mutations"),
                "fe_carrier_max_frac": window.get("max_carrier_frac"),
                "fe_carrier_min_count": window.get("min_carrier_count"),
                "fe_n_folds": manifest.get("fold_assignment", {}).get("n_splits"),
                "random_state": manifest.get("fold_assignment", {}).get("random_state"),
                "catalogue": manifest.get("catalogue", {}).get("name"),
            }
        )

        # --- metrics (the outcomes) ---
        cat = manifest.get("catalogue", {})
        mlflow.log_metrics(
            {
                "n_samples": manifest.get("n_samples", 0),
                "n_R": manifest.get("n_R", 0),
                "n_S": manifest.get("n_S", 0),
                "pct_R": manifest.get("pct_R", 0.0),
                "n_columns_total": manifest.get("n_columns_total", 0),
                "n_raw_cols": sections.get("raw", 0),
                "n_cov_cols": sections.get("cov", 0),
                "catalogue_FP": cat.get("FP_count", 0),
                "catalogue_FN": cat.get("FN_count", 0),
                "build_wall_seconds": build.get("wall_time_seconds", 0.0),
            }
        )
        return run.info.run_id


def _log_causal_concordance(manifest: dict, manifest_path: Path) -> str:
    """Log a causal-concordance run: selection params (level, pool, rule,
    seeds) -> metrics (n_concordant, catalogue overlap, wall time)."""
    with mlflow.start_run(run_name=_run_name(manifest, "causal_concordance")) as run:
        mlflow.set_tag("stage", "causal_concordance")
        mlflow.set_tag("drug", manifest.get("drug"))
        mlflow.set_tag("candidate_level", manifest.get("candidate_level"))
        mlflow.set_tag("manifest_path", str(manifest_path))

        seeds = manifest.get("seeds", {})
        mlflow.log_params(
            {
                "drug": manifest.get("drug"),
                "candidate_level": manifest.get("candidate_level"),
                "candidate_pool_size": manifest.get("candidate_pool_size"),
                "ebm_top_k": manifest.get("ebm_top_k"),
                "concordance_rule": manifest.get("concordance_rule"),
                "seed_econml": seeds.get("econml"),
                "seed_dowhy": seeds.get("dowhy_refuter"),
                "seed_ebm": seeds.get("ebm"),
            }
        )

        bench = manifest.get("catalogue_benchmark", {})
        sens = manifest.get("concordance_sensitivity", {})
        metrics = {
            "n_concordant": manifest.get("n_concordant", 0),
            "candidate_pool_size": manifest.get("candidate_pool_size", 0),
            "wall_seconds": manifest.get("wall_time_seconds", 0.0),
        }
        if bench:
            metrics.update(
                {
                    "concordant_with_catalogue_signal": bench.get("n_concordant_with_catalogue_signal", 0),
                    "concordant_novel": bench.get("n_concordant_novel", 0),
                }
            )
        if sens:
            # The concordance-rule sensitivity — CATE is the discriminating
            # vote, so cate_required/strict_and collapse to the true drivers
            # while voting_ge_2 admits lineage-marker noise (design/15 R1).
            metrics.update(
                {
                    "concordant_voting_ge_2": sens.get("voting_ge_2", 0),
                    "concordant_strict_and_3": sens.get("strict_and_3", 0),
                    "concordant_cate_required": sens.get("cate_required", 0),
                }
            )
        mlflow.log_metrics(metrics)

        # The per-candidate result parquet (if present) is the drill-down
        # artifact the dashboard links to.
        for candidate in (
            manifest_path.parent / f"causal_concordance_RIF_{manifest.get('candidate_level')}_level.parquet",
        ):
            if candidate.exists():
                mlflow.log_artifact(str(candidate), artifact_path="results")
        return run.info.run_id


def _log_h2o_model(manifest: dict, manifest_path: Path) -> str:
    """Log an H2O stacked-ensemble run: feature-set + budget -> CV AUC etc.
    This is the row where feature engineering finally meets a real model
    metric (design/03, Paper 1)."""
    with mlflow.start_run(run_name=f"h2o:{manifest.get('drug')}:{manifest.get('feature_set')}") as run:
        mlflow.set_tag("stage", "h2o_model")
        mlflow.set_tag("drug", manifest.get("drug"))
        mlflow.set_tag("feature_set", manifest.get("feature_set"))
        mlflow.set_tag("manifest_path", str(manifest_path))

        mlflow.log_params(
            {
                "drug": manifest.get("drug"),
                "feature_set": manifest.get("feature_set"),
                "n_features": manifest.get("n_features"),
                "max_runtime_secs": manifest.get("max_runtime_secs"),
                "leader_model_id": manifest.get("cv_metrics", {}).get("leader_model_id"),
                "leader_is_stacked_ensemble": manifest.get("cv_metrics", {}).get("leader_is_stacked_ensemble"),
                "random_state": manifest.get("random_state"),
            }
        )
        cv = manifest.get("cv_metrics", {})
        surr = manifest.get("surrogate", {})
        metrics = {
            "cv_auc": cv.get("cv_auc"),
            "cv_logloss": cv.get("cv_logloss"),
            "n_leaderboard_models": cv.get("n_leaderboard_models"),
            "n_features": manifest.get("n_features"),
            "surrogate_fidelity_r2": surr.get("fidelity_r2"),
            "wall_seconds": manifest.get("wall_time_seconds"),
        }
        if cv.get("cv_aucpr") is not None:
            metrics["cv_aucpr"] = cv.get("cv_aucpr")
        mlflow.log_metrics({k: v for k, v in metrics.items() if v is not None})

        for name in ("leaderboard.csv", "surrogate_terms.csv"):
            p = manifest_path.parent / name
            if p.exists():
                mlflow.log_artifact(str(p), artifact_path="h2o")
        return run.info.run_id


def _log_cv_evaluation(manifest: dict, manifest_path: Path) -> list[str]:
    """Log the honest CV evaluation — one MLflow run per feature set, with
    the random / lineage-LOO / country-LOO AUCs. The lineage-LOO drop is
    the honesty number; the concordant-beats-denovo-under-lineage-CV flip is
    the headline (design/03)."""
    run_ids = []
    for fs, r in manifest.get("results", {}).items():
        with mlflow.start_run(run_name=f"cv:{manifest.get('drug')}:{fs}") as run:
            mlflow.set_tag("stage", "cv_evaluation")
            mlflow.set_tag("drug", manifest.get("drug"))
            mlflow.set_tag("feature_set", fs)
            mlflow.set_tag("manifest_path", str(manifest_path))
            mlflow.log_params({
                "drug": manifest.get("drug"), "feature_set": fs,
                "n_features": r.get("n_features"), "model": manifest.get("model"),
            })
            rnd = r.get("random", {}).get("overall_auc")
            lin = r.get("lineage_loo", {}).get("mean_group_auc")
            cty = r.get("country_loo", {}).get("mean_group_auc")
            metrics = {
                "auc_random": rnd,
                "auc_lineage_loo": lin,
                "auc_country_loo": cty,
                "lineage_leakage_inflation": (rnd - lin) if (rnd is not None and lin is not None) else None,
                "n_features": r.get("n_features"),
            }
            # per-held-out-lineage AUCs — the fairness/equity breakdown
            for lname, v in r.get("lineage_loo", {}).get("per_group", {}).items():
                if v.get("auc") is not None:
                    metrics[f"auc_heldout_{lname}"] = v["auc"]
            mlflow.log_metrics({k: v for k, v in metrics.items() if v is not None})
            run_ids.append(run.info.run_id)
    return run_ids


def _classify(manifest: dict) -> str | None:
    if "mart_version" in manifest:
        return "feature_mart"
    if manifest.get("stage", "").startswith("cv-evaluation"):
        return "cv_evaluation"
    if manifest.get("stage", "").startswith("h2o") or "feature_set" in manifest and "cv_metrics" in manifest:
        return "h2o_model"
    if manifest.get("stage", "").startswith("R1") or "candidate_level" in manifest:
        return "causal_concordance"
    return None


def backfill(results_root: Path, tracking_uri: str) -> list[str]:
    """Log every manifest.json / *.metadata.json under results_root to MLflow."""
    mlflow.set_tracking_uri(tracking_uri)
    mlflow.set_experiment(EXPERIMENT_NAME)

    manifests = sorted(results_root.rglob("manifest.json")) + sorted(
        results_root.rglob("*.metadata.json")
    )
    logged = []
    for path in manifests:
        try:
            manifest = json.loads(path.read_text())
        except json.JSONDecodeError:
            continue
        stage = _classify(manifest)
        if stage == "feature_mart":
            logged.append(_log_feature_mart(manifest, path))
        elif stage == "causal_concordance":
            logged.append(_log_causal_concordance(manifest, path))
        elif stage == "h2o_model":
            logged.append(_log_h2o_model(manifest, path))
        elif stage == "cv_evaluation":
            logged.extend(_log_cv_evaluation(manifest, path))
    return logged


def main() -> None:
    parser = argparse.ArgumentParser(description="Backfill FE->metric runs into MLflow (garden template).")
    parser.add_argument("--results", type=Path, default=Path("analysis/results"))
    parser.add_argument("--tracking-uri", default=DEFAULT_TRACKING_URI)
    args = parser.parse_args()

    logged = backfill(args.results, args.tracking_uri)
    print(f"Logged {len(logged)} runs to MLflow experiment '{EXPERIMENT_NAME}' "
          f"at {args.tracking_uri}")
    print(f"Browse:  mlflow ui --backend-store-uri {args.tracking_uri}   (http://localhost:5000)")


if __name__ == "__main__":
    main()
