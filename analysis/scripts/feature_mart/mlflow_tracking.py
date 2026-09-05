"""Log an H2O training run to MLflow, without letting tracking break training.

Two properties this module is built around.

Tracking is optional and must stay that way. A training run that took an hour of
cluster time is the valuable thing; a tracking server that is down, unreachable
from a Nomad task, or simply not configured must not lose it. Every function
here swallows its own failures and says so on stderr, and the caller's
manifest.json is written whether or not any of this worked. The manifest remains
the source of truth; MLflow is a view over it.

Tracking is written from inside the run, not backfilled afterwards. A backfill
records only what succeeded: the runs worth inspecting most closely are the ones
that failed halfway, and those never reach a manifest. Logging as the run
proceeds means a crashed run still leaves its parameters, its partial metrics
and its failure status behind.

The three H2O layers map onto MLflow's shape rather than being flattened into
it. Each parent run is one drug and one feature set; every base model and every
ensemble becomes a nested child run, so the leaderboard is navigable in the UI
instead of living in a CSV attached to a single blob.
"""

from __future__ import annotations

import contextlib
import json
import os
import sys
from pathlib import Path
from typing import Any

EXPERIMENT = os.environ.get("MLFLOW_EXPERIMENT", "mtb-resistotyper-ml/h2o")


def _warn(msg: str) -> None:
    print(f"[mlflow] {msg}", file=sys.stderr, flush=True)


def enabled() -> bool:
    return bool(os.environ.get("MLFLOW_TRACKING_URI"))


@contextlib.contextmanager
def run(drug: str, feature_set: str, params: dict[str, Any]):
    """A parent run for one (drug, feature set), or a no-op if tracking is off.

    Yields a handle whose methods are safe to call unconditionally: when
    tracking is unavailable every one of them does nothing, so the training
    script needs no `if tracking:` branches and cannot drift into a state where
    the tracked and untracked paths differ.
    """
    if not enabled():
        yield _Null()
        return
    try:
        import mlflow
    except ImportError:
        _warn("mlflow is not installed in this environment; not tracking")
        yield _Null()
        return

    try:
        mlflow.set_tracking_uri(os.environ["MLFLOW_TRACKING_URI"])
        mlflow.set_experiment(EXPERIMENT)
        active = mlflow.start_run(run_name=f"{drug}-{feature_set}")
    except Exception as exc:
        _warn(f"could not start a run ({type(exc).__name__}: {exc}); not tracking")
        yield _Null()
        return

    handle = _Live(mlflow)
    try:
        handle.params({**params, "drug": drug, "feature_set": feature_set})
        # The provenance that makes a metric reproducible, recorded as tags so it
        # is filterable in the UI rather than buried in a parameter blob.
        handle.tags({
            "pipeline": "mtb-resistotyper-ml",
            "nextflow_run": os.environ.get("NXF_UUID", ""),
            "nomad_alloc": os.environ.get("NOMAD_ALLOC_ID", ""),
            "container_tag": os.environ.get("MTB_CONTAINER_TAG", ""),
            "git_revision": os.environ.get("MTB_GIT_REVISION", ""),
        })
        yield handle
        handle.status("FINISHED")
    except Exception:
        # A training failure is a tracked outcome, not an absence of one.
        handle.status("FAILED")
        raise
    finally:
        with contextlib.suppress(Exception):
            mlflow.end_run()


class _Null:
    """Every method a no-op, so the caller never branches on availability."""

    def params(self, _: dict) -> None: ...
    def tags(self, _: dict) -> None: ...
    def metrics(self, _: dict, step: int | None = None) -> None: ...
    def child(self, name: str, params: dict, metrics: dict) -> None: ...
    def artifact(self, path: Path, subdir: str | None = None) -> None: ...
    def status(self, _: str) -> None: ...


class _Live:
    def __init__(self, mlflow) -> None:
        self._m = mlflow

    def params(self, d: dict) -> None:
        with contextlib.suppress(Exception):
            # MLflow rejects a parameter over 500 characters, and a base-model
            # list can exceed that. Truncate rather than lose the whole run.
            self._m.log_params({k: str(v)[:490] for k, v in d.items() if v is not None})

    def tags(self, d: dict) -> None:
        with contextlib.suppress(Exception):
            self._m.set_tags({k: str(v)[:490] for k, v in d.items() if v})

    def metrics(self, d: dict, step: int | None = None) -> None:
        with contextlib.suppress(Exception):
            self._m.log_metrics(
                {k: float(v) for k, v in d.items()
                 if isinstance(v, (int, float)) and v == v}, step=step)

    def child(self, name: str, params: dict, metrics: dict) -> None:
        """One nested run per base model or ensemble, so the leaderboard is browsable."""
        try:
            with self._m.start_run(run_name=name, nested=True):
                self._m.log_params({k: str(v)[:490] for k, v in params.items() if v is not None})
                self._m.log_metrics({k: float(v) for k, v in metrics.items()
                                     if isinstance(v, (int, float)) and v == v})
        except Exception as exc:
            _warn(f"child run {name!r} not logged ({type(exc).__name__})")

    def artifact(self, path: Path, subdir: str | None = None) -> None:
        try:
            if Path(path).is_dir():
                self._m.log_artifacts(str(path), artifact_path=subdir)
            elif Path(path).exists():
                self._m.log_artifact(str(path), artifact_path=subdir)
        except Exception as exc:
            _warn(f"artifact {path} not logged ({type(exc).__name__}: {exc})")

    def status(self, s: str) -> None:
        with contextlib.suppress(Exception):
            self._m.end_run(status=s)


def log_manifest(handle, manifest: dict, out_dir: Path) -> None:
    """Project a finished manifest onto the run: metrics, children, artefacts.

    Kept separate from the manifest's construction so that the two never
    disagree about what a run achieved: MLflow is fed the same object that is
    written to disk, not a parallel summary assembled beside it.
    """
    handle.metrics({
        "best_cv_auc": manifest["best_model"]["cv_auc"],
        "best_base_auc": manifest["best_base_auc"],
        "stacking_lift": manifest["stacking_lift_over_best_base"],
        "n_features": manifest["n_features"],
        "wall_time_seconds": manifest["wall_time_seconds"],
    })
    handle.params({
        "best_model": manifest["best_model"]["which"],
        "best_is_ensemble": manifest["best_model"]["is_ensemble"],
        "best_base_family": manifest["best_base_family"],
        "sort_metric": manifest["sort_metric"],
        "balance_classes": manifest["balance_classes"],
    })
    for b in manifest.get("layer2_base_models", []):
        handle.child(f"base:{b['family']}",
                     {"family": b["family"], "model_id": b["model_id"]},
                     {"cv_auc": b["cv_auc"], "cv_aucpr": b["cv_aucpr"],
                      "cv_logloss": b["cv_logloss"], "train_secs": b["train_secs"]})
    for e in manifest.get("layer3_ensembles", []):
        handle.child(f"ensemble:{e['name']}",
                     {"metalearner": e["metalearner"],
                      "n_base_models": e["n_base_models"],
                      "base_models": ",".join(e.get("base_models", [])),
                      "metalearner_weights": json.dumps(e.get("metalearner_weights"))},
                     {"cv_auc": e["cv_auc"], "cv_aucpr": e["cv_aucpr"],
                      "train_secs": e["train_secs"]})
    # The whole output directory: MOJO, leaderboard, GLM catalogue, manifest. The
    # MOJO is what makes the run redeployable rather than merely described.
    handle.artifact(out_dir)
