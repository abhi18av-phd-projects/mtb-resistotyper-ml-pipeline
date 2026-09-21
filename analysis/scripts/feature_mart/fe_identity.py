"""The identity of the feature-engineering technique a model was built with.

Two identifiers, because "which FE" and "which data" vary independently:

  fe_version  a semantic version, declared by a person, bumped deliberately.
              MAJOR: features are encoded or named differently, or the fold
                     scheme changes -- a model trained on one cannot read the
                     other's input.
              MINOR: the feature space changes but the encoding does not --
                     new covariates or steps, a new selection rule, a changed
                     default such as top_n.
              PATCH: output is unchanged on a fixed reference release.
  fe_id       "fe-" + a hash of the complete FE configuration: the version,
              the mart-builder settings, the determinant-gene map, the FE
              steps that ran and the selection settings. Measured, not
              declared. It deliberately EXCLUDES the database, so the same
              technique applied to CRyPTIC v2.1.2 and v3.4.0 carries the same
              fe_id -- which is what lets a cross-release comparison insist on
              equal FE. It also excludes the code hashes: those are recorded in
              fe_config, but hashing them would give two builds that produce
              identical features different ids after any edit, even a comment.
              Code changes that matter reach fe_id through FE_VERSION, which
              the lock forces to move.
  mart_id     "mart-" + a hash of fe_id and the database checksum: this FE on
              this database.

The version is kept honest by FE_LOCK.json: `--check` fails when the FE code has
changed but FE_VERSION has not, which is the failure that left the previous,
informal FE version (v1.0.2, carried as the mart filename label) unchanged
across every later change.

The configuration is assembled from blocks the PRODUCERS record -- build_mart
writes its settings, determinant-gene map and code hash into the mart sidecar;
causal.py writes its settings and code hash into the selection manifest -- so a
consumer never re-derives them from code it may not share. For artefacts built
before those blocks existed, `code_root` reconstructs them from the source tree
at the revision that built them.
"""

from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
import sys
from pathlib import Path

FE_VERSION = "1.0.4"

HERE = Path(__file__).resolve().parent               # analysis/scripts/feature_mart
SCRIPTS = HERE.parent                                 # analysis/scripts
LOCK = HERE / "FE_LOCK.json"

# Every file whose behaviour changes the features a model sees. A change to any
# of them is an FE change, whatever the parameters say.
FE_CODE = ["feature_mart/build_mart.py", "feature_mart/causal.py"]
FE_CODE_GLOBS = ["fe_steps/*.py"]

MART_SETTINGS = ("mutation_selector", "top_n_mutations", "max_carrier_frac",
                 "min_carrier_count", "gene_level_flags", "fold_strategy", "n_folds",
                 "random_state")
SELECTION_SETTINGS = ("candidate_level", "concordance_rule", "ebm_top_k",
                      "cate_n_estimators", "seeds")


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def code_sha256(root: Path = SCRIPTS, include_steps: bool = False) -> dict[str, str]:
    """sha256 of the FE source files under an analysis/scripts root.

    The identity records the two producers (mart builder, selection); the lock
    also covers every FE step, so a change to a step forces a version bump even
    though a step only enters fe_config when it actually runs.
    """
    files = [root / f for f in FE_CODE]
    if include_steps:
        for g in FE_CODE_GLOBS:
            files += sorted(root.glob(g))
    return {str(f.relative_to(root)): _sha(f) for f in files
            if f.exists() and f.name != "__init__.py"}


def _canon(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), default=str)


def _determinant_genes_from(root: Path) -> dict:
    """The determinant-gene map as defined by build_mart.py in a source tree."""
    spec = importlib.util.spec_from_file_location("_bm_for_fe", root / "feature_mart/build_mart.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return {k: sorted(v) for k, v in sorted(mod.DETERMINANT_GENES.items())}


def mart_block(meta: dict, code_root: Path | None = None) -> dict:
    """FE settings of the mart stage, from its sidecar (or reconstructed)."""
    if isinstance(meta.get("fe"), dict) and "mart" in meta["fe"]:
        return meta["fe"]["mart"]
    # Retro: sidecars written before the fe block carry the same settings under
    # build.carrier_rate_window and fold_assignment.
    b = dict(meta.get("build") or {})
    crw = b.get("carrier_rate_window") or {}
    fa = meta.get("fold_assignment") or {}
    settings = {
        "mutation_selector": crw.get("mutation_selector"),
        "top_n_mutations": b.get("top_n_mutations"),
        "max_carrier_frac": crw.get("max_carrier_frac"),
        "min_carrier_count": crw.get("min_carrier_count"),
        "gene_level_flags": crw.get("gene_level_flags"),
        "fold_strategy": fa.get("strategy") or crw.get("fold_strategy"),
        "n_folds": fa.get("n_splits"),
        "random_state": fa.get("random_state"),
    }
    if code_root is None:
        raise SystemExit("this mart predates the fe block; pass code_root, the "
                         "analysis/scripts tree at the revision that built it")
    return {"settings": settings,
            "determinant_genes": _determinant_genes_from(code_root)}


def selection_block(manifest: dict) -> dict:
    """FE settings of the concordance-selection stage (configured, not derived:
    the realised candidate-pool size depends on the data and is left out)."""
    if isinstance(manifest.get("fe"), dict) and "selection" in manifest["fe"]:
        return manifest["fe"]["selection"]
    return {"settings": {k: manifest.get(k) for k in SELECTION_SETTINGS}}


def steps_block(pre: str | None, fold: str | None) -> dict:
    """The optional FE steps that ran, each with its registry version."""
    def resolve(spec: str | None) -> list[dict]:
        names = [n.strip() for n in (spec or "").split(",") if n.strip()]
        if not names:
            return []
        sys.path.insert(0, str(SCRIPTS.parent.parent))
        from analysis.scripts.fe_steps import registry
        registry.load_builtin_steps()
        return [{"name": n, "version": getattr(registry.get(n), "version", "1.0.0")}
                for n in names]
    return {"pre_mart": resolve(pre), "per_fold": resolve(fold)}


def identity(mart_meta: dict, selection_manifest: dict | None = None,
             pre_steps: str | None = None, fold_steps: str | None = None,
             code_root: Path | None = None, fe_version: str | None = None) -> dict:
    """{fe_version, fe_id, fe_config} for one model's FE."""
    recorded = (mart_meta.get("fe") or {}).get("code_sha256")
    code = recorded if recorded else code_sha256(code_root or SCRIPTS)
    if selection_manifest and (selection_manifest.get("fe") or {}).get("code_sha256"):
        code = {**code, **selection_manifest["fe"]["code_sha256"]}
    version = fe_version or (mart_meta.get("fe") or {}).get("fe_version") or FE_VERSION
    hashed = {
        "fe_version": version,
        "mart": mart_block(mart_meta, code_root),
        "selection": selection_block(selection_manifest or {}),
        "steps": steps_block(pre_steps, fold_steps),
    }
    return {"fe_version": version,
            "fe_id": "fe-" + hashlib.sha256(_canon(hashed).encode()).hexdigest()[:12],
            "fe_config": {**hashed, "code_sha256": dict(sorted(code.items()))}}


def mart_id(fe_id: str, db_sha256: str | None) -> str:
    """This FE on this database."""
    return "mart-" + hashlib.sha256(f"{fe_id}|{db_sha256 or ''}".encode()).hexdigest()[:12]


def producer_block(stage: str, settings: dict, extra: dict | None = None) -> dict:
    """What a producer writes into its own manifest: version, settings, code hash.

    `stage` is "mart" or "selection". The code hash covers the producer's own
    source file, which is what the running process can actually vouch for.
    """
    src = {"mart": "feature_mart/build_mart.py", "selection": "feature_mart/causal.py"}[stage]
    body = {"settings": settings, **(extra or {})}
    return {"fe_version": FE_VERSION, stage: body,
            "code_sha256": {src: _sha(SCRIPTS / src)}}


def _check() -> int:
    now = code_sha256(include_steps=True)
    if not LOCK.exists():
        print(f"no {LOCK.name}; run with --lock to create it"); return 1
    lock = json.loads(LOCK.read_text())
    changed = sorted(k for k in set(now) | set(lock["code_sha256"])
                     if now.get(k) != lock["code_sha256"].get(k))
    if changed and lock["fe_version"] == FE_VERSION:
        print(f"FE code changed but FE_VERSION is still {FE_VERSION}:")
        for k in changed:
            print(f"  {k}")
        print("Bump FE_VERSION (see FE_CHANGELOG.md), then run --lock.")
        return 1
    if lock["fe_version"] != FE_VERSION and not changed:
        print(f"FE_VERSION moved to {FE_VERSION} with no FE code change; run --lock "
              "only if the change is intended")
    print(f"ok: FE {FE_VERSION}, {len(now)} files" + (" (lock stale: run --lock)" if changed else ""))
    return 0


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--check", action="store_true", help="fail if FE code changed without a version bump")
    p.add_argument("--lock", action="store_true", help="record the current FE code for FE_VERSION")
    a = p.parse_args()
    if a.lock:
        LOCK.write_text(json.dumps({"fe_version": FE_VERSION,
                                    "code_sha256": code_sha256(include_steps=True)},
                                   indent=2) + "\n")
        print(f"locked FE {FE_VERSION}"); return 0
    return _check()


if __name__ == "__main__":
    raise SystemExit(main())
