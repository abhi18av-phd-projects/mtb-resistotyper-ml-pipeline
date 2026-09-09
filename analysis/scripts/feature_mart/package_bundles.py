"""Turn train_logistic's output into deployable model bundles.

`train_logistic` writes `<DRUG>_logistic.json` plus a run manifest: enough to
score with, not enough to deploy. A bundle is the artefact contract the runner
consumes -- `model.json`, `feature_schema.json`, `model_card.json` -- and the
card is the part that matters, because it carries the measured operating range
inside the artefact so a model cannot be loaded without also receiving the
statement of where it works.

Until now that packaging was done by hand, which is why the shipped bundles
carry `slim-2026.05` as their release: a label invented at packaging time
rather than read from the data. Here the release is read from the mart metadata
the model was actually trained on, in CRyPTIC's own convention.

Usage:
    python -m analysis.scripts.feature_mart.package_bundles \
        --models logistic_models --marts marts --out bundles
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path

# The boundary the shipped bundles were tiered on: MXF 0.927 and EMB 0.909 are
# "usable", ETH 0.891 and KAN 0.883 "moderate". Stated here rather than left
# implicit, because analysis/scripts/feature_mart/tier_report.py uses 0.88 for
# the same word and the two disagree between 0.88 and 0.90 -- a model can be
# "usable" in one artefact and "moderate" in the other. Recorded in every card
# it produces so the reader can see which rule was applied.
TIERS = [(0.90, "usable"), (0.80, "moderate"), (0.0, "at chance")]

# raw__katG_S315T -> (katG, S315T). Gene names carry underscores (PE_PGRS7),
# so the split is on the FIRST underscore after the prefix only when the rest
# parses as a mutation; otherwise the whole tail is the mutation and the gene
# is whatever preceded it in the mart's own naming.
_FEATURE = re.compile(r"^raw__(?P<gene>.+?)_(?P<mutation>[A-Za-z0-9*\-]+)$")


def _tier(auc: float | None) -> str:
    if auc is None:
        return "undeclared"
    for floor, name in TIERS:
        if auc >= floor:
            return name
    return "at chance"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mart_metadata(marts: Path, drug: str) -> dict:
    for pattern in (f"feature_mart_{drug}_*vFULL.metadata.json",
                    f"feature_mart_{drug}_vFULL.metadata.json"):
        hits = sorted(glob.glob(str(marts / pattern)))
        if hits:
            return json.loads(Path(hits[0]).read_text())
    return {}


def _schema(drug: str, features: list[str], reference: str) -> dict:
    out = []
    for name in features:
        if not name.startswith("raw__"):
            # Covariates and engineered columns are inputs too; naming them
            # "mutation" would make the contract lie about what they are.
            out.append({"name": name, "kind": "covariate", "dtype": "float32",
                        "fill": "mean"})
            continue
        m = _FEATURE.match(name)
        entry = {"name": name, "kind": "mutation", "dtype": "int8", "fill": 0}
        if m:
            entry |= {"gene": m.group("gene"), "mutation": m.group("mutation")}
        out.append(entry)
    return {"schema_version": "1.0.0", "drug": drug, "feature_set": "concordant+covariates",
            "n_features": len(out), "nomenclature": "GARC", "reference": reference,
            "features": out}


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--models", type=Path, required=True,
                   help="train_logistic --out directory")
    p.add_argument("--marts", type=Path, required=True)
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--reference", default="NC_000962.3")
    p.add_argument("--licence", default="EPL-2.0")
    a = p.parse_args()

    manifest = json.loads((a.models / "manifest.json").read_text())
    per_drug = manifest.get("per_drug", {})
    a.out.mkdir(parents=True, exist_ok=True)
    written = []

    for model_path in sorted(a.models.glob("*_logistic.json")):
        model = json.loads(model_path.read_text())
        drug = model["drug"]
        summary = per_drug.get(drug, {})
        meta = _mart_metadata(a.marts, drug)
        release = meta.get("cryptic_version")
        if not release:
            raise SystemExit(
                f"{drug}: the mart metadata names no cryptic_version, so the "
                "bundle cannot state which release trained it. Rebuild the mart "
                "with a build_mart that reads the release from the database.")

        d = a.out / drug
        d.mkdir(exist_ok=True)
        (d / "model.json").write_text(json.dumps(model, indent=2) + "\n")
        (d / "feature_schema.json").write_text(json.dumps(
            _schema(drug, model["features"], a.reference), indent=2) + "\n")

        auc = summary.get("honest_nested_loo_auc")
        card = {
            "schema_version": "1.0.0",
            "drug": drug,
            "model": model.get("model", "standardscaler+l2_logistic"),
            "artefact": {
                "file": "model.json",
                "sha256": _sha256(d / "model.json"),
                "bytes": (d / "model.json").stat().st_size,
                "inference": ("stdlib Python; sigmoid((x - scaler_mean)/scaler_scale "
                              ". coef + intercept)"),
            },
            "operating_range": {
                "primary_metric": "nested lineage-leave-one-out cross-validation",
                "auc": auc,
                "tier": _tier(auc),
                "tier_thresholds": {name: floor for floor, name in TIERS},
                "per_lineage_auc": summary.get("per_lineage"),
                "stability_std": summary.get("stability_std"),
                "comparator_auc_lightgbm": summary.get("baseline_lightgbm"),
                "caveat": ("Measured on the CRyPTIC compendium. A deploying "
                           "population whose lineage or geographic mix differs from "
                           "that cohort may see lower performance; the per-lineage "
                           "figures are given so that this can be judged rather "
                           "than assumed."),
            },
            "training": {
                "n_features": summary.get("n_features", len(model["features"])),
                "n_concordant": summary.get("n_concordant"),
                "tuned_C": model.get("C"),
            },
            "data": {
                "source": "CRyPTIC compendium",
                "release": release,
                "n_isolates": meta.get("n_samples"),
                "n_resistant": meta.get("n_R"),
                "pct_resistant": meta.get("pct_R"),
                "licence": "CC-BY-4.0",
            },
            "provenance": {
                "packaged_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "mart": meta.get("mart_version"),
                "packaged_by": "analysis.scripts.feature_mart.package_bundles",
            },
            "intended_use": ("Layer 2. Fires only where the curated WHO catalogue "
                             "returns Unknown or Fail, and never overrides a "
                             "catalogue call. Research use only: not prospectively "
                             "validated, no regulatory clearance."),
            "licence": a.licence,
        }
        (d / "model_card.json").write_text(json.dumps(card, indent=2) + "\n")
        written.append((drug, release, auc, card["operating_range"]["tier"]))

    for drug, release, auc, tier in written:
        print(f"{drug:<5} {release:<14} auc={auc} tier={tier}")
    print(f"wrote {len(written)} bundles to {a.out}")


if __name__ == "__main__":
    main()
