"""Package train_logistic's output as a self-contained, depositable model release.

`train_logistic` writes `<DRUG>_logistic.json` plus a run manifest: enough to
score with, not enough to deploy or to cite. This assembles one directory per
CRyPTIC release that can be zipped and uploaded as a single Zenodo record, and
that the serving layer installs unchanged:

    mtb-resistotyper-ml-models-<release>/
      README.md            what the release is, its operating range, how to use it
      CITATION.cff         how to cite it
      .zenodo.json         deposit metadata
      LICENSE              the models' licence
      MANIFEST.json        release, cohort, pipeline revision, images, per-drug summary
      SHA256SUMS           a checksum for every other file
      models/<DRUG>/       the artefact contract the runner consumes:
        model.json           scaler, coefficients, intercept, tuned C
        feature_schema.json  the exact input contract
        model_card.json      provenance and the measured operating range
      provenance/
        training_manifest.json   nested lineage-LOO detail, per drug
        marts/<DRUG>.json        cohort and mart provenance sidecars
        selection/<DRUG>.json    causal-concordance selection manifests

The feature marts themselves are not included: they are the pipeline's deposit
artefact (--create_duckdb), a separate record. What is here is everything needed
to USE a model, and everything needed to know how it was made.

Until this existed the shipped bundles were assembled by hand, which is how
"slim-2026.05" became a release label -- invented at packaging time rather than
read from the data. The release now comes from the marts the models were fitted
on, and a mart that names no release is refused.
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path

# The boundary the shipped bundles were tiered on: MXF 0.927 and EMB 0.909 are
# "usable", ETH 0.891 and KAN 0.883 "moderate". Stated here rather than left
# implicit, because tier_report.py uses 0.88 for the same word; recorded in
# every card so a reader can see which rule was applied.
TIERS = [(0.90, "usable"), (0.80, "moderate"), (0.0, "at chance")]

CREATORS = [{"name": "Sharma, Abhinav", "orcid": "0000-0002-6402-6993"}]


def _tier(auc: float | None) -> str:
    if auc is None:
        return "undeclared"
    for floor, name in TIERS:
        if auc >= floor:
            return name
    return "at chance"


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _mart_file(marts: Path, drug: str, suffix: str) -> Path | None:
    hits = sorted(glob.glob(str(marts / f"feature_mart_{drug}_*vFULL{suffix}")))
    return Path(hits[0]) if hits else None


def _schema(drug: str, features: list[str], reference: str,
            fills: dict[str, float] | None = None) -> dict:
    """The input contract, in the vocabulary the runner's vectorize() accepts.

    That vocabulary is closed: `mutation`, `lineage`, `coverage`, and anything
    else raises "unknown feature kind". This first wrote coverage columns as
    kind "covariate" with fill "mean" -- a schema that looked reasonable, passed
    through packaging silently, and made every prediction from a bundle carrying
    coverage features fail at serving time. A coverage feature names the
    instance covariate it reads (`source`) and the value to use when that
    covariate is absent (`fill`, the training-set median), exactly as the
    hand-assembled bundles already shipped do.
    """
    fills = fills or {}
    out = []
    for name in features:
        if name.startswith("cov__lineage_") and name != "cov__lineage_raw":
            out.append({"name": name, "kind": "lineage", "dtype": "int8",
                        "lineage": "lineage" + name.rsplit("_L", 1)[-1]})
        elif name.startswith("cov__"):
            if name not in fills:
                raise SystemExit(
                    f"{drug}: no training median for {name}; a coverage feature "
                    "without a numeric fill cannot be scored when the covariate "
                    "is absent from an instance")
            out.append({"name": name, "kind": "coverage", "dtype": "float32",
                        "source": name[len("cov__"):], "fill": fills[name]})
        else:
            entry = {"name": name, "kind": "mutation", "dtype": "int8"}
            gene, _, mutation = name[len("raw__"):].partition("_")
            if gene and mutation:
                entry |= {"gene": gene, "mutation": mutation}
            out.append(entry)
    return {"schema_version": "1.0.0", "drug": drug,
            "feature_set": "concordant+covariates", "n_features": len(out),
            "nomenclature": "GARC", "reference": reference, "features": out}


def _readme(release: str, rows: list[dict], missing: list[str],
            revision: str, images: str) -> str:
    table = "\n".join(
        f"| {r['drug']} | {r['auc']:.4f} | {r['tier']} | {r['n_features']} | "
        f"{r['n_isolates']:,} | {r['n_resistant']:,} |" for r in rows)
    caveat = ""
    if missing:
        caveat = (
            "\n**Features this release could not provide.** CRyPTIC "
            f"{release} does not publish {', '.join(missing)}. The marts carry "
            "no such column, so these models are fitted without the covariates "
            "derived from them -- a genuine difference in feature space from "
            "releases that publish them, not an omission of this packaging. "
            "Compare operating ranges across releases with that in mind.\n")
    return f"""# mtb-resistotyper-ml-models {release}

Per-drug L2-logistic resistance models for *Mycobacterium tuberculosis*, fitted
on **CRyPTIC {release}** and evaluated under nested lineage-leave-one-out
cross-validation.

**Research use only. Not a diagnostic.** These are Layer 2 models: intended to
fire only where the curated WHO catalogue returns Unknown or Fail, never to
override it. They have not been prospectively validated.

## Operating range

AUC under nested lineage-leave-one-out: trained on three of the four major
lineages, tested on the fourth, with the regularisation strength tuned inside
the training lineages only. Tier: usable >= 0.90, moderate >= 0.80.

| drug | AUC | tier | features | isolates | resistant |
|---|---|---|---|---|---|
{table}
{caveat}
Per-lineage AUC, stability and the tuned `C` are in each `model_card.json`.

## Layout

    models/<DRUG>/model.json           scaler, signed coefficients, intercept
    models/<DRUG>/feature_schema.json  the exact input contract
    models/<DRUG>/model_card.json      provenance and the operating range
    provenance/                        training, cohort and selection records
    MANIFEST.json                      the release in one machine-readable file
    SHA256SUMS                         verify with `shasum -a 256 -c SHA256SUMS`

## Using the models

The runner, `mtb-resistotyper-ml`, discovers any directory of bundles; point it
at `models/`. Scoring needs nothing beyond the standard library:

    p(R) = sigmoid( ((x - scaler_mean) / scaler_scale) . coef + intercept )

## Provenance

Trained by the mtb-resistotyper-ml pipeline (`{revision}`) in containers
{images}. The CRyPTIC compendium is CC-BY-4.0.
"""


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--models", type=Path, required=True, help="train_logistic --out")
    p.add_argument("--marts", type=Path, required=True)
    p.add_argument("--causal", type=Path, default=None,
                   help="causal output root: <DRUG>/manifest.json is copied to provenance")
    p.add_argument("--out", type=Path, required=True,
                   help="parent directory; the release directory is created inside it")
    p.add_argument("--reference", default="NC_000962.3")
    p.add_argument("--licence", default="EPL-2.0")
    p.add_argument("--licence-file", type=Path, default=None)
    p.add_argument("--pipeline-revision", default="unrecorded")
    p.add_argument("--images", default="unrecorded")
    a = p.parse_args()

    manifest = json.loads((a.models / "manifest.json").read_text())
    per_drug = manifest.get("per_drug", {})
    stamp = datetime.now(timezone.utc).isoformat(timespec="seconds")

    built = []
    for model_path in sorted(a.models.glob("*_logistic.json")):
        model = json.loads(model_path.read_text())
        drug = model["drug"]
        meta_path = _mart_file(a.marts, drug, ".metadata.json")
        meta = json.loads(meta_path.read_text()) if meta_path else {}
        if not meta.get("cryptic_version"):
            raise SystemExit(
                f"{drug}: the mart metadata names no cryptic_version, so the "
                "release cannot state which data trained it.")
        built.append((drug, model, per_drug.get(drug, {}), meta, meta_path))

    releases = {m["cryptic_version"] for _, _, _, m, _ in built}
    if len(releases) != 1:
        raise SystemExit(f"one release per package; these marts span {sorted(releases)}")
    release = releases.pop()

    root = a.out / f"mtb-resistotyper-ml-models-{release}"
    if root.exists():
        shutil.rmtree(root)
    (root / "models").mkdir(parents=True)
    (root / "provenance" / "marts").mkdir(parents=True)
    (root / "provenance" / "selection").mkdir(parents=True)

    rows, missing = [], set()
    for drug, model, summary, meta, meta_path in built:
        d = root / "models" / drug
        d.mkdir()
        (d / "model.json").write_text(json.dumps(model, indent=2) + "\n")

        # Training-set medians for the coverage features, from the mart the
        # model was fitted on: an instance that omits a covariate is scored at
        # the centre of the training distribution rather than at zero.
        cov = [f for f in model["features"]
               if f.startswith("cov__") and not f.startswith("cov__lineage_")]
        fills = {}
        if cov:
            import pandas as pd
            mart = _mart_file(a.marts, drug, ".parquet")
            if mart is None:
                raise SystemExit(f"{drug}: no mart to take coverage medians from")
            df = pd.read_parquet(mart, columns=cov)
            fills = {c: float(df[c].median()) for c in cov}
        (d / "feature_schema.json").write_text(json.dumps(
            _schema(drug, model["features"], a.reference, fills), indent=2) + "\n")

        auc = summary.get("honest_nested_loo_auc")
        missing |= set(meta.get("missing_genome_covariates") or [])
        missing |= set(meta.get("missing_mutation_columns") or [])
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
                "C_per_outer_fold": summary.get("C_per_outer_fold"),
            },
            "data": {
                "source": "CRyPTIC compendium",
                "release": release,
                "n_isolates": meta.get("n_samples"),
                "n_resistant": meta.get("n_R"),
                "pct_resistant": meta.get("pct_R"),
                "missing_genome_covariates": meta.get("missing_genome_covariates"),
                "missing_mutation_columns": meta.get("missing_mutation_columns"),
                "licence": "CC-BY-4.0",
            },
            "provenance": {
                "packaged_utc": stamp,
                "mart": meta.get("mart_version"),
                "pipeline_revision": a.pipeline_revision,
                "images": a.images,
                "packaged_by": "analysis.scripts.feature_mart.package_bundles",
            },
            "intended_use": ("Layer 2. Fires only where the curated WHO catalogue "
                             "returns Unknown or Fail, and never overrides a "
                             "catalogue call. Research use only: not prospectively "
                             "validated, no regulatory clearance."),
            "licence": a.licence,
        }
        (d / "model_card.json").write_text(json.dumps(card, indent=2) + "\n")

        if meta_path:
            shutil.copy2(meta_path, root / "provenance" / "marts" / f"{drug}.json")
        if a.causal and (a.causal / drug / "manifest.json").exists():
            shutil.copy2(a.causal / drug / "manifest.json",
                         root / "provenance" / "selection" / f"{drug}.json")
        rows.append({"drug": drug, "auc": auc or float("nan"), "tier": card["operating_range"]["tier"],
                     "n_features": card["training"]["n_features"],
                     "n_isolates": meta.get("n_samples") or 0,
                     "n_resistant": meta.get("n_R") or 0})

    shutil.copy2(a.models / "manifest.json", root / "provenance" / "training_manifest.json")
    rows.sort(key=lambda r: -r["auc"])
    missing_l = sorted(missing)

    (root / "MANIFEST.json").write_text(json.dumps({
        "release": release,
        "packaged_utc": stamp,
        "pipeline_revision": a.pipeline_revision,
        "images": a.images,
        "evaluation": "nested lineage-leave-one-out cross-validation",
        "tier_thresholds": {name: floor for floor, name in TIERS},
        "missing_features_in_release": missing_l,
        "models": rows,
    }, indent=2) + "\n")
    (root / "README.md").write_text(_readme(release, rows, missing_l,
                                            a.pipeline_revision, a.images))
    if a.licence_file and a.licence_file.exists():
        shutil.copy2(a.licence_file, root / "LICENSE")
    title = (f"mtb-resistotyper-ml-models {release}: per-drug resistance models "
             f"for M. tuberculosis trained on CRyPTIC {release}")
    (root / ".zenodo.json").write_text(json.dumps({
        "title": title, "upload_type": "dataset", "version": release,
        "license": a.licence, "creators": CREATORS,
        "keywords": ["antimicrobial resistance", "tuberculosis", "machine learning",
                     "model distribution", "lineage-aware evaluation", "CRyPTIC"],
        "description": (f"Per-drug L2-logistic resistance models for Mycobacterium "
                        f"tuberculosis trained on CRyPTIC {release} and evaluated under "
                        "nested lineage-leave-one-out cross-validation. Each bundle "
                        "carries model.json, feature_schema.json and model_card.json "
                        "with the measured operating range. Layer 2 artefacts: they "
                        "fire only where the WHO catalogue returns Unknown or Fail. "
                        "Research use only; not prospectively validated."),
    }, indent=2) + "\n")
    (root / "CITATION.cff").write_text(
        "cff-version: 1.2.0\n"
        'message: "If you use these models, please cite the accompanying article."\n'
        f'title: "{title}"\n'
        f'version: "{release}"\n'
        "authors:\n  - family-names: Sharma\n    given-names: Abhinav\n"
        '    orcid: "https://orcid.org/0000-0002-6402-6993"\n'
        f"license: {a.licence}\ntype: dataset\n")

    files = sorted(f for f in root.rglob("*") if f.is_file() and f.name != "SHA256SUMS")
    (root / "SHA256SUMS").write_text("".join(
        f"{_sha256(f)}  {f.relative_to(root)}\n" for f in files))

    for r in rows:
        print(f"{r['drug']:<5} {release:<12} auc={r['auc']:.4f} tier={r['tier']}")
    print(f"wrote {len(rows)} models -> {root}")


if __name__ == "__main__":
    main()
