# mtb-resistotyper-ml-pipeline

Feature engineering and model training for *Mycobacterium tuberculosis* drug-resistance
prediction, as two Nextflow pipelines.

```bash
# build the database from a CRyPTIC release (rarely)
nextflow run ./cryptic-db --cryptic_src /path/to/cryptic-tables-v3.4.0

# then engineer features and train (often)
nextflow run . -params-file experiments/baseline.yml
```

## The rule the pipeline exists to enforce

Every feature-engineering step declares whether it reads the phenotype label, and that
property alone decides where it may run.

| | runs | valid for |
|---|---|---|
| **label-free** — region filters, class filters, annotation joins | once, before the mart | every fold |
| **label-using** — association ranking, concordance selection | inside the training fold | that fold only |

Run a label-using step before the mart and it picks the column universe while looking at the
lineage held out for testing. `registry.assert_placement` refuses it, so a misplacement is a
configuration error rather than a mart that looks fine and is not.

That is not hypothetical. Concordance selection previously fitted on the whole mart —
`causal.py` said so in its own docstring — and every lineage-leave-one-out number in the
project was computed with features chosen while looking at the held-out lineage. Doing it
properly is 11 drugs × 4 lineages = 44 independent runs of an EBM, a causal forest and a
refutation pass: unaffordable in a loop on one machine, unremarkable distributed. That is what
this pipeline is for.

## An experiment is a params file

```bash
nextflow run . -params-file experiments/drop-pe-ppe.yml
```

Committing the file makes an experiment a reviewable artefact rather than a shell command
somebody remembers. Each varies exactly one thing against `baseline.yml`, owns its `outdir`,
and states its prediction in the header — a prediction written before the run is evidence, one
written after is a description. See `experiments/README.md`.

Sweeping several FE configurations against identical downstream folds:

```bash
nextflow run . --fe_configs conf/fe-sweep.example.csv
```

## Registered FE steps

```bash
python -m analysis.scripts.fe_steps.run_step --list --mart x --out y --steps x --phase pre-mart
```

| step | layer | phase | source |
|---|---|---|---|
| `region_filter` | universe | pre-mart | PE/PPE, ESX, IS/phage — measured population-structure proxies |
| `class_filter` | universe | pre-mart | mutation classes; synonymous is the measured lineage proxy |
| `callability_filter` | universe | pre-mart | Marin 2022 EBR/RLC low-callability regions |
| `rank_association` | rank | per-fold | chi² or mutual information on the binary label |
| `rank_mic` | rank | per-fold | Mann-Whitney effect on the graded MIC — Pei 2024 |
| `select_cross_model` | select | per-fold | random forest ∩ MLP top-N — Ghosh 2026 |

Adding one means writing a `run()` and calling `register()`. Nothing in `main.nf` changes,
because the workflow only ever passes a string and a phase.

## Storage on the abc cluster

Paths follow ADR-0060. **A project is a prefix in its group's bucket, never a bucket of its
own** — the shared pipeline cache and the workbench home are both keyed on the *group* bucket,
so a project bucket silently opts out of both.

```
s3://${MTB_GROUP_BUCKET}/
  references/cryptic/cryptic-slim.duckdb        curated, read-only, shared
  pipelines/mtb-resistotyper-ml/workdir         group-shared cache — runner-only
  users/abhi/results/<experiment>               visible as ~/abc/me/<group>/results/
```

`workDir` sits in the **group-shared** area on purpose. The deterministic resume session UUID
keys off that path, so a colleague re-running this pipeline reuses the work a previous run did
rather than repeating it. Put it under a user or a project and that sharing quietly stops.

Override any of it:

```bash
nextflow run . -profile nomad   --group_bucket su-other-group --user_prefix users/someone-else
```

MinIO credentials come from the environment, never the repository, and the **API port is
`:9000`** — not the console `:9001`:

```bash
export AWS_ACCESS_KEY_ID=...  AWS_SECRET_ACCESS_KEY=...  ABC_MINIO_ENDPOINT=...
nextflow run . -profile nomad -params-file experiments/baseline.yml
```

## Profiles

| profile | what it does |
|---|---|
| `standard` | local, no containers, the virtual environments on your machine |
| `containers` | local Docker, using the three images |
| `nomad` | the abc cluster, containers required |
| `test` | one drug, one lineage, short H2O runtime |

## Containers

Three images, because the stacks are mutually exclusive — `econml`/`dowhy` pin numpy 1.x while
this project's H2O environment is numpy 2.x. See `containers/README.md`. Each image carries
`analysis/` at `/opt/mtb`, which is what lets this repository stay pure Nextflow: a Nomad task
needs no shared filesystem and no checkout to find its scripts.

**The consequence to remember:** changing an analysis script means rebuilding the image. The
workflow is versioned by its git tag, the analysis code by the image tag, and the two are
joined in `nextflow.config`.

## Licence

EPL-2.0.
