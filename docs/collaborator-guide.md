# Working on the features and the models

For a collaborator whose remit is **feature engineering and machine learning**, not the
infrastructure underneath. It assumes you have never run this pipeline and do not want to
learn Nextflow, Nomad or MinIO in order to ask a scientific question.

The split it is written against: **the pipeline dimension is somebody else's** — does a thing
run reproducibly, at scale, on the cluster — and **the science dimension is yours**: which
features the model sees, which learner sees them, which hyperparameters, which metric, and
what the resulting number means. Most stages need both. The design of that split is
`design/14-team-roadmap-abhi-ana.md` in the research repository; this document is how you
actually drive the thing.

---

## 1. What is yours

| | yours | not yours |
|---|---|---|
| **Features** | which techniques run, in which order, with which parameters; adding a new one | the mart builder's internals, the parquet layout, the staging |
| **Models** | learner, ensemble shape, hyperparameters, runtime budget, sort metric, class balancing | the JVM, the container, how a task gets scheduled |
| **Evaluation** | which CV regime is primary, what a per-lineage delta means, which drugs are reportable | the fold machinery, the deposit format |
| **Interpretation** | everything | — |

The practical consequence: **you should never need to edit `main.nf`, `nextflow.config` or
anything in `modules/`.** If a question you want to ask requires editing one of those, that is
a gap in the pipeline rather than work for you — say so and it gets closed. Everything in this
guide is a params file, a step file, or a flag.

---

## 2. Access you need

Ask for these explicitly; none of them is automatic.

1. **Read access to two GitHub repositories** — `mtb-resistotyper-ml-pipeline` (this one, the FE
   and training pipelines) and `mtb-resistotyper-ml` (the runner and the deployment surfaces).
   Both are private.
2. **A `site.env`** — copy `site.env.example` and have the real values filled in. It holds the
   bucket, the object-store endpoint, the Nomad namespace and node pools, and the abc CLI
   context. It is gitignored on purpose: everything in it describes *where* the pipeline runs,
   never *what* it computes.
3. **An abc cluster account and context**, if you intend to run anything at full scale. Local
   Docker is enough for everything in §4 and §5 and for a `-profile test` run; it is not enough
   for a full campaign (11 drugs × 4 lineages of EBM plus a causal forest plus a refutation
   pass).
4. **Object-store credentials** for the group bucket, as environment variables. Never in a file
   in the repository.
5. **A database** — see §6. You do not have to build one.

Two things that are *not* access problems and look like them: a worker job that is accepted and
never placed (wrong node pool — `MTB_WORKER_POOL` in `site.env`), and a database fetched at
0.4 MB/s (the head was handed the public edge instead of the internal object-store address —
`MTB_S3_ENDPOINT`). `scripts/launch.sh` carries the corrections for both; its header explains
each one.

---

## 3. The one rule

**Every feature-engineering step declares whether it reads the phenotype label, and that
property alone decides where it may run.**

```
  ┌────────── label-free · runs ONCE, before the mart ──────────┐
  │  L0 universe      L1 annotate     L2 derive     L3 pivot     │
  │  shrink the       attach priors   add columns   → THE MART   │
  │  candidate set                                  fold-agnostic│
  └──────────────────────────────────────────────────────────────┘
  ══════════════════════ LABEL BOUNDARY ══════════════════════
  ┌──────── label-using · runs INSIDE each training fold ───────┐
  │  L4 rank          L5 select               L6 train          │
  │  narrow columns   concordance & friends   H2O               │
  └─────────────────────────────────────────────────────────────┘
```

L0–L3 are the expensive part and run once per drug and FE configuration. L4–L5 are the
fan-out. The mart is the substrate; what varies by fold is which subset of it the model sees.

Run a label-using step before the mart and the column universe is chosen while looking at the
lineage being held out for testing. `registry.assert_placement` refuses that, so a
misplacement is a configuration error at startup rather than a mart that looks fine and is not.

This is not hypothetical and it is worth knowing before you trust any older number from this
project: concordance selection previously fitted on the whole mart, and every
lineage-leave-one-out figure produced before the fix was computed with features chosen while
looking at the held-out lineage. Doing it properly is 11 drugs × 4 lineages = 44 independent
runs of an EBM, a causal forest and a refutation pass. That is what the pipeline exists for.

A mart carries **`label_safe`** in its note, true only when every step in its chain was
label-free. A mart that fails it is valid for a deployment refit and invalid for evaluation,
and says so about itself rather than relying on anyone's memory.

---

## 4. Day one

### A run that does nothing

Every process has a stub, so the whole wiring can be exercised in seconds with no container,
no database and no cluster:

```bash
nextflow run . -stub-run --db /dev/null --drugs RIF
```

Do this after any change you make. What it catches is the class of fault that otherwise costs a
full campaign to find — a tuple published to a path that does not exist, two campaign arms fed
into one channel, manifests paired by directory instead of by fold. All three shipped, and all
three were only reachable by running the real work first.

### A run that does something

```bash
nextflow run . -profile containers,test --db /path/to/cryptic.duckdb
```

`test` is one drug, one lineage, a 25-candidate pool and a 60-second H2O budget: a couple of
minutes, and it still exercises every process. Run it before you run anything real.

**Use `-profile containers`, not `standard`.** The `standard` profile runs against virtual
environments on the host at paths that only resolve when this repository is checked out *inside*
the research monorepo as `pipeline/`. In a standalone clone they do not exist and the first task
fails with an interpreter that is not there. The containers carry `analysis/` at `/opt/mtb`, so
they need no checkout and no shared filesystem.

---

## 5. Changing feature engineering

### The knobs that already exist

These are the mart-build parameters. Set them in a params file (§5.4), not on the command line,
once you care about the answer.

| param | default | what it decides |
|---|---|---|
| `fe_mutation_selector` | `top_n` | `top_n` \| `chi2` \| `mi` \| `all` |
| `fe_top_n` | 1000 | how many mutation columns survive ranking |
| `fe_max_carrier_frac` | 0.95 | upper carrier-rate bound — above this a mutation is a lineage marker, not a resistance variant |
| `fe_min_carrier_count` | 10 | lower bound; below this nothing can be estimated |
| `fe_gene_level_flags` | true | add gene-level aggregate columns |
| `fe_fold_strategy` | `stratified` | `stratified` \| `lineage` \| `country` — use `lineage` for the honest regime |
| `fe_n_folds` | 5 | more folds for low-prevalence drugs is an open question |
| `candidate_level` | `mutation` | `mutation` \| `gene` |
| `candidate_pool_size` | null (300) | candidates entering the concordance vote; this dominates selection wall time |

The carrier-rate window is the one to sanity-check first on a drug you know well. Its job is to
stop the top-N picking up lineage markers instead of resistance variants, and it has only been
checked carefully on a few drugs.

### The step chains

Techniques beyond the built-in mart parameters are **registered steps**, composed as ordered
chains, one chain each side of the label boundary:

```bash
nextflow run . --db ... \
  --fe_pre_steps  'region_filter:families=pe_ppe,esx' \
  --fe_fold_steps 'rank_association:method=chi2;top_n=1000'
```

A step may carry parameters as `name:k=v;k=v`. The currently registered set:

| step | layer | phase | what it does |
|---|---|---|---|
| `region_filter` | universe | pre-mart | drops PE/PPE, ESX, IS/phage — measured population-structure proxies |
| `class_filter` | universe | pre-mart | keeps or drops mutation classes; synonymous is the measured lineage proxy |
| `callability_filter` | universe | pre-mart | drops Marin 2022 EBR/RLC low-callability regions |
| `rank_association` | rank | per-fold | chi² or mutual information on the binary label |
| `rank_mic` | rank | per-fold | Mann-Whitney effect on the graded MIC |
| `select_cross_model` | select | per-fold | random forest ∩ MLP top-N |

To print the live registry rather than trusting this table:

```bash
python -m analysis.scripts.fe_steps.run_step --list --mart x --out y --steps x --phase pre-mart
```

The four dummy arguments are a wart — `--list` prints and exits, but argparse requires them
anyway. Ignore them.

### Adding a step

One contract for every technique, which is what keeps this from being pipeline work:

```
mart.parquet + gene_annotations.parquet + params
    →  mart.parquet + a note recording what the step did and whether it was label-free
```

Write a `run()` and call `register()`. Nothing in `main.nf` changes, because the workflow only
ever passes a string and a phase. `analysis/scripts/fe_steps/class_filter.py` is the shortest
complete example; copy its shape.

```python
def run(mart: pd.DataFrame, *, your_param: str = "default", **_) -> tuple[pd.DataFrame, dict]:
    ...
    return kept, {"step": "your_step", ...}   # the note; say what you did, in numbers

register(Step(
    name="your_step",
    layer="universe",        # universe | annotate | derive | rank | select
    label_free=True,         # the declaration the boundary is enforced on
    summary="one line, for the registry listing",
    run=run,
))
```

Three things the note should carry, because they are what makes a result readable six months
later: how many features went in, how many came out, and the fraction dropped. Add the step
name to `_BUILTIN_STEPS` in `registry.py` so it is imported at registry load.

Note that steps are imported **one at a time** and a step whose dependency is missing is
reported on stderr and left out of the registry rather than ending the run. A dependency
missing for a step nobody asked for is a reduced menu, not a failure.

### An experiment is a params file

```bash
nextflow run . -params-file experiments/drop-pe-ppe.yml
```

Committing the file is what makes an experiment a reviewable artefact rather than a shell
command somebody remembers. `experiments/baseline.yml` is the thing every other file copies.
The rules that make two runs comparable, from `experiments/README.md`:

1. **Vary one thing.** Everything else is copied verbatim from the baseline. A file that
   changes the FE *and* the drug list answers no question.
2. **Name the question, not the settings.** `does-callability-reach-beyond-pe-ppe`, not
   `callability-v2`.
3. **Give every experiment its own `outdir`.** Two result sets in one directory cannot be
   diffed.
4. **State the question and the prediction in the header comment.** A prediction written before
   the run is evidence; one written after is a description.
5. **Pin the database.** Two runs against different databases are not a comparison.

To put several FE configurations against identical downstream folds in one run:

```bash
nextflow run . --fe_configs conf/fe-sweep.example.csv
```

### What you must bump when you change FE code

Every mart, deposit and model records which feature engineering produced it, and the scheme
only works if the declaration moves when the code does.

```bash
python -m analysis.scripts.feature_mart.fe_identity --check   # fails if FE code moved and FE_VERSION did not
# bump FE_VERSION, add an entry to analysis/scripts/feature_mart/FE_CHANGELOG.md, then:
python -m analysis.scripts.feature_mart.fe_identity --lock
```

MAJOR when features are encoded or named differently or the fold scheme changes; MINOR when the
feature space changes but the encoding does not; PATCH when output is unchanged on a fixed
reference release. `docs/fe-versioning.md` has the full scheme and what each identifier
deliberately excludes.

**The analysis code ships inside the container image, so changing a step means rebuilding the
image** before it can run on the cluster. That is a tag push, not a manual build — see the
README's *Containers* section. Locally, `-profile containers` with a rebuilt local image, or
`-profile standard` inside the monorepo layout, both pick up an edit without a release.

---

## 6. The data

Two CRyPTIC compendium databases, and the choice is a real one:

| | size | has `variants` | for |
|---|---|---|---|
| `cryptic-slim.duckdb` | ~2.8 GB | no | everything below; mutation-level, laptop-feasible |
| `cryptic.duckdb` | ~19 GB | yes | genome-coordinate features; a deliberate later phase |

Two release editions exist, **v2.1.2** and **v3.4.0**, and they are not interchangeable: the
RIF mart on v3.4.0 holds 14,490 isolates and 1,817 columns against v2.1.2's 14,602 and 2,028.
Some columns exist only in one — the release-tag column the cohort extraction needs is published
only by the v3.4.0 build.

You do not build a database. They are curated, read-only and shared, at
`s3://<group-bucket>/references/cryptic/`. On the cluster the `nomad` profile resolves `db`
there from `db_name` alone, so a committed experiment names a file and not a bucket. Locally,
ask for a copy of the slim database and pass `--db` a path.

### Starting from a deposit instead

`--from_duckdb <file>` skips the database assembly, feature engineering and causal selection
entirely, unpacks their results back into the layout the training stages expect, and starts at
training:

```bash
nextflow run . -profile containers --from_duckdb mtb-resistotyper-fe.duckdb
```

This is the fastest path to **working on models without touching features at all** — no cluster,
no database, and none of the twelve hours that produced the deposit. Two things to know:

- The deposit holds **both campaign arms and they are not interchangeable.** Units with a
  held-out lineage carry a selection made inside that training fold and feed evaluation only;
  the `none` unit carries a selection refit on all data and feeds the shipped model only.
  Crossing them produces deployment models fitted on fold-restricted features and a "held-out"
  score computed on rows whose features already saw them.
- A deposit rerun does **not** redo feature engineering, so this run's `fe_*` settings do not
  describe the mart. Only `--fe_name` is used, as a grouping label.

---

## 7. Changing the model

Training is H2O: AutoML, then base models, then stacked ensembles, with the JVM booted inside
the task and torn down with it. The parameters the pipeline exposes:

| param | default | what it decides |
|---|---|---|
| `h2o_entrypoint` | `train_h2o_stacked` | which training script runs; `train_h2o` is the older single-AutoML one |
| `h2o_feature_set` | `concordant` | `denovo` \| `concordant` \| `all` — which arm's features the learner sees |
| `h2o_runtime_secs` | 900 | the AutoML budget, per drug per fold |
| `h2o_sort_metric` | `AUC` | **see below** |
| `h2o_balance_classes` | false | class weighting |

**The sort metric is the decision worth making deliberately.** AUC weights both classes
equally, which is the wrong instinct where resistance is rare — BDQ, LZD and CFZ sit at roughly
1–5% — because a classifier that never predicts R still scores well on it. AUCPR is the honest
choice there. The default is AUC because it is comparable across drugs, not because it is right
for every drug.

The heap is derived from the task's own memory allocation rather than set independently, so the
JVM cannot overcommit what the scheduler promised it. If you need a bigger budget, that is a
resource label in `nextflow.config` and therefore a conversation, not a flag.

Anything beyond these — the meta-learner choice, the blender architecture, calibration, the
base-learner sweep, `exclude_algos` — lives in
`analysis/scripts/feature_mart/train_h2o_stacked.py`. That file is yours. Editing it means
bumping the H2O image tag before it runs on the cluster, for the same reason as FE.

### Tracking

Off by default, because a pipeline that fails because a dashboard is unreachable has its
priorities backwards. Point it at a server and runs record themselves as they happen:

```bash
nextflow run . --mlflow_uri https://mlflow.example.org   # or set MTB_MLFLOW_URI
```

Each run is tagged with the container tag, the git revision and the session UUID, so a number in
the UI traces back to the run and the image that made it. `manifest.json` in the output stays
the source of truth either way.

---

## 8. Reading the output

```
<outdir>/
  database/                        integrity checks
  cohort/                          the labelled substrate
  marts/<drug>/<fe-name>/          the mart and its sidecar — the reproducibility unit
  selection/<drug>/held-out-<L>/   concordance vote and the per-fold selection manifest
  evaluation/<drug>/               cross-validation manifests
  models/<drug>/held-out-<L>/      H2O artefacts and h2o_manifest.json
  tiers/                           the per-drug tier report — start here
  publish/                         the DuckDB deposit, if --create_duckdb was passed
  reports/                         Nextflow's own report, timeline, trace and DAG
```

`held-out-none` is the deployment refit, not a fold. The four `held-out-lineageN` directories
are the honest evaluation. Do not average them with the first.

Start at `tiers/`. It is the per-drug summary the rest of the tree supports.

---

## 9. Running on the cluster

```bash
scripts/launch.sh experiments/<name>.yml
```

Not a bare `abc pipeline run`. The script carries four corrections this deployment needs — the
nf-nomad plugin pin, the head and worker node pools, and the internal object-store endpoint —
and each of them, omitted, fails in a way that looks like something else. The header says how.

Two things to expect. **Submission is bounded to three concurrent tasks** on purpose: left
unbounded, Nextflow submitted five 6-CPU selections at once onto a node shared with eighteen
platform services, and the placement failure arrives dressed as a task failure, so Nextflow
retried and each retry added another blocked reservation to a queue that was already full. And
**the work directory is group-shared**, so the resume session key is shared too: a colleague
re-running the same pipeline reuses what your run did rather than repeating it.

A full campaign is on the order of twelve hours. Run `-profile test` first. Always.

---

## 10. When it breaks

| symptom | cause |
|---|---|
| a step is "unknown" but you can see the file | it is not in `_BUILTIN_STEPS`, or its import failed — check stderr for `[fe_steps] step ... unavailable` |
| configuration refuses a step you placed | it reads the label and you put it pre-mart; move it to `fe_fold_steps` |
| "No nodes were eligible for evaluation" | worker node pool with no clients — `MTB_WORKER_POOL` |
| `[NOMAD] Cannot determine exit status` | usually a *placement* failure wearing a task failure's clothes — the node is full |
| a task succeeded and produced nothing | a failed allocation whose exit code was never read; set `--nomad_delete_on_completion false` and the job survives for inspection |
| the database transfer is glacially slow | the head is using the public edge, not the internal endpoint |
| your edit to a script had no effect on the cluster | the image was not rebuilt; the analysis code lives inside it |
| the FE identity check fails | you changed FE code without bumping `FE_VERSION` — that is the lock working |

Worker jobs are deleted on completion by default, which is right for a healthy cluster and
wrong the moment something fails, because the failure is deleted with the job. The first thing
to do when debugging is turn that off.

---

## 11. The technique backlog

The architecture exists so that each remaining technique is a configuration line rather than a
fortnight. What has and has not been built, by layer:

| layer | technique | state |
|---|---|---|
| L0 | carrier-rate window | shipped in `build_mart` |
| L0 | `region_filter` — PE/PPE, ESX, IS/phage | implemented |
| L0 | `callability_filter` — EBR/RLC mask | implemented |
| L0 | `class_filter` — mutation classes | implemented |
| L0 | tier-1/2 gene filter | not built; needs the tier list |
| L1 | functional category · essentiality · PFAM | annotations built and shipped in `assets/` |
| L1 | LoF refinement — true LoF vs tolerated frameshift | not built |
| L1 | ESM-2 zero-shot | evaluation script exists, not a registered step |
| L1 | ΔΔG / biophysical | not built; highest effort |
| L2 | gene burden · AA property · variant cardinality | scripts exist, not registered steps |
| L4 | `top_n` prevalence | shipped |
| L4 | `rank_association` — chi², mutual information | implemented |
| L4 | `rank_mic` — Mann-Whitney on graded MIC | implemented |
| L4 | mRMR | not built |
| L5 | concordance — EBM · CATE · FWL | shipped, fold-internal |
| L5 | `select_cross_model` — RF ∩ MLP | implemented |
| L5 | Boruta · knockoffs · stability selection | not built |
| L5 | DoubleML · TMLE · causalml · causal-learn | not built |
| L5 | EBM `init_score` = catalogue | not built |

Whole-genome association methods — pyseer, treeWAS, bugwas, tree-based homoplasy — do not fit
this diagram at all. They need alignments and a phylogeny, so they are an upstream pipeline
rather than a layer.

**A warning worth having up front.** Five independent methods — invariance, phylogenetic
homoplasy, ESM-2, functional-category burden and curated resistance-gene burden — already
validate the concordance selection while adding no predictive signal. The prior on any further
technique improving AUC is low, and the project's own conclusion is that the FE negative space
is largely settled. The value here is not that one of these will find lift; it is that a
well-designed negative result is publishable in work whose argument is about honest evaluation
rather than about a leaderboard. Design an experiment that can produce a reportable negative,
and `experiments/README.md`'s rule 4 — state the prediction before the run — is what makes it
one.

---

## 12. Where the rest of it is

In this repository:

- `README.md` — the three entry branches, storage, profiles, containers, catalogue builds, the
  deposit format
- `experiments/README.md` — the rules that make experiments comparable
- `docs/fe-versioning.md` — `fe_version`, `fe_id`, `mart_id` and what each excludes
- `analysis/scripts/feature_mart/FE_CHANGELOG.md` — what changed in FE, and how it was verified
- `containers/README.md` — why there are four images and what each carries
- `site.env.example` — every deployment-specific setting, with the failure mode of each

In the research repository (`PHD-pub-mtb-resistotyper-ml/design/`), which is where the
methodology arguments live rather than the operating instructions:

- `02-feature-engineering.md`, `09-feature-mart.md` — the feature specification
- `03-modeling-pipeline.md` — the H2O / EBM / causal architecture and the stacking decision
  matrix
- `14-team-roadmap-abhi-ana.md` — the lane split, in detail
- `15-causal-concordance-and-multidrug.md` — the concordance method
- `19`–`25` — the per-technique literature and feasibility notes behind §11's table
- `23-fe-layers-and-the-label-boundary.md` — the layer diagram and the boundary argument

The platform context — why this runs the way it does on a federated African cluster — is
Chapter 6 of the dissertation, which points back here for the FE detail rather than restating
it.
