# Submitting this pipeline to the abc cluster

Operating instructions for a collaborator and, specifically, for **an agent working on the
collaborator's behalf**. It starts from nothing installed and assumes you have been handed
your own cluster credentials but have never submitted a job here.

Read [`collaborator-guide.md`](collaborator-guide.md) first if you have not: it covers what the
pipeline computes and how to change it. This document covers installing and configuring the
CLI, proving the cluster path works before trusting it, getting a run submitted, iterating
without paying for a full campaign each time, and which parameters are worth varying.

**§§1–4 are done once.** §§5–11 are the working loop.

> **For an agent.** §11 is a hard-rules list. Read it before running anything. Two of its rules
> exist because they have already been broken in this project: a command that prints a live
> credential into a transcript burns that credential, and switching the shared CLI context
> yanks it out from under every other session on the machine.

---

## 1. Install and configure the CLI

### Install

The CLI repository is **public**, so this needs no credential:

```bash
# into the current directory
curl -fsSL -H "Accept: application/vnd.github.raw+json" \
  "https://api.github.com/repos/abc-cluster/abc-cluster-cli/contents/scripts/install-abc.sh?ref=main" | sh

# or straight to /usr/local/bin/abc
curl -fsSL -H "Accept: application/vnd.github.raw+json" \
  "https://api.github.com/repos/abc-cluster/abc-cluster-cli/contents/scripts/install-abc.sh?ref=main" | sh -s -- --sudo
```

The installer is fetched through the GitHub **Contents API** rather than
`raw.githubusercontent.com` on purpose: `ref=main` on the CDN can serve an outdated snapshot.

```bash
abc version          # expect v0.1.82 or newer
abc self-update      # in place, thereafter
```

**Require v0.1.82 or newer.** Releases before it have two defects that bite this pipeline
directly: `abc data fetch` uploaded with the wrong credentials, and the pipeline head ran with
a memory default low enough to be OOM-killed into a restart loop that reported itself as
something else. Both were fixed on 2026-10-31.

**Never hand-roll a binary.** Install from a release or `abc self-update`; do not cross-compile
one into `/tmp` and copy it around. A forked toolchain makes it impossible to say afterwards
which CLI submitted a run.

For the local rungs of §7 you also need **Nextflow** (the pipeline manifest requires ≥ 24.10.0;
the pinned nf-nomad plugin releases require ≥ 26.04.3), a **JDK**, and **Docker**. A cluster run
needs none of them locally — the head job carries its own Nextflow and the work runs in
containers.

### Initialise local config

```bash
abc config init      # creates ~/.abc/config.yaml with a placeholder "default" context
abc config list      # every key, with the active context marked
abc config get <key>
```

Three things worth knowing about the local state:

- **`~/.abc/config.yaml`** holds every context, including access tokens. It is a secret file.
  Do not cat it, copy it into a chat, or commit it.
- **`~/.abc/db/local.db`** is a local SQLite holding run and data records; `abc report` reads it.
  Inspect with `abc localdb`.
- **Every `abc` invocation is audited locally** — timestamp, context, verb, redacted argv, exit
  code, duration. `ABC_NO_AUDIT=1` disables it. This is local bookkeeping, not surveillance, and
  it is what lets you reconstruct what you ran last week.

## 2. What you need handed to you

Eight values. None of them can be derived, and the last two are secrets.

| | variable | example | secret |
|---|---|---|---|
| 1 | context name | `seedling-ana` | no |
| 2 | API endpoint | `https://seedling.abc-cluster.cloud` | no |
| 3 | group bucket | `su-mbhg-bioinformatics` | no |
| 4 | Nomad namespace | `su-mbhg-bioinformatics` | no |
| 5 | head node pool | `platform` | no |
| 6 | worker node pool | `platform` — **not `compute`**; see §9 | no |
| 7 | internal S3 endpoint | `http://<internal-ip>:9000` — the API port `:9000`, never the console `:9001` | no |
| 8 | API access token · MinIO access key + secret | — | **yes** |

Values 3–7 are the same for everyone working on this pipeline and are already in
`site.env.example` as comments. Values 1, 2 and 8 are yours.

**What your MinIO key and secret are for, and what they are not for.** They are for reaching
object storage *directly* — `mc`, the AWS CLI, DuckDB's httpfs extension, a notebook reading a
parquet file. **A cluster run does not need them in your environment.** `abc pipeline run`
injects `AWS_ACCESS_KEY_ID`, `AWS_SECRET_ACCESS_KEY` and a scoped `NOMAD_TOKEN` into the head
job itself, from your identity on the server side. If you find yourself exporting AWS keys to
make a submission work, something else is wrong — check §9 before adding credentials to the
problem.

---

## 3. One-time setup

```bash
abc auth context add seedling-ana \
    --endpoint https://seedling.abc-cluster.cloud \
    --access-token "$(cat ~/.abc-token)"      # read from a file, never typed inline
```

Two notes.

**`abc auth context add` makes the new context active.** On your own machine that is harmless —
it is the only way to install the context. It is called out because the standing rule in this
project is *never* to run `abc auth context use`: the active context is global mutable state in
`~/.abc/config.yaml`, shared by every concurrent session on a machine, and switching it breaks
the others. Adding one context once is fine. Switching between contexts is not.

**Pass the context per command anyway:**

```bash
ABC_CLI_CONTEXT=seedling-ana abc <subcommand> ...
```

That is the convention the pipeline's launcher requires (it reads `ABC_CLI_CONTEXT` from
`site.env` and refuses to run without it), and it is the habit that keeps working the day you
have a second context.

Verify, without printing anything sensitive:

```bash
ABC_CLI_CONTEXT=seedling-ana abc auth whoami
ABC_CLI_CONTEXT=seedling-ana abc data list                 # should show your group's bucket
```

If `abc data list` shows no bucket, or a run later says *"your user has no group bucket (group:
none)"*, you are not a member of the group. That is an account state only an administrator can
change — ask, and until then pass `--work-dir` explicitly as §6 does.

### `site.env`

```bash
cd mtb-resistotyper-ml-pipeline
cp site.env.example site.env
$EDITOR site.env          # fill in values 1 and 3–7 from §2
```

`site.env` is gitignored. Everything in it describes *where* the pipeline runs; nothing in it
describes what it computes — that is what `experiments/*.yml` is for, which is why those files
carry no hostnames or buckets and can be re-run anywhere.

**Never commit `site.env`, and never paste its contents into a chat, an issue or a log.**

---

## 4. Prove it works before touching the pipeline

Four checks, escalating. Run them in order, the first time and after any change to credentials,
context or CLI version. Each one isolates a different layer, so the first failure tells you
where to look. **Do not skip to §6 with any of them failing** — every symptom in §9 is cheaper
to diagnose here than inside a twelve-hour run.

### 4.1 — The CLI can reach the cluster

```bash
abc doctor --context seedling-ana
```

`--context` runs the checks against a named context **without switching the active one**, which
is the whole reason to prefer it. Three groups run in order: config loads and has an address and
token; Nomad is reachable and healthy, and MinIO is probed when an S3 endpoint is configured; and
a minimal `raw_exec` batch job is submitted as a probe, waited on for up to 60 s, then purged.
Exit 0 means all three passed.

```bash
abc doctor --context seedling-ana --skip-job     # groups 1 and 2 only, if the cluster is busy
abc doctor --context seedling-ana --json         # machine-readable
abc doctor --context seedling-ana --bundle       # a REDACTED support bundle, safe to share
```

`--bundle` is the one diagnostic output here that is deliberately redacted and therefore safe to
send to someone. Prefer it over pasting raw command output when asking for help.

### 4.2 — Storage works, both directions

There is no storage flag on the verification job; prove it with a round-trip instead.

```bash
C=seedling-ana
ABC_CLI_CONTEXT=$C abc data list                              # your buckets
echo "hello $(date -u +%FT%TZ)" > /tmp/abc-probe.txt
ABC_CLI_CONTEXT=$C abc data push /tmp/abc-probe.txt \
    s3://<group-bucket>/users/<you>/abc-probe.txt
ABC_CLI_CONTEXT=$C abc data stat s3://<group-bucket>/users/<you>/abc-probe.txt
ABC_CLI_CONTEXT=$C abc data pull s3://<group-bucket>/users/<you>/abc-probe.txt /tmp/back.txt
diff /tmp/abc-probe.txt /tmp/back.txt && echo "round-trip ok"
ABC_CLI_CONTEXT=$C abc data delete s3://<group-bucket>/users/<you>/abc-probe.txt
```

`delete` is recoverable — it moves the object to the group's `trash/`. `remove` and `purge` are
not; do not reach for either during a check.

If `abc data list` shows nothing, see the group-membership note in §3.

### 4.3 — The cluster will run a job for you: `hello-cluster`

```bash
ABC_CLI_CONTEXT=seedling-ana abc job run hello-cluster \
    --cores 1 --mem 512M --node-pool platform --name probe-hello
```

`hello-cluster` is a **built-in verification workload** — a randomised `stress-ng` pass across
CPU, memory and I/O. It is the canonical proof that the client→scheduler→node path works for
your identity.

**Pass `--cores 1 --mem 512M`.** It defaults to 4 cores, and on this deployment the only
registered client in `platform` also runs eighteen of the platform's own services; a 4-core
reservation can sit unplaced behind them. Pass `--node-pool platform` explicitly for the reason
§6 gives — the context default may name a pool with no clients at all.

```bash
ABC_CLI_CONTEXT=seedling-ana abc job status probe-hello
ABC_CLI_CONTEXT=seedling-ana abc job logs   probe-hello
ABC_CLI_CONTEXT=seedling-ana abc job trace  probe-hello     # if it failed
ABC_CLI_CONTEXT=seedling-ana abc job stop   probe-hello
```

Two fallbacks, in order, if it never starts:

- `--sleep=120s` holds the allocation open before the work begins, so you can exec in and look
  around rather than inferring from a job that has already gone.
- `--driver docker` if the default `exec` driver hangs. Exec-driver hangs have been seen on this
  fleet's nodes (a kernel/Landlock interaction); `raw_exec` and `docker` work. **Report it rather
  than silently switching** — which driver works is cluster state worth knowing, not a
  preference.

### 4.4 — Nextflow works end to end: `nextflow-io/hello`

The first check that exercises the actual pipeline path — plugin resolution, node pools, the S3
work directory and the head job — on a workload that takes seconds instead of hours.

```bash
ABC_CLI_CONTEXT=seedling-ana abc pipeline run nextflow-io/hello \
    --head-pool platform --worker-pool platform \
    --work-dir s3://<group-bucket>/users/<you>/hello-workdir/ \
    --name probe-nf-hello --wait --logs
```

`--work-dir` is explicit because auto-derivation needs group membership; if yours is in place you
can drop it. Everything else is the minimum this deployment needs, and it is the same shape the
pipeline's own launcher produces.

When this prints four greetings, the path the real pipeline uses is proven. **That is the point
at which §6 is safe.** If it fails, the failure is in the deployment and not in the pipeline —
match it against §9 before changing anything in this repository.

## 5. The workbench

A code-server IDE on the cluster with your home directory bind-mounted, which survives
disconnection.

```bash
ABC_CLI_CONTEXT=seedling-ana abc workbench start --cores 4 --mem 8192
ABC_CLI_CONTEXT=seedling-ana abc workbench url      # prints the URL and password
ABC_CLI_CONTEXT=seedling-ana abc workbench status
ABC_CLI_CONTEXT=seedling-ana abc workbench stop     # releases the reservation
```

Use it for reading results, iterating on an analysis script, and running notebooks against data
already in the bucket. **Do not launch the pipeline from inside it** unless you have a reason —
a submission from your laptop and a submission from the workbench produce the same job, and the
laptop one is easier to interrupt.

`abc workbench url` prints a password. It is a credential: do not echo it anywhere that is
recorded. `abc workbench connect` mints a token for VS Code or Positron and has the same
property.

**Stop the session when you are done.** It holds cores and memory on a node that is also
running the platform's own services, and §9's queue problem gets worse while it is up.

---

## 6. Submitting a run

**Use the launcher. Not a bare `abc pipeline run`.**

```bash
scripts/launch.sh experiments/<name>.yml
```

The launcher carries four corrections this deployment needs, and each one, omitted, fails in a
way that looks like a different problem. Its header documents all four; the summary:

| correction | what happens without it |
|---|---|
| `--plugin nf-nomad@<pinned>` | the CLI resolves an older executor regardless of what the profile pins, and the profile's pin is silently ignored |
| `--head-pool` / `--worker-pool` | the worker pool default names a pool with no clients, so every worker job is accepted and then never placed: *"No nodes were eligible for evaluation"* |
| `NF_MINIO_ENDPOINT` forwarded | the head is handed the **public edge** for object storage, so a database on the same physical host is fetched Stellenbosch → Johannesburg → Stellenbosch. Measured 0.4 MB/s against 85 MB/s on a 10 GB file |
| `--env MTB_*` forwarded | `nextflow.config` resolves these with `System.getenv` *inside the head job*, so sourcing `site.env` locally is not enough. Unforwarded, the run reads `s3://su-example-group/`, which does not exist |

Extra flags pass straight through:

```bash
scripts/launch.sh experiments/baseline.yml --wait --logs
scripts/launch.sh experiments/baseline.yml --dry-run     # print the job spec, submit nothing
```

**`--dry-run` prints the generated job specification, and that specification contains your
injected credentials in clear text.** It is genuinely useful for checking pools, namespace and
work-dir — but pipe it through a filter, never into a transcript or a shared log:

```bash
scripts/launch.sh experiments/baseline.yml --dry-run 2>&1 \
  | grep -E 'node_pool|namespace|work-dir|NF_MINIO_ENDPOINT'
```

---

## 7. Iterating without paying for a campaign

A full campaign is 11 drugs × 4 lineages of an EBM, a causal forest and a refutation pass, plus
training: **on the order of twelve hours.** You will want to change something and run again.
These are the four rungs, cheapest first. Climb from the bottom every time.

### Rung 1 — a stub run (seconds, no cluster)

```bash
nextflow run . -stub-run --db /dev/null --drugs RIF
nextflow run . -stub-run --build_catalogue true --db /dev/null
```

Every process defines a stub, so the whole wiring runs with no container, no database and no
scheduler. This catches the class of fault that is otherwise only reachable by running the real
work first — a tuple published to a path built from the wrong field, two campaign arms merged
into one channel, manifests paired by directory instead of by fold. All three of those shipped
in this project, and each cost a full campaign to find.

**Run this after every change, including one you are sure about.**

### Rung 2 — the test profile (minutes, local Docker)

```bash
nextflow run . -profile containers,test --db /path/to/cryptic-slim.duckdb
```

One drug, one lineage, a 25-candidate pool, a 60-second H2O budget. Exercises every process
with real data. Use `-profile containers`, not `standard` — `standard` points at virtual
environments that only exist in the research-monorepo checkout layout.

### Rung 3 — the test profile on the cluster (minutes)

```bash
scripts/launch.sh experiments/<name>.yml --param drugs=RIF --param lineages=lineage2 --wait --logs
```

This is the rung that catches *deployment* faults rather than pipeline faults: pools, namespace,
endpoint, credentials, image pulls. Nothing below it can.

### Rung 4 — the real run

```bash
scripts/launch.sh experiments/<name>.yml
```

### The checkpoint that actually saves you time

**`--resume` does not cache across launches here. Do not rely on it.**

Every `abc pipeline run` generates a fresh `NXF_UUID`, a per-run work directory and a fresh
`NXF_CLOUDCACHE_PATH` — all three keyed on the new run id — and sets
`NXF_IGNORE_RESUME_HISTORY=true`. The `--resume` flag only appends `-resume <that fresh UUID>`,
so there is no prior session for Nextflow to match against. This was measured: a `--resume`
rerun re-executed the whole pipeline with **zero** cached processes.

What `-resume` does buy is recovery *within* one run, if the head job is preempted or restarted
mid-flight. That is worth having and is not cross-launch caching.

The real checkpoint is the **deposit**:

```bash
# once: write a self-contained DuckDB holding the state after FE and before training
scripts/launch.sh experiments/baseline.yml --param create_duckdb=mtb-fe-baseline.duckdb

# thereafter: start at training, locally, in minutes
nextflow run . -profile containers --from_duckdb mtb-fe-baseline.duckdb
```

`--from_duckdb` skips database assembly, feature engineering and causal selection entirely,
unpacks their results into the layout the training stages expect, and starts at training. If
your iteration is about **the learner and not the features**, this is the loop to live in: no
cluster, no twelve hours.

Two properties of a deposit that will produce a wrong answer if ignored:

- It holds **both campaign arms, and they are not interchangeable.** Units with a held-out
  lineage carry a selection made inside that training fold and feed evaluation only; the `none`
  unit carries a selection refit on all data and feeds the shipped model only. Crossing them
  gives you deployment models fitted on fold-restricted features and a "held-out" score
  computed on rows whose features already saw the held-out lineage.
- A deposit rerun **does not redo feature engineering**, so this run's `fe_*` settings do not
  describe the mart. Only `--fe_name` is read, as a grouping label for the tier report. Changing
  `fe_top_n` against a deposit changes nothing and silently mislabels the result.

### Three stage checkpoints worth knowing

The pipeline has three top-level entry branches, and a run takes exactly one:

| branch | flag | runs | stops before |
|---|---|---|---|
| train | *(none)* | database → cohort → mart → FE → concordance → training → evaluation | nothing |
| continue from a deposit | `--from_duckdb <file>` | training → evaluation | database build, FE and selection — unpacked, not recomputed |
| build a catalogue | `--build_catalogue true` | database → cohort → catomatic → comparison | FE, training and evaluation entirely |

---

## 8. Parameters worth varying

Set these in a committed `experiments/*.yml`, not on the command line, as soon as you care about
the answer. `experiments/baseline.yml` is the file every other one copies; the comparability
rules are in `experiments/README.md` and the first is **vary one thing**.

### Scope — the cost dial

| param | default | note |
|---|---|---|
| `drugs` | 11 drugs | the single biggest cost multiplier |
| `lineages` | `lineage1..4` | the other one; 11 × 4 is the campaign |
| `candidate_pool_size` | null (300) | dominates selection wall time — every candidate costs a CATE fit plus a refutation pass. The `test` profile sets 25 |

### Feature engineering

| param | default | decides |
|---|---|---|
| `fe_mutation_selector` | `top_n` | `top_n` \| `chi2` \| `mi` \| `all` |
| `fe_top_n` | 1000 | surviving mutation columns |
| `fe_max_carrier_frac` | 0.95 | upper carrier-rate bound — above it a mutation is a lineage marker, not a resistance variant |
| `fe_min_carrier_count` | 10 | lower bound |
| `fe_gene_level_flags` | true | gene-level aggregate columns |
| `fe_fold_strategy` | `stratified` | `stratified` \| `lineage` \| `country` — **use `lineage` for the honest regime** |
| `fe_n_folds` | 5 | |
| `fe_pre_steps` | `''` | label-free step chain, e.g. `region_filter:families=pe_ppe,esx` |
| `fe_fold_steps` | `''` | label-using step chain, e.g. `rank_association:method=chi2;top_n=1000` |
| `candidate_level` | `mutation` | `mutation` \| `gene` |

A step may carry parameters as `name:k=v;k=v`. Which chain a step belongs in is not a
preference: a step that reads the phenotype label and runs pre-mart chooses the column universe
while looking at the held-out lineage, and the registry refuses it at startup rather than
producing a mart that looks fine and is not.

### Training

| param | default | decides |
|---|---|---|
| `h2o_entrypoint` | `train_h2o_stacked` | which training script runs |
| `h2o_feature_set` | `concordant` | `denovo` \| `concordant` \| `all` |
| `h2o_runtime_secs` | 900 | AutoML budget, per drug per fold |
| `h2o_sort_metric` | `AUC` | AUC weights both classes equally, which is wrong where resistance is rare — BDQ, LZD and CFZ sit near 1–5%, and a classifier that never predicts R still scores well. **AUCPR is the honest choice there** |
| `h2o_balance_classes` | false | |
| `mlflow_uri` | unset | tracking is opt-in and non-fatal; a dashboard being down must not cost cluster time |

### Artefacts and debugging

| param | default | decides |
|---|---|---|
| `create_duckdb` | null | write the deposit — §7's checkpoint |
| `from_duckdb` | null | start from one |
| `nomad_delete_on_completion` | true | **set false while debugging** — a completed worker job is deleted with its allocation, events and logs, so a failure deletes its own evidence |
| `skip_db_check` | false | |

Catalogue-build parameters (`catalogue_*`) are documented in the README. Four of them will
silently produce a *plausible but wrong* catalogue if changed carelessly; read that section
before touching them.

---

## 9. What to expect, and what goes wrong

**Submission is bounded to three concurrent tasks on purpose.** Left unbounded, Nextflow
submitted five 6-CPU selections at once onto a node shared with eighteen of the platform's own
services, and Nomad answered *"cpu exhausted on 1 node(s)"*. That arrives as `[NOMAD] Cannot
determine exit status` — a placement failure wearing a task failure's clothes — so Nextflow
retried, and each retry added another blocked reservation to a queue that was already full.

**The work directory is group-shared**, which means the resume session key is too: a colleague
re-running the same pipeline reuses what your run did rather than repeating it. Putting the work
directory under a user or a project prefix quietly stops that sharing.

**The compute pool may be empty.** At the time of writing only `aither` (pool `platform`) and a
transient node are registered; `nomad01`–`04` are absent. That is why the worker pool is
`platform` and not `compute`, and it is a cluster state that changes — if placement starts
failing, check which pools have clients before changing anything in this repository.

| symptom | cause |
|---|---|
| "No nodes were eligible for evaluation" | worker pool with no clients — `MTB_WORKER_POOL` |
| "your user has no group bucket (group: none)" | your identity is not in the group; pass `--work-dir` explicitly, and ask to be added |
| `[NOMAD] Cannot determine exit status` | usually placement, not the task — the node is full |
| a task "succeeded" and produced nothing | a failed allocation whose exit code was never read; set `nomad_delete_on_completion=false` and inspect the surviving job |
| the run reads `s3://su-example-group/` | `MTB_*` not forwarded into the head job — use the launcher |
| database transfer at sub-MB/s | the head is on the public edge; `NF_MINIO_ENDPOINT` not forwarded |
| your script edit had no effect | the analysis code ships **inside the container image**. Changing a step means rebuilding and re-tagging the image before it runs on the cluster |
| the FE identity check fails | FE code changed without `FE_VERSION` being bumped — that is the lock working. Bump, add a `FE_CHANGELOG.md` entry, then `--lock` |

Watching and stopping a run:

```bash
ABC_CLI_CONTEXT=seedling-ana abc pipeline runs                 # recent runs
ABC_CLI_CONTEXT=seedling-ana abc pipeline logs <run-name>      # overview, per-task output
ABC_CLI_CONTEXT=seedling-ana abc pipeline stop <run-name>      # the run and its worker jobs
```

**Verify that a stop actually stopped.** A known defect in this deployment is that a stopped
run can keep working — the head goes away and workers carry on. After a stop, check
`abc pipeline runs` and the node's own job list before assuming the resources are released,
and say so if they are not.

---

## 10. Where the output goes

```
<outdir>/
  database/                        integrity checks
  cohort/                          the labelled substrate
  marts/<drug>/<fe-name>/          the mart and its sidecar — the reproducibility unit
  selection/<drug>/held-out-<L>/   concordance vote and per-fold selection manifest
  evaluation/<drug>/               cross-validation manifests
  models/<drug>/held-out-<L>/      H2O artefacts and h2o_manifest.json
  tiers/                           per-drug tier report — start here
  publish/                         the deposit, if create_duckdb was set
  reports/                         Nextflow's report, timeline, trace and DAG
```

On the cluster `outdir` defaults to `s3://<group-bucket>/<user-prefix>/results/<experiment>`,
which the workbench surfaces at `~/abc/me/<group>/results/`.

`held-out-none` is the **deployment refit, not a fold.** The four `held-out-lineageN`
directories are the honest evaluation. Do not average them with the first.

---

## 11. Hard rules for an agent

1. **Never run a command that prints a live credential into your transcript or any log.** That
   includes an unfiltered `--dry-run` (the job spec carries an injected MinIO secret and a Nomad
   token in clear text), `abc workbench url`, `abc workbench connect`, `cat site.env`, and
   `echo` of any token variable. Filter the output or redirect it to a file. A credential that
   reaches a transcript has to be rotated, and rotating it is someone else's afternoon.
2. **Never run `abc auth context use`**, or edit `active_context` in `~/.abc/config.yaml`. Pass
   `ABC_CLI_CONTEXT=<name>` on every command instead.
3. **Never commit `site.env`**, a token file, or anything under `work/`.
4. **Stage only your own files** — `git add <paths>`, never `git add -A`. Sibling sessions leave
   uncommitted work in these trees.
5. **Treat `~/.abc/config.yaml` as a secret file.** It holds every context's access token. Do
   not print it, copy it into a message, or put it in a repository.
6. **Install the CLI from a release or `abc self-update`.** Never cross-compile a binary into
   `/tmp` and copy it around — a forked toolchain makes it impossible to say which CLI submitted
   a given run.
7. **Run §4 before §6, and after any credential or version change.** Every symptom in §9 is
   cheaper to diagnose against `hello-cluster` than inside a twelve-hour campaign.
8. **Climb §7 from rung 1 every time.** Do not submit a twelve-hour campaign to find out whether
   a parameter name is spelled correctly. A stub run costs seconds and would have caught it.
9. **Do not relaunch to fix a typo.** Every launch is a fresh run with no cache; get the
   configuration right before resubmitting. Minimising rerun cost is a configuration discipline
   here, not a scheduler feature.
10. **Stop the workbench when you are not using it.** It holds cores on the same node the
   pipeline needs.
11. **Do not edit `main.nf`, `nextflow.config` or anything in `modules/`** to ask a scientific
   question. If a question requires it, that is a pipeline gap — say so. Everything you need is
   a params file, a step file, or a flag.
12. **Report a failure with its evidence.** Set `nomad_delete_on_completion=false`, keep the
   allocation, and quote what it actually said. A paraphrased Nomad error is almost always
   mis-diagnosed — four of the symptoms in §9 arrive wearing another failure's clothes.

### Entry point for a fresh session

Put this in a `CLAUDE.md` beside your checkout so it is read automatically:

```markdown
# Working notes

This repo is the mtb-resistotyper-ml feature-engineering and training pipeline.

- Operating instructions: `docs/cluster-submission.md`. Read §11 (hard rules) before
  running anything, §§1-4 to set up and verify the CLI, and §7 before submitting to the
  cluster.
- Before any cluster work in a fresh environment, run §4: `abc doctor --context <ctx>`,
  a storage round-trip, `abc job run hello-cluster`, then `abc pipeline run
  nextflow-io/hello`. Do not submit the real pipeline with any of those failing.
- What the pipeline computes, and how to change the features or the learner:
  `docs/collaborator-guide.md`.
- Cluster context is `ABC_CLI_CONTEXT=<my-context>`, passed per command. Never run
  `abc auth context use`.
- Submit with `scripts/launch.sh experiments/<name>.yml`, never a bare `abc pipeline run`.
- Never print a credential into the transcript: no unfiltered `--dry-run`, no
  `abc workbench url`, no `cat site.env`.
```
