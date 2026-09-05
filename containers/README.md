# Pipeline container images

Three images, because the stacks are mutually exclusive.

| image | carries | used by |
|---|---|---|
| `mtb-fe` | duckdb · pandas · pyarrow | cohort, mart build, label-free FE steps, tier report |
| `mtb-causal` | econml · dowhy · interpret · lightgbm · sklearn 1.4 · **numpy 1.26** | concordance selection, CV evaluation, label-using FE steps |
| `mtb-h2o` | h2o 3.46 · JRE 17 · sklearn 1.7 · **numpy 2.2** | AutoML, base models, stacked ensembles, MOJO export |

`mtb-causal` and `mtb-h2o` cannot be merged: econml 0.15 and lightgbm 4.6 are built against
the numpy 1.x ABI, while this project's H2O environment is numpy 2.2. That is the same split
that forced three interpreters on the laptop, made explicit.

Versions are pinned to the environments the published results were produced in, so a cluster
run reproduces a laptop run rather than approximating it. Both non-trivial images assert their
imports at **build** time — `mtb-causal` imports every voter, `mtb-h2o` completes a JVM
handshake — because an ABI break is far cheaper to find in a build than six hours into a
selection job.

## The analysis code ships inside the image

Each image `COPY`s `analysis/` to `/opt/mtb` and sets `PYTHONPATH`. That is what lets the
pipeline repository stay pure Nextflow: a Nomad task needs no shared filesystem and no
checkout to find `analysis.scripts.feature_mart.build_mart`.

The consequence to remember: **changing an analysis script means rebuilding the image.** The
pipeline is versioned by its git tag; the analysis code is versioned by the image tag, and the
two are joined in `nextflow.config`.

## Building

```bash
./build.sh                 # all three
./build.sh causal          # just one
TAG=v0.2.0 ./build.sh      # override
```

Run it anywhere with a docker daemon and push access to the registry in `MTB_REGISTRY`
(see `site.env.example`). Images are published to GitHub Container Registry so the cluster
can pull them. An earlier arrangement used a host-local naming convention with no registry
behind it, which does not survive Nomad's
docker driver resolves the tag from the local image cache.
