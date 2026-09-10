# Metro map

`metro_map.svg` is the pipeline drawn as a transit map: one line per route
through it, one station per process. `metro_map.mmd` is its source.

The map is derived from the workflow, not drawn from memory. Nextflow writes the
process DAG of any run to `<outdir>/reports/dag.mmd`, and a stub run produces it
in seconds with no container, database or cluster, so each entry branch can be
captured on its own:

    nextflow run main.nf -stub-run --db /dev/null --drugs RIF,INH --create_duckdb x.duckdb
    nextflow run main.nf -stub-run --from_duckdb <any file>
    nextflow run main.nf -stub-run --build_catalogue true --db /dev/null --catalogue_reference <csv>

Those three DAGs were combined by hand into `metro_map.mmd`. The DAG is identical
for every drug, so lines are **routes**, not drugs: the evaluation arm, the
deployment arm, re-entry from a deposit, and the catalogue build. Open circles
are steps that are off by default; the square is the integrity gate.

Render with [nf-metro](https://github.com/seqeralabs/nf-metro) (2.0.0; Python 3.11+):

    pip install nf-metro
    nf-metro render docs/images/metro_map.mmd -o docs/images/metro_map.svg --mode light

**Keep it honest when the workflow changes.** Re-run the stub DAG for the branch
you touched and check the map still matches it. Two points to watch: the map ends
in *H2O models* because that is what `TRAIN_H2O` produces -- the deployable
L2-logistic bundles are made by `bin/train-release-bundles.sh`, outside the
workflow -- and *Label-free FE* is drawn as optional because `FE_PRE`'s output is
discarded unless `--fe_pre_steps` is set.
