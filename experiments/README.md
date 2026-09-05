# Experiments

One file per question. Locally: `nextflow run .. -params-file experiments/<name>.yml`.
On the cluster: **`scripts/launch.sh experiments/<name>.yml`**, which carries four
corrections a bare `abc pipeline run` needs on this deployment — the nf-nomad
plugin pin, the worker node pool, and the internal MinIO endpoint. Each of them
fails in a way that looks like something else; the script's header says how.

A params file is the whole configuration of a run, so committing one makes an experiment a
reviewable artefact rather than a shell command somebody remembers. Six months later the file
says exactly what was varied and what was held fixed; a terminal history does not.

**The rules that make these comparable:**

1. **Vary one thing.** Everything else is copied verbatim from `baseline.yml`. A file that
   changes the FE *and* the drug list answers no question.
2. **Name the question, not the settings.** `does-callability-reach-beyond-pe-ppe`, not
   `callability-v2`.
3. **Give every experiment its own `outdir`.** Runs must not overwrite each other, and two
   result sets in one directory cannot be diffed.
4. **State the question and the prediction in the header comment.** A prediction written before
   the run is evidence; one written after is a description.
5. **Pin the database.** The `db` path fixes which database the answer is about — full or slim,
   and which release. Two runs against different databases are not a comparison.
