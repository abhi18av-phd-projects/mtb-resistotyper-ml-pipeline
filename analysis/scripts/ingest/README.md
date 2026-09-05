# ingest — the serving input contract (design/16 Stage 0)

The **ideal input** we ask users to send for drug-resistance prediction:
the CRyPTIC/gnomonicus **genotype + sample QC covariates**, in **GARC**
nomenclature. Deliberately **no catalogue interpretation** (no `effects`,
`predictions`, or R/S calls) — producing that is the model's job; shipping it
would leak the label.

## Files

| File | Role |
|---|---|
| `input_spec.schema.json` | The contract — JSON Schema (Draft 2020-12). Source of truth. |
| `reference_input.py` | Generates conformant instances from the CRyPTIC DB; self-validates. |
| `../../results/ingest/examples/*.input.json` | Grounded examples = Stage-1 round-trip fixtures. |

## The shape

```jsonc
{
  "schema_version": "1.0.0",
  "sample_id": "...",
  "nomenclature": "GARC",              // canon; other formats translate to GARC upstream
  "covariates": {                      // 100% populated in the CRyPTIC DB
    "lineage": "lineage2",             // -> cov__lineage_L1..L4
    "median_coverage": 116.5,          // -> cov__median_coverage
    "breadth": 99.08                   // -> cov__tb_breadth
  },
  "variants": [                        // ALL called variants (panel-wide), not just resistance genes
    {"gene": "rpoB", "mutation": "S450L", "frs": null, "coverage": 43.0, "is_minor": false}
  ]
}
```

- **Absence = reference**: only called variants are listed; anything not present
  is treated as wild-type (0) by the vectoriser.
- **`frs` is often null** even in the source DB — the de-novo model's `mean_frs`
  summary is computed over non-null values only.
- **Panel-wide variants are expected**: phylogenetic/synonymous calls (e.g.
  `Rv0008c S145P`) feed the de-novo model's summary features and the
  lineage-correlated columns. The concordant (deployable) model ignores all but
  a handful.

## Known edge case for the vectoriser (Stage 1)

**`X`-wildcard calls** (`rpoB S450X`, `katG S315X`) appear at ~1-2% of key
codons: the caller saw a change at that codon but couldn't resolve the residue
(often low coverage). GARC-valid, so the input carries them faithfully — but
`S450X` ≠ the model feature `raw__rpoB_S450L`. A **failed call at rpoB-450 is
not confirmed wild-type**, so the vectoriser must decide: treat `X` as 0
(pragmatic, current default), as missing, or flag it in the reasoning. This is a
Stage-1 correctness decision, recorded here so it isn't silently defaulted.

## Generate

```bash
.pixi/envs/default/bin/python -m analysis.scripts.ingest.reference_input \
    --db analysis/databases/duckdb/cryptic-slim.duckdb --auto \
    --out analysis/results/ingest/examples
# or a specific isolate: --uniqueid site.00.subj.1000347.lab.H111540004.iso.1
```
