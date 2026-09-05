"""L0 — drop features in regions short reads cannot call reliably.

design/22 evaluated Marin et al. 2022 (Bioinformatics, doi 10.1093/bioinformatics/btac023),
which released per-position empirical base-level recall and a refined
low-confidence region set for H37Rv, measured from 36 isolates sequenced on both
Illumina and PacBio. On NC_000962.3, so the coordinates match CRyPTIC directly.

The reason this is worth a step rather than a footnote: we already know PE/PPE
features collapse under lineage-leave-one-out, but not WHY, and the two candidate
explanations imply different things. If they collapse because PE/PPE variants
track lineage, only PE/PPE is affected. If they collapse because those regions are
where short-read calls fail -- and mappability failure is itself lineage-dependent,
since divergence from the reference drives it -- then low-callability positions
OUTSIDE PE/PPE are contaminated too, and no gene-name filter reaches them.

This step tests that. Both explanations justify the PE/PPE filter, so it is safe
either way; only the second predicts anything about the rest of the genome.

Label-free: callability is a property of the reference and the sequencing
technology, measured on isolates unrelated to this cohort.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from analysis.scripts.fe_steps.registry import Step, register


def load_regions(path: Path) -> list[tuple[int, int]]:
    """A BED of low-confidence intervals -> sorted (start, end) in 1-based inclusive."""
    spans = []
    for line in Path(path).read_text().splitlines():
        if not line.strip() or line.startswith(("#", "track", "browser")):
            continue
        parts = line.split()
        if len(parts) < 3:
            continue
        spans.append((int(parts[1]) + 1, int(parts[2])))  # BED is 0-based half-open
    return sorted(spans)


def run(mart: pd.DataFrame, *, regions: str | None = None,
        annotations: pd.DataFrame | None = None, **_) -> tuple[pd.DataFrame, dict]:
    if not regions:
        raise ValueError(
            "callability_filter needs --regions pointing at a low-confidence BED "
            "(RLC_Regions.H37Rv.bed from farhat-lab/mtb-illumina-wgs-evaluation). "
            "design/22 evaluated the resource; it is a reference checkout and is "
            "deliberately not vendored into this repo.")
    spans = load_regions(Path(regions))
    if annotations is None:
        raise ValueError(
            "callability_filter needs the gene annotation map for coordinates; "
            "pass --annotations gene_annotations.parquet")

    coord_cols = [c for c in ("start", "end", "gene_start", "gene_end")
                  if c in annotations.columns]
    if len(coord_cols) < 2:
        raise ValueError(
            f"annotations carry no coordinate columns (have {list(annotations.columns)}); "
            "the gene map must supply start/end for a positional filter")
    start_col, end_col = coord_cols[0], coord_cols[1]
    key = annotations.columns[0]

    def overlaps(lo: int, hi: int) -> bool:
        return any(not (hi < s or lo > e) for s, e in spans)

    low_conf_genes = {
        str(r[key]) for _, r in annotations.iterrows()
        if pd.notna(r[start_col]) and pd.notna(r[end_col])
        and overlaps(int(r[start_col]), int(r[end_col]))
    }

    removed = [c for c in mart.columns
               if c.startswith("raw__")
               and c[len("raw__"):].partition("_")[0] in low_conf_genes]
    kept = mart.drop(columns=removed)

    n_mut = sum(1 for c in mart.columns if c.startswith("raw__"))
    note = {
        "step": "callability_filter",
        "regions_file": str(regions),
        "n_low_confidence_spans": len(spans),
        "bp_low_confidence": sum(e - s + 1 for s, e in spans),
        "genes_overlapping": len(low_conf_genes),
        "dropped": len(removed),
        "mutation_features_before": n_mut,
        "mutation_features_after": n_mut - len(removed),
        "resolution": "gene-level overlap; per-position masking needs variant coordinates",
    }
    return kept, note


register(Step(
    name="callability_filter",
    layer="universe",
    label_free=True,
    summary="drop features in empirically low-callability regions (Marin 2022 EBR/RLC)",
    run=run,
))
