"""Compare a built catalogue against a reference, per drug.

This is the acceptance test of Chapter 6 §6.4 made executable. It reports what
Table 6.3 reports -- entries, agreement, resistance-call recovery -- so the
table regenerates rather than being transcribed once.

Two things it refuses to do. It does not claim identity: the catomatic paper
records that its own cohort is "seven samples fewer than used to build WHOv1;
the reason for this is unknown", so agreement is reported as a figure with a
delta and never as a pass/fail on exactness. And it does not silently reconcile
prediction vocabularies: the published catomatic catalogue is `RUS` while the
catalogue this tool serves is `RFUS`, and a comparison that quietly mapped one
onto the other would hide the difference it exists to measure.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

# GARC1 piezo CSVs carry these columns; only three are needed for a comparison.
KEY = "MUTATION"
CALL = "PREDICTION"


def _load(path: Path) -> dict[str, str]:
    """mutation -> prediction, for one drug's rows of a piezo catalogue."""
    out: dict[str, str] = {}
    with path.open(newline="") as fh:
        for row in csv.DictReader(fh):
            mut = (row.get(KEY) or "").strip()
            if not mut:
                continue
            out[mut] = (row.get(CALL) or "").strip()
    return out


def _drug_rows(path: Path, drug: str) -> dict[str, str]:
    """Same, restricted to one drug when the file carries many."""
    out: dict[str, str] = {}
    with path.open(newline="") as fh:
        reader = csv.DictReader(fh)
        has_drug = reader.fieldnames and "DRUG" in reader.fieldnames
        for row in reader:
            if has_drug and (row.get("DRUG") or "").strip().upper() != drug.upper():
                continue
            mut = (row.get(KEY) or "").strip()
            if mut:
                out[mut] = (row.get(CALL) or "").strip()
    return out


def compare(built: Path, reference: Path, drug: str, out_dir: Path) -> dict:
    a = _drug_rows(built, drug)
    b = _drug_rows(reference, drug)

    shared = sorted(set(a) & set(b))
    agree = [m for m in shared if a[m] == b[m]]
    disagree = [{"mutation": m, "built": a[m], "reference": b[m]}
                for m in shared if a[m] != b[m]]

    # Resistance-call recovery is the number the chapter quotes as "42 of 42":
    # of the reference's R calls, how many did the build also call R. It is
    # reported separately from overall agreement because a catalogue that
    # agrees on thousands of S calls and misses an R call is not "99% correct"
    # in any sense a laboratory cares about.
    ref_r = [m for m, c in b.items() if c == "R"]
    recovered_r = [m for m in ref_r if a.get(m) == "R"]

    report = {
        "drug": drug,
        "built_file": built.name,
        "reference_file": reference.name,
        "n_entries_built": len(a),
        "n_entries_reference": len(b),
        "n_shared": len(shared),
        "n_agree": len(agree),
        "pct_agree": round(100.0 * len(agree) / len(shared), 2) if shared else None,
        "n_disagree": len(disagree),
        # Two readings, both reported, because they answer different questions
        # and the published tables are not always explicit about which they
        # quote. The COUNT compares how many R calls each catalogue makes; the
        # RECOVERY asks how many of the reference's specific R mutations this
        # build also called R. A catalogue can match on count while disagreeing
        # on which mutations are resistant.
        "resistance_calls_built": len([m for m, c in a.items() if c == "R"]),
        "resistance_calls_reference": len(ref_r),
        "resistance_calls_recovered": len(recovered_r),
        "resistance_recovery": f"{len(recovered_r)} of {len(ref_r)}",
        "only_in_built": sorted(set(a) - set(b))[:50],
        "only_in_reference": sorted(set(b) - set(a))[:50],
        "n_only_in_built": len(set(a) - set(b)),
        "n_only_in_reference": len(set(b) - set(a)),
        "call_distribution_built": dict(Counter(a.values())),
        "call_distribution_reference": dict(Counter(b.values())),
        "disagreements": disagree[:200],
    }

    # Vocabulary mismatch is reported, not resolved. RUS vs RFUS is a real
    # difference between the published catalogue and the one this tool serves.
    vocab_a, vocab_b = set(a.values()), set(b.values())
    if vocab_a != vocab_b:
        report["prediction_vocabulary_differs"] = {
            "built": sorted(vocab_a), "reference": sorted(vocab_b),
            "note": ("The two catalogues do not use the same prediction values. "
                     "Agreement below is computed on exact string equality; "
                     "reconcile the vocabularies deliberately before reading it "
                     "as a reproduction figure."),
        }

    out_dir.mkdir(parents=True, exist_ok=True)
    (out_dir / f"comparison_{drug}.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--built", type=Path, required=True)
    p.add_argument("--reference", type=Path, required=True)
    p.add_argument("--drug", required=True)
    p.add_argument("--out", type=Path, default=Path("."))
    a = p.parse_args()
    r = compare(a.built, a.reference, a.drug, a.out)
    print(f"{r['drug']:<5} entries={r['n_entries_built']} "
          f"(ref {r['n_entries_reference']}) shared={r['n_shared']} "
          f"agree={r['n_agree']} ({r['pct_agree']}%) "
          f"R recovered={r['resistance_recovery']}")
    if "prediction_vocabulary_differs" in r:
        print("       WARNING: prediction vocabularies differ "
              f"{r['prediction_vocabulary_differs']['built']} vs "
              f"{r['prediction_vocabulary_differs']['reference']}")


if __name__ == "__main__":
    main()
