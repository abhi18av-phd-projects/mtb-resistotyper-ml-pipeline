#!/usr/bin/env python3
"""Checksum every input file before anything reads it, and check none is missing.

A database is reproducible only if the inputs are identified by content. A
release directory that has been edited in place is otherwise indistinguishable
from one that has not, and every downstream provenance claim inherits that
uncertainty.

Recording what is present does not establish that everything is present, and the
difference is not academic. The v3.4.0 source directory was short one file,
VARIANTS.parquet, 1.29 GB and the largest in the release. Nothing complained: the
manifest faithfully checksummed the sixteen files that had arrived, the build
guarded the missing table as a release difference and skipped its stage, and the
database went to object storage without a 59-million-row table. The absence was
then read as evidence about the release rather than about the download, and it
survived into a thesis chapter.

So --expect compares the directory against the release's own file list, fetched
from Zenodo, and fails on anything missing or the wrong size. A checksum answers
"is this file what it claims to be". Only a manifest from the publisher answers
"are these all the files".
"""
from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

CHUNK = 1 << 20


def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(CHUNK), b""):
            h.update(block)
    return h.hexdigest()



def _expected_files(spec: str) -> dict[str, int]:
    """The release's own file list: {name: size}.

    `spec` is a Zenodo record id or a path to its API JSON. The publisher's
    manifest is the only thing that can say what a complete release looks like.
    """
    import json as _json
    if Path(spec).exists():
        rec = _json.loads(Path(spec).read_text())
    else:
        import urllib.request
        with urllib.request.urlopen(f"https://zenodo.org/api/records/{spec}") as r:
            rec = _json.load(r)
    return {f["key"]: f["size"] for f in rec["files"]}


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--release", default="")
    ap.add_argument("--out", type=Path, required=True)
    ap.add_argument("--expect", default="",
                    help="Zenodo record id, or a path to its API JSON, to check "
                         "completeness against. Strongly recommended.")
    ap.add_argument("--allow-missing", default="DATA_SCHEMA.pdf,RELEASE_NOTES.md",
                    help="comma-separated names that need not be present; "
                         "documentation, not data")
    a = ap.parse_args()

    files = sorted(
        p for p in a.src.rglob("*")
        if p.is_file() and p.suffix.lower() in {".parquet", ".csv", ".tsv", ".gz"}
    )
    if not files:
        raise SystemExit(f"no parquet/CSV inputs under {a.src}")

    release = a.release.strip()
    if not release:
        # cryptic-tables-v3.4.0 -> v3.4.0
        m = re.search(r"v\d+(?:\.\d+)*", a.src.resolve().name)
        release = m.group(0) if m else a.src.resolve().name

    if a.expect:
        expected = _expected_files(a.expect)
        optional = {x.strip() for x in a.allow_missing.split(",") if x.strip()}
        present = {p.name for p in files}
        # A gzipped lookup decompressed in place still satisfies its entry: the
        # data arrived, in a form the builder reads.
        present |= {n + ".gz" for n in present}
        problems = []
        for name, size in sorted(expected.items()):
            if name in optional:
                continue
            if name not in present:
                problems.append(f"MISSING  {name}  ({size:,} bytes)")
                continue
            local = a.src / name
            if local.exists() and local.stat().st_size != size:
                problems.append(
                    f"SIZE     {name}  local {local.stat().st_size:,} != release {size:,}")
        if problems:
            raise SystemExit(
                "source directory does not match the release:\n  "
                + "\n  ".join(problems)
                + "\n\nRefusing to build. A database assembled from an incomplete "
                  "release is not that release, and the gap is invisible afterwards: "
                  "a table absent because it was never downloaded looks exactly like "
                  "a table the release does not ship.")
        print(f"completeness: {len(expected)} files in the release, all accounted for")

    manifest = {
        "release": release,
        "checked_against": a.expect or None,
        "source_dir": str(a.src.resolve()),
        "n_files": len(files),
        "total_bytes": sum(p.stat().st_size for p in files),
        "files": {str(p.relative_to(a.src)): sha256(p) for p in files},
    }
    a.out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{len(files)} input file(s), release {release}")


if __name__ == "__main__":
    main()
