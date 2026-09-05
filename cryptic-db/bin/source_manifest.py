#!/usr/bin/env python3
"""Checksum every input file before anything reads it.

A database is reproducible only if the inputs are identified by content. A
release directory that has been edited in place is otherwise indistinguishable
from one that has not, and every downstream provenance claim inherits that
uncertainty.
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


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", type=Path, required=True)
    ap.add_argument("--release", default="")
    ap.add_argument("--out", type=Path, required=True)
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

    manifest = {
        "release": release,
        "source_dir": str(a.src.resolve()),
        "n_files": len(files),
        "total_bytes": sum(p.stat().st_size for p in files),
        "files": {str(p.relative_to(a.src)): sha256(p) for p in files},
    }
    a.out.write_text(json.dumps(manifest, indent=2) + "\n")
    print(f"{len(files)} input file(s), release {release}")


if __name__ == "__main__":
    main()
