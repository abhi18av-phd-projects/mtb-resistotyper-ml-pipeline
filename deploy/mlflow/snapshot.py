#!/usr/bin/env python3
"""Persist mlflow.db back to object storage, periodically and on shutdown.

The inherited entrypoint restored the sqlite database on start and never wrote
it back, so the comment claiming it was "restored/persisted via the entrypoint"
was half true: every redeploy silently destroyed the entire training record
while leaving the artefacts behind in S3, orphaned from the rows that referenced
them. A tracking server that loses its history on restart is worse than none,
because the loss is invisible until someone looks for a run that is gone.

Snapshots use SQLite's own backup API rather than copying the file. A running
MLflow server writes continuously, and `cp` of a live sqlite database yields a
torn page as readily as a valid one; `Connection.backup` takes a consistent
snapshot of a database that is being written to, which is exactly the situation
here.
"""

from __future__ import annotations

import os
import signal
import sqlite3
import sys
import tempfile
import threading
import time
from pathlib import Path

DB = Path(os.environ.get("MLFLOW_DB_PATH", "/data/mlflow.db"))
BUCKET = os.environ["MLFLOW_DB_BUCKET"]
KEY = os.environ.get("MLFLOW_DB_OBJECT", "mlflow/mlflow.db")
INTERVAL = int(os.environ.get("MLFLOW_SNAPSHOT_SECONDS", "60"))
_stop = threading.Event()


def _log(m: str) -> None:
    print(f"[snapshot] {m}", file=sys.stderr, flush=True)


def snapshot(reason: str) -> None:
    if not DB.exists():
        return
    try:
        import boto3
        with tempfile.NamedTemporaryFile(suffix=".db", delete=False) as fh:
            tmp = Path(fh.name)
        src = sqlite3.connect(f"file:{DB}?mode=ro", uri=True)
        dst = sqlite3.connect(str(tmp))
        with dst:
            src.backup(dst)
        src.close(); dst.close()
        size = tmp.stat().st_size
        boto3.client("s3", endpoint_url=os.environ["ABC_MINIO_ENDPOINT"]).upload_file(
            str(tmp), BUCKET, KEY)
        tmp.unlink(missing_ok=True)
        _log(f"{reason}: wrote {size} bytes to s3://{BUCKET}/{KEY}")
    except Exception as exc:
        # Never fatal. A failed snapshot must not take down a tracking server
        # that is otherwise serving, and the next tick will try again.
        _log(f"{reason}: FAILED ({type(exc).__name__}: {exc})")


def _on_signal(signum, _frame) -> None:
    _log(f"signal {signum}: final snapshot before shutdown")
    snapshot("shutdown")
    _stop.set()


def main() -> None:
    signal.signal(signal.SIGTERM, _on_signal)
    signal.signal(signal.SIGINT, _on_signal)
    _log(f"every {INTERVAL}s -> s3://{BUCKET}/{KEY}")
    while not _stop.wait(INTERVAL):
        snapshot("periodic")


if __name__ == "__main__":
    main()
