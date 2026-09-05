# MLflow tracking server for the FE->metric dashboard, deployed as an
# abc-app on the cluster.
#
# Tracks H2O training runs from the Nextflow pipeline. sqlite backend + MinIO
# artifact store, both restored/persisted via the entrypoint, since the abc-app
# platform has no persistent volumes (verified during the crypticdb-llm deploy,
# 2026-07).
#
# The store moved from a project bucket (su-mtb-resistotyper-ml) to a PREFIX in
# the group bucket. ADR-0060: a project is a prefix, never a bucket of its own,
# because the shared pipeline cache and the workbench home are both keyed on the
# group bucket and a project bucket silently opts out of both.

FROM python:3.11-slim

WORKDIR /app

RUN apt-get update \
 && apt-get install -y --no-install-recommends ca-certificates \
 && rm -rf /var/lib/apt/lists/*

# mlflow (server + UI) + boto3 for the MinIO restore/snapshot + the S3
# artifact store. Pin mlflow to the version the local backfill used so the
# sqlite schema matches.
# anyio is pinned, and the pin is load-bearing. anyio 4.10 introduced a lazy
# import shim that does not expose `from_thread` as a module attribute, which is
# exactly what starlette's WSGI middleware reaches for when MLflow serves a
# request: every response 500s with
#   AttributeError: module 'anyio' has no attribute 'from_thread'
# The July build worked only because the then-current anyio resolved compatibly.
# An unpinned transitive dependency is a deploy that breaks on a day nothing
# changed.
RUN pip install --no-cache-dir "mlflow==3.14.0" "boto3>=1.34" "anyio==4.9.0"

COPY entrypoint.sh /entrypoint.sh
COPY snapshot.py /snapshot.py
RUN chmod +x /entrypoint.sh

EXPOSE 5000
ENV PYTHONUNBUFFERED=1
ENTRYPOINT ["/entrypoint.sh"]
