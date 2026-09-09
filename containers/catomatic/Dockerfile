# mtb-catomatic — algorithmic catalogue construction.
#
# A fourth image rather than a fourth stack in an existing one. catomatic pulls
# piezo and a statistics stack; the FE image runs on every drug of every mart
# build and should not grow to carry a dependency only the catalogue path uses.
#
# catomatic is consumed AS PUBLISHED, from PyPI, with its own command-line
# interface. That is deliberate and it is the claim: the chapter argues that a
# published, notebook-bound analysis survives being decomposed into declared
# cluster jobs. Vendoring or rewriting it would forfeit exactly that.
FROM python:3.11.13-slim

LABEL org.opencontainers.image.title="mtb-catomatic" \
      org.opencontainers.image.description="Catalogue construction with catomatic" \
      org.opencontainers.image.source="https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-pipeline" \
      org.opencontainers.image.licenses="EPL-2.0"

ARG GIT_REVISION=unknown
LABEL org.opencontainers.image.revision="${GIT_REVISION}"

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# Pinned, like every other image here. catomatic is MIT and its dependencies
# already include piezo, so the catalogue it emits is readable by the same
# library that reads the WHO catalogue in Layer 1.
RUN pip install --no-cache-dir \
        catomatic==0.2.1 \
        duckdb==1.3.0 \
        pandas==2.3.0 \
        pyarrow==18.1.0

WORKDIR /opt/mtb
COPY analysis/ /opt/mtb/analysis/
ENV PYTHONPATH=/opt/mtb

ENV MTB_PYTHON=/usr/local/bin/python
