# mtb-fe — the label-free layers: cohort, mart build, universe filters, reporting.
#
# Deliberately the smallest of the three. Nothing here reads the phenotype label
# for ranking, so nothing here needs a modelling stack; adding one would only
# make the image that runs on every drug slower to pull.
#
# Versions are pinned to the pixi environment the marts were actually built in,
# so a cluster run reproduces a laptop run rather than approximating it.
FROM python:3.11.13-slim

LABEL org.opencontainers.image.title="mtb-fe" \
      org.opencontainers.image.description="Feature-mart construction and label-free FE steps" \
      org.opencontainers.image.source="https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-pipeline" \
      org.opencontainers.image.licenses="EPL-2.0"
# The revision is passed at build time so a published image can be traced to
# the commit that produced it. An image tag pins an analysis version; the
# revision says which.
ARG GIT_REVISION=unknown
LABEL org.opencontainers.image.revision="${GIT_REVISION}"

RUN apt-get update && apt-get install -y --no-install-recommends \
        procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*

RUN pip install --no-cache-dir \
        duckdb==1.3.0 \
        pandas==2.3.0 \
        pyarrow==18.1.0 \
        numpy==2.2.6 \
        scipy==1.15.3

# The analysis code travels IN the image, so the pipeline repo stays pure
# Nextflow and a task needs no shared filesystem to find its scripts.
WORKDIR /opt/mtb
COPY analysis/ /opt/mtb/analysis/
ENV PYTHONPATH=/opt/mtb

# The modules invoke `python -m analysis.scripts...`; make that the default.
ENV MTB_PYTHON=/usr/local/bin/python
