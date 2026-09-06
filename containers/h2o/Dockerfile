# mtb-h2o — AutoML, base models, stacked ensembles.
#
# H2O boots its own JVM inside the task and is torn down with it, so this image
# carries a JRE and nothing else needs one. A long-lived H2O service would be
# shared mutable state across tasks; per-task keeps the heap bounded by the
# allocation Nomad actually reserved.
#
# H2O_MEM is set by the process from task.memory, not here — an image-level
# default would be a promise the scheduler never made.
FROM eclipse-temurin:17-jre-jammy

LABEL org.opencontainers.image.title="mtb-h2o" \
      org.opencontainers.image.description="H2O AutoML training and MOJO export" \
      org.opencontainers.image.source="https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-pipeline" \
      org.opencontainers.image.licenses="EPL-2.0"
# The revision is passed at build time so a published image can be traced to
# the commit that produced it. An image tag pins an analysis version; the
# revision says which.
ARG GIT_REVISION=unknown
LABEL org.opencontainers.image.revision="${GIT_REVISION}"

RUN apt-get update && apt-get install -y --no-install-recommends \
        python3.11 python3-pip python3.11-venv procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# A virtual environment rather than the system interpreter. --break-system-packages
# does not exist in this base image's pip (Ubuntu 22.04 predates it), and a venv
# answers the question that flag exists to override rather than overriding it.
ENV VIRTUAL_ENV=/opt/venv
RUN python3.11 -m venv "$VIRTUAL_ENV"
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

RUN pip install --no-cache-dir --upgrade pip \
 && pip install --no-cache-dir \
        numpy==2.2.6 \
        pandas==2.3.0 \
        pyarrow==18.1.0 \
        duckdb==1.3.0 \
        scikit-learn==1.7.1 \
        h2o==3.46.0.7 \
        interpret==0.6.9

# setuptools < 81 on purpose. Something in the causal stack still imports
# pkg_resources, which setuptools REMOVED in v81; a fresh environment therefore
# fails at import with ModuleNotFoundError while an older machine keeps working
# because its setuptools predates the removal. Pinning here makes the image
# reproduce the environment the results were produced in rather than whatever
# setuptools ships the day it is built.
RUN pip install --no-cache-dir "setuptools<81"

# The tracking client, pinned to the server's version so the sqlite schema and
# the REST contract cannot drift apart between the two. anyio is pinned for the
# same reason it is pinned in the server image: 4.10's lazy-import shim breaks
# starlette's WSGI bridge, and an unpinned transitive dependency is a run that
# fails on a day nothing changed.
RUN pip install --no-cache-dir "mlflow==3.14.0" "anyio==4.9.0"

# H2O is useless without a working JVM handshake; prove it at build time.
RUN java -version && python -c "import h2o, numpy; from interpret.glassbox import ExplainableBoostingRegressor; print('h2o', h2o.__version__, 'numpy', numpy.__version__, 'interpret ok')"

WORKDIR /opt/mtb
COPY analysis/ /opt/mtb/analysis/
ENV PYTHONPATH=/opt/mtb
ENV JAVA_HOME=/opt/java/openjdk
ENV MTB_PYTHON=/opt/venv/bin/python
