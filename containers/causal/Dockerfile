# mtb-causal — the label-using layers: concordance selection, evaluation, ranking.
#
# The expensive image, and the one the 44 fold-internal selection jobs pull.
# econml and dowhy are why this cannot share an image with H2O: the two stacks
# pin incompatible numpy and scikit-learn, which is the same split that forced
# separate interpreters on the laptop.
#
# Versions match ~/.cache/mtb-causal-venv, the environment every published
# concordance result was produced in. numpy is held below 2 because econml 0.15
# and lightgbm 4.6 were built against the 1.x ABI here.
FROM python:3.11.8-slim

LABEL org.opencontainers.image.title="mtb-causal" \
      org.opencontainers.image.description="Causal-concordance selection, CV evaluation, label-using FE steps" \
      org.opencontainers.image.source="https://github.com/abhi18av-phd-projects/mtb-resistotyper-ml-pipeline" \
      org.opencontainers.image.licenses="EPL-2.0"
# The revision is passed at build time so a published image can be traced to
# the commit that produced it. An image tag pins an analysis version; the
# revision says which.
ARG GIT_REVISION=unknown
LABEL org.opencontainers.image.revision="${GIT_REVISION}"

RUN apt-get update && apt-get install -y --no-install-recommends \
        build-essential libgomp1 procps ca-certificates \
    && rm -rf /var/lib/apt/lists/*

# A CONSTRAINTS file, not an install-order trick. The previous version installed
# numpy 1.26.4 first and commented that "nothing below silently pulls 2.x over
# it". Something did: the built image carried numpy 2.4.6, because pip will
# happily upgrade an already-installed pin to satisfy a later requirement.
#
# The consequence was not an import error at build time, because the build's own
# check imports numpy, scipy, sklearn, lightgbm, econml, dowhy and interpret --
# none of which touch the broken path. It surfaced six hours into a selection
# job as `np.obj2sctype was removed in the NumPy 2.0 release`, raised from inside
# shap, which arrives transitively and was never pinned at all.
#
# PIP_CONSTRAINT applies to every resolution in this layer, including transitive
# ones, so the pin cannot be walked over.
RUN printf 'numpy==1.26.4\nshap==0.43.0\n' > /tmp/constraints.txt
ENV PIP_CONSTRAINT=/tmp/constraints.txt
RUN pip install --no-cache-dir numpy==1.26.4 \
 && pip install --no-cache-dir \
        shap==0.43.0 \
        scipy==1.15.3 \
        pandas==2.2.3 \
        pyarrow==24.0.0 \
        duckdb==1.5.4 \
        scikit-learn==1.4.2 \
        lightgbm==4.6.0 \
        interpret==0.7.3 \
        econml==0.15.0 \
        dowhy==0.13

# setuptools < 81 on purpose. Something in the causal stack still imports
# pkg_resources, which setuptools REMOVED in v81; a fresh environment therefore
# fails at import with ModuleNotFoundError while an older machine keeps working
# because its setuptools predates the removal. Pinning here makes the image
# reproduce the environment the results were produced in rather than whatever
# setuptools ships the day it is built.
RUN pip install --no-cache-dir "setuptools<81"

# Fail the BUILD rather than a six-hour selection job: importing every voter
# here is the cheapest place to catch an ABI break.
# Import shap too, and assert the numpy major version. The previous check passed
# on an image that could not run a selection: it imported everything except the
# library that actually broke, and asserted no version at all.
# EVALUATE_CV runs in this image and logs the honest-evaluation metrics. Without
# mlflow here the tracking degrades to a no-op and the lineage-held-out results
# -- the numbers the campaign exists to produce -- never reach the dashboard.
RUN pip install --no-cache-dir "mlflow-skinny==3.14.0"

RUN python -c "import numpy, scipy, sklearn, lightgbm, econml, dowhy, interpret, shap, mlflow; \
assert numpy.__version__.startswith('1.'), f'numpy pin lost: {numpy.__version__}'; \
print('voters import OK:', numpy.__version__, sklearn.__version__, 'shap', shap.__version__, 'mlflow', mlflow.__version__)"

WORKDIR /opt/mtb
COPY analysis/ /opt/mtb/analysis/
ENV PYTHONPATH=/opt/mtb
ENV MTB_PYTHON=/usr/local/bin/python
