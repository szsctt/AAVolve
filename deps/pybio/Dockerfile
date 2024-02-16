FROM mambaorg/micromamba:1.5.6
USER root
COPY env.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

ENV PATH $PATH:/opt/conda/bin/

WORKDIR /app
