FROM mambaorg/micromamba:1.5.6
USER root
COPY ./deps/pybio/env.yaml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

COPY ./scripts /app/scripts
COPY ./tests /app/tests
COPY ./data /app/data
COPY .pytest.ini /app
ENV PATH $PATH:/opt/conda/bin/

WORKDIR /app
ENTRYPOINT ["pytest"]