FROM mambaorg/micromamba:1.5.6
USER root

RUN apt-get update && apt-get install -y wget && apt-get clean

COPY env.yml /tmp/env.yaml
RUN micromamba install -y -n base -f /tmp/env.yaml && \
    micromamba clean --all --yes

RUN wget -O quarto-1.4.551-linux-amd64.tar.gz https://github.com/quarto-dev/quarto-cli/releases/download/v1.4.551/quarto-1.4.551-linux-amd64.tar.gz &&\
    tar -xvf /opt/quarto/quarto-1.4.551-linux-amd64.tar.gz -C /opt/quarto/ &&\
    rm /opt/quarto/quarto-1.4.551-linux-amd64.tar.gz

ENV PATH $PATH:/opt/conda/bin/:/opt/quarto/quarto-1.4.551-linux-amd64/bin/

WORKDIR /app
