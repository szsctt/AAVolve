# syntax=docker/dockerfile:1
# https://docs.docker.com/build/building/multi-stage/
# sudo docker build . -t szsctt/c3poa
# sudo docker run --rm -it szsctt/c3poa /bin/bash
# sudo docker push szsctt/c3poa

#FROM mambaorg/micromamba:latest
#
#RUN micromamba install --yes --name base --channel conda-forge --channel bioconda \
#	cython python=3.7 numpy setuptools wheel scipy tqdm racon blat mappy &&\
#	micromamba clean --all --yes
#
#ARG MAMBA_DOCKERFILE_ACTIVATE=1
#RUN pip install pyabpoa==1.4 editdistance

# racon
FROM ubuntu:20.04 as racon-builder

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update
RUN apt install build-essential wget unzip libpng-dev git cmake -y

RUN git clone --recursive https://github.com/isovic/racon.git racon
RUN cd racon &&\
mkdir build &&\
cd build &&\
cmake -DCMAKE_BUILD_TYPE=Release .. &&\
make    
    
# conk, python and C3POa
FROM python:3.7 as pybuilder

ARG DEBIAN_FRONTEND=noninteractive
RUN apt-get update && apt install \
        build-essential git cmake libz-dev wget tar libsimde-dev -y

# build conk wheel
RUN python3 -m pip install --user --upgrade wheel setuptools Cython
RUN git clone https://github.com/rvolden/conk
RUN cd conk &&\
      python3 setup.py sdist bdist_wheel

# build pyabpoa wheel
RUN pip install pyabpoa==1.4.0

# build other wheels
RUN  pip wheel -w wheels numpy scipy pyabpoa==1.4.0 mappy Cython tqdm editdistance

# clone C3POa
RUN wget https://github.com/christopher-vollmers/C3POa/archive/refs/tags/v3.1.tar.gz &&\
    tar -xvf v3.1.tar.gz 

# final image
FROM python:3.7-slim

RUN apt-get update && apt install -y \
        pigz wget libcurl4 \
    && rm -rf /var/lib/apt/lists/*

RUN mkdir -p deps/ C3POa/blat C3POa/racon/build/bin/
COPY --from=pybuilder C3POa-3.1 C3POa
RUN chmod +x C3POa/C3POa.py
COPY --from=racon-builder /racon/build/bin/racon C3POa/racon/build/bin/
COPY --from=pybuilder /conk/dist/conk-1.0.0-cp37-cp37m-linux_x86_64.whl deps/
COPY --from=pybuilder /wheels deps/wheels
#COPY --from=pybuilder /root/.local /root/.local
RUN wget -O C3POa/blat/blat http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/blat/blat &&\
    chmod +x C3POa/blat/blat
ENV PATH="${PATH}:/deps/:C3POa/blat:C3POa/racon/build/bin/:C3POa"

RUN pip install --no-cache-dir deps/conk-1.0.0-cp37-cp37m-linux_x86_64.whl
RUN pip install --no-cache-dir deps/wheels/*
