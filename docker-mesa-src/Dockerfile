# Set it from docker-compose / docker build
ARG PYTHON_VERSION=3.10
FROM python:${PYTHON_VERSION}-bullseye

ARG MESASDKVERSION=22.6.1
ARG MESAVERSION=r22.05.1
ARG ZENODOURL=https://zenodo.org/record/6547951/files/mesa-${MESAVERSION}.zip?download=1

RUN ln -sfv /bin/bash /bin/sh

RUN apt-get update &&\
    apt-get install -y binutils bzip2 libc-dev libx11-dev libz-dev make wget &&\
    apt-get clean -y &&\
    rm -rf /var/lib/apt/lists/* &&\
    truncate -s 0 /var/log/*log

RUN wget --no-verbose -U "" -O /mesasdk.tar.gz http://www.astro.wisc.edu/~townsend/resource/download/mesasdk/mesasdk-x86_64-linux-${MESASDKVERSION}.tar.gz &&\
    tar -xvf /mesasdk.tar.gz
    # Keep it for the future when it will be removed from the web!
    # rm /mesasdk.tar.gz

RUN wget --no-verbose -O/mesa.zip ${ZENODOURL} && \
    unzip /mesa.zip && \
    rm /mesa.zip && \
    mv /mesa-${MESAVERSION} /mesa

COPY etc_ld.so.conf.d_mesasdk.conf /etc/ld.so.conf.d/mesasdk.conf
RUN ldconfig
