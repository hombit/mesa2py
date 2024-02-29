FROM python:3.11-bookworm

ARG MESASDKVERSION=20.3.2
ARG MESAVERSION=12778

RUN ln -sfv /bin/bash /bin/sh

RUN apt-get update &&\
    apt-get install -y binutils bzip2 libc-dev libx11-dev libz-dev make wget &&\
    apt-get clean -y &&\
    rm -rf /var/lib/apt/lists/* &&\
    truncate -s 0 /var/log/*log

RUN wget --no-verbose -U "" -O /mesasdk.tar.gz http://user.astro.wisc.edu/~townsend/resource/download/mesasdk/mesasdk-x86_64-linux-${MESASDKVERSION}.tar.gz &&\
    tar -xvf /mesasdk.tar.gz

# BE CAREFUL, this is not a univeral link !!!
RUN wget --no-verbose -U "" -O /mesa.zip https://zenodo.org/records/3698354/files/mesa-r${MESAVERSION}.zip?download=1
RUN unzip /mesa.zip && rm mesa.zip
RUN mv mesa-r${MESAVERSION} mesa

COPY etc_ld.so.conf.d_mesasdk.conf /etc/ld.so.conf.d/mesasdk.conf
RUN ldconfig
