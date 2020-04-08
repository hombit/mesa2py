FROM python:3.7-stretch

ARG MESASDKVERSION=20.3.2
ARG MESAVERSION=12778

RUN ln -sfv /bin/bash /bin/sh

RUN apt-get update &&\
    apt-get install -y binutils bzip2 libc-dev libx11-dev libz-dev make subversion wget &&\
    apt-get clean -y &&\
    rm -rf /var/lib/apt/lists/* &&\
    truncate -s 0 /var/log/*log

RUN wget --no-verbose -U "" -O /mesasdk.tar.gz http://www.astro.wisc.edu/~townsend/resource/download/mesasdk/mesasdk-x86_64-linux-${MESASDKVERSION}.tar.gz &&\
    tar -xvf /mesasdk.tar.gz
#     rm /mesasdk.tar.gz

RUN svn co -r ${MESAVERSION} https://subversion.assembla.com/svn/mesa^mesa/trunk /mesa

COPY etc_ld.so.conf.d_mesasdk.conf /etc/ld.so.conf.d/mesasdk.conf
RUN ldconfig
