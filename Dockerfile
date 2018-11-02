FROM python as builder

ARG MESASDKVERSION=20180822
ARG MESAVERSION=10398

RUN ln -sfv /bin/bash /bin/sh

RUN apt-get update &&\
    apt-get install -y binutils bzip2 libc-dev libx11-dev libz-dev make subversion wget pbzip2 pigz &&\
    apt-get clean -y &&\
    rm -rf /var/lib/apt/lists/* &&\
    truncate -s 0 /var/log/*log

RUN ln -sf $(which pbzip2) $(which bzip2)

RUN wget --no-verbose -U "" -O /mesasdk.tar.gz http://www.astro.wisc.edu/~townsend/resource/download/mesasdk/mesasdk-x86_64-linux-${MESASDKVERSION}.tar.gz &&\
    tar -I pigz -xvf /mesasdk.tar.gz &&\
    rm /mesasdk.tar.gz

RUN svn co -r ${MESAVERSION} svn://svn.code.sf.net/p/mesa/code/trunk /mesa

# Add -fPIC to utils/makefile_header
RUN sed -i'' -e 's/FCstatic =/FCstatic = -fPIC/' /mesa/utils/makefile_header
RUN echo "SPECIAL_C_FLAGS = -fPIC" >> /mesa/utils/makefile_header
# Add -fPIC to crlibm/make/makefile
RUN sed -i'' -e 's/COMPILE_C = $(CC)/COMPILE_C = $(CC) -fPIC/' /mesa/crlibm/make/makefile

#WORKDIR /mesa
#RUN export MESASDK_ROOT=/mesasdk &&\
#    export MESA_DIR=/mesa &&\
#    source /mesasdk/bin/mesasdk_init.sh &&\
#    ./install_mods_in_parallel

RUN pip install numpy cython

COPY call_mesa_script.sh setup.py /mesa2py/
WORKDIR /mesa2py
RUN export MESASDK_ROOT=/mesasdk &&\
    export MESA_DIR=/mesa &&\
    python setup.py buildmesa

COPY . /mesa2py
RUN export MESASDK_ROOT=/mesasdk &&\
    export MESA_DIR=/mesa &&\
    python setup.py install
