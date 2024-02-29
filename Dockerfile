# FROM mesa_source
FROM ghcr.io/hombit/mesa-src:12778-python3.11-bookworm

RUN apt-get update &&\
    apt-get clean -y &&\
    rm -rf /var/lib/apt/lists/* &&\
    truncate -s 0 /var/log/*log

RUN pip3 install numpy cython

RUN sed -i'' -e 's/FCstatic =/FCstatic = -fPIC/' /mesa/utils/makefile_header
RUN echo "SPECIAL_C_FLAGS = -fPIC" >> /mesa/utils/makefile_header

env MESASDK_ROOT=/mesasdk
env MESA_DIR=/mesa

COPY call_mesa_script.sh setup.py /mesa2py/
WORKDIR /mesa2py
RUN source $MESASDK_ROOT/bin/mesasdk_init.sh
RUN python3 setup.py buildmesa

RUN pip3 install pytest pytest-subtests parameterized

COPY opacity.* /mesa2py/
RUN python3 setup.py install

COPY test.py /mesa2py/
RUN pytest test.py
