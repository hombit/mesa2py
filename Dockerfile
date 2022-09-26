# Set it from docker-compose / docker build
ARG PYTHON_VERSION=3.10
FROM ghcr.io/hombit/mesa-src:22.05-python${PYTHON_VERSION}-bullseye

RUN pip3 install numpy cython

RUN sed -i'' -e 's/FCstatic =/FCstatic = -fPIC/' /mesa/utils/makefile_header
RUN echo "SPECIAL_C_FLAGS = -fPIC" >> /mesa/utils/makefile_header

env MESASDK_ROOT=/mesasdk
env MESA_DIR=/mesa

RUN useradd -ms /bin/bash mesa
RUN chown -R mesa:mesa /mesa
USER mesa

COPY call_mesa_script.sh setup.py /mesa2py/
WORKDIR /mesa2py
RUN python3 setup.py buildmesa

RUN pip3 install pytest pytest-subtests parameterized

COPY opacity.* /mesa2py/
RUN python3 setup.py install

COPY test.py /mesa2py/
RUN pytest test.py
