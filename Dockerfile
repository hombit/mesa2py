# Set it from docker-compose / docker build
ARG PYTHON_VERSION=3.10
FROM ghcr.io/hombit/mesa-src:22.05-python${PYTHON_VERSION}-bullseye

RUN sed -i'' -e 's/FCstatic =/FCstatic = -fPIC/' /mesa/utils/makefile_header
RUN echo "SPECIAL_C_FLAGS = -fPIC" >> /mesa/utils/makefile_header

ENV MESASDK_ROOT=/mesasdk
ENV MESA_DIR=/mesa

RUN useradd -ms /bin/bash mesa
RUN chown -R mesa:mesa /mesa
RUN mkdir /mesa2py && chown -R mesa:mesa /mesa2py

USER mesa

RUN pip3 install --user numpy cython pytest pytest-subtests parameterized

COPY build_mesa.sh /mesa2py/
WORKDIR /mesa2py
RUN ./build_mesa.sh

VOLUME /mesa/data/eosDT_data/cache
VOLUME /mesa/data/kap_data/cache

COPY setup.py /mesa2py/
COPY --chown=mesa:mesa src/ /mesa2py/src/
RUN source $MESASDK_ROOT/bin/mesasdk_init.sh && \
    python3 setup.py install --user

COPY test.py /mesa2py/
RUN python -mpytest test.py -s
