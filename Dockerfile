ARG PYTHON_VERSION=3.10
FROM python:${PYTHON_VERSION}-bullseye

ENV MESASDK_ROOT=/mesasdk
ENV MESA_DIR=/mesa
ENV MESA2PY_SRC=/mesa2py

# It is easier to deal with bash
RUN ln -sf /bin/bash /bin/sh

# Install system requirements
RUN apt-get update &&\
    apt-get install -y binutils bzip2 libc-dev libx11-dev libz-dev make wget tcsh &&\
    apt-get clean -y &&\
    rm -rf /var/lib/apt/lists/* &&\
    truncate -s 0 /var/log/*log

# Get mesa-sdk saved into the image build by ./mesa-sdk/Dockerfile
COPY --from=ghcr.io/hombit/mesa-sdk:22.6.1 /mesasdk ${MESASDK_ROOT}
# Add MESA-SDK libraries to the system library path
RUN echo "${MESASDK_ROOT}/lib" > /etc/ld.so.conf.d/mesasdk.conf &&\
    echo "${MESASDK_ROOT}/math-slots/default/lib" >> /etc/ld.so.conf.d/mesasdk.conf &&\
    ldconfig

# MESA (reasonably) asks to not build itself under the root user
RUN useradd -ms /bin/bash mesa
# Create directories we need to be owned by the user while we are root
RUN mkdir -p ${MESA_DIR} ${MESA2PY_SRC} &&\
    chown -R mesa:mesa ${MESA_DIR} ${MESA2PY_SRC}
USER mesa

COPY --chown=mesa:mesa ./mesa ${MESA_DIR}

# Patch makefiles to compile position-independent static libraries
RUN sed -i'' -e 's/FCstatic =/FCstatic = -fPIC/' ${MESA_DIR}/utils/makefile_header &&\
    echo "SPECIAL_C_FLAGS = -fPIC" >> ${MESA_DIR}/utils/makefile_header

# We need kap module and its dependencies only
RUN export INSTALL_KAP_LINE= &&\
    head -n$(grep -n 'do_one_parallel kap' "${MESA_DIR}/install" | cut -f1 -d:) "${MESA_DIR}/install" > "/tmp/install_head" &&\
    tail -n+$(grep -n 'do_one' "${MESA_DIR}/install" | cut -f1 -d: | tail -n1) "${MESA_DIR}/install" | tail -n+2 > "/tmp/install_tail" &&\
    cat "/tmp/install_head" "/tmp/install_tail" > "${MESA_DIR}/install" &&\
    rm /tmp/install_head /tmp/install_tail &&\
    export -n INSTALL_KAP_LINE

# Build MESA
COPY build_mesa.sh /build_mesa.sh
RUN /build_mesa.sh

#VOLUME /mesa/data/eosDT_data/cache
#VOLUME /mesa/data/kap_data/cache

# Build mesa2py
RUN pip3 install --user numpy cython pytest pytest-subtests parameterized
COPY setup.py ${MESA2PY_SRC}/setup.py
COPY --chown=mesa:mesa src/ ${MESA2PY_SRC}/src/
WORKDIR ${MESA2PY_SRC}
RUN source $MESASDK_ROOT/bin/mesasdk_init.sh && \
    python3 setup.py install --user

COPY test.py ${MESA2PY_SRC}/test.py
RUN python -mpytest test.py -s
