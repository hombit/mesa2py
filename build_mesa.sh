#!/bin/bash

source ${MESASDK_ROOT}/bin/mesasdk_init.sh

function build_tool {
  CUR_DIR=$(pwd)

  NAME=$1
  TARGZFILE=$2
  DIRNAME="${TARGZFILE%.tar.gz}"

  cd ${MESA_DIR}/utils
  tar -xzvf ${TARGZFILE} -C /tmp
  cd /tmp/${DIRNAME}

  ./configure
  make $NAME
  mv -v ${NAME} ${MESA_DIR}/bin

  rm -rfv /tmp/${DIRNAME}
  cd ${CUR_DIR}
}

build_tool ndiff ndiff-2.00.tar.gz
build_tool makedepf90 makedepf90-2.8.8.tar.gz

cd ${MESA_DIR}
./install