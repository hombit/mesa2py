MESASDK_ROOT = /Applications/mesasdk
MESA_DIR = ./mesa

FC = ${MESASDK_ROOT}/bin/gfortran

FCFLAGS = -I ${MESA_DIR}/include -O3 -fopenmp
LDFLAGS = -L ${MESASDK_ROOT}/lib -lgomp -llapack -lblas -L ${MESA_DIR}/lib -lnet -leos -lkap -lrates -lchem -linterp_2d -linterp_1d -lnum -lf2crlibm -lmtx -lconst -lutils -lcrlibm

all: main

run: all
	MESA_DIR=$(realpath ${MESA_DIR}) ./main

main: opacity.o main.f90
	${FC} ${FCFLAGS} opacity.o main.f90 ${LDFLAGS} -o main

%.o: %.f90
	${FC}  ${FCFLAGS} -c $<
