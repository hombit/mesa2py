MESASDK_ROOT = /Applications/mesasdk
MESA_DIR = ./mesa

FC = ${MESASDK_ROOT}/bin/gfortran
CC = ${MESASDK_ROOT}/bin/gcc
CPP = ${MESASDK_ROOT}/bin/cpp

CCFLAGS = -O3 -g
FCFLAGS = ${CCFLAGS} -I ${MESA_DIR}/include -fopenmp
LDFLAGS = -g -L ${MESASDK_ROOT}/lib -lgomp -llapack -lblas -L ${MESA_DIR}/lib -lnet -leos -lkap -lrates -lchem -linterp_2d -linterp_1d -lnum -lf2crlibm -lmtx -lconst -lutils -lcrlibm

all: main_fort main_c

run: all
	MESA_DIR=$(realpath ${MESA_DIR}) ./main_fort
	MESA_DIR=$(realpath ${MESA_DIR}) ./main_c

main_fort: opacity.o main_fort.f90
	${FC} ${FCFLAGS} opacity.o main_fort.f90 ${LDFLAGS} -o main_fort

main_c: opacity.o main_c.o
	${FC} ${FCFLAGS} opacity.o main_c.o ${LDFLAGS} -o main_c

%.o: %.f90
	${CPP} $< .pp.$<
	${FC} ${FCFLAGS} -c .pp.$< -o $@

%.o: %.c opacity.h macros.h
	${CC} ${CCFLAGS} -c $<
