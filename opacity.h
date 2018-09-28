#ifndef OPACITY_H
#define OPACITY_H


#include "macros.h"

typedef struct {
	int EOS_HANDLER, KAP_HANDLER;
	double XA[_SPECIES];
	double Y, ABAR, ZBAR, Z2BAR, YE;
	int* NET_ISO, CHEM_ID;
	double X, Z, Zfrac_C, Zfrac_N, Zfrac_O, Zfrac_Ne;
} Opacity;

Opacity init_Opacity();

void shutdown_Opacity(Opacity*);

void eos_PT(Opacity* op, double Pgas, double T,
	double* Rho, double* log10Rho,
	double* dlnRho_dlnPgas_const_T, double* dlnRho_dlnT_const_Pgas,
	double res[_SPECIES],
	double d_dlnRho_const_T[_SPECIES], double d_dlnT_const_Rho[_SPECIES],
	double d_dabar_const_TRho[_SPECIES], double d_dzbar_const_TRho[_SPECIES],
	int* ierr);


#endif  // OPACITY_H
