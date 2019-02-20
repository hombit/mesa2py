#ifndef OPACITY_H
#define OPACITY_H


#include "macros.h"

int NUM_EOS_RESULTS;

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
	double*,
	int* ierr);

void kap_DT(Opacity* op, double Rho, double T, double lnfree_e,
	double* kappa, double* dlnkap_dlnRho, double* dlnkap_dlnT, int* ierr);

#endif  // OPACITY_H
