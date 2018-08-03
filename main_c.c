#include <stdio.h>

#include "opacity.h"


int main() {
    Opacity op = init_Opacity();
    double Pgas = 4e5;
    double T = 4e4;
    double Rho, log10Rho;
    double dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas;
    double res[_SPECIES], d_dlnRho_const_T[_SPECIES], d_dlnT_const_Rho[_SPECIES];
    double d_dabar_const_TRho[_SPECIES], d_dzbar_const_TRho[_SPECIES];
    int ierr = 0;
    eos_PT(&op, Pgas, T,
        &Rho, &log10Rho,
        &dlnRho_dlnPgas_const_T, &dlnRho_dlnT_const_Pgas,
	    res, d_dlnRho_const_T, d_dlnT_const_Rho,
	    d_dabar_const_TRho, d_dzbar_const_TRho,
	    &ierr);
    shutdown_Opacity(&op);
    printf("Rho\t%.17lg\n", Rho);
    return 0;
}

