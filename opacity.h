#ifndef OPACITY_H
#define OPACITY_H

int NUM_EOS_RESULTS;
int NUM_CHEM_ISOS_POINTER;
int SOLSIZE;

typedef struct {
	int EOS_HANDLER, KAP_HANDLER;
	int SPECIES;
	double X, Y, Z, XC, XN, XO, XNe, ABAR, ZBAR, Z2BAR, YE;
	int* NET_ISO;
	int* CHEM_ID;
	double* XA;
} Opacity;

void get_sol_x(double*);
void get_sol_chem_id(int*);
void init_mesa();
int get_num_chem_isos();

void init_Opacity(Opacity*);

void shutdown_Opacity(Opacity*);

void eos_PT(Opacity* op, double Pgas, double T,
	double* Rho, double* log10Rho,
	double* dlnRho_dlnPgas_const_T, double* dlnRho_dlnT_const_Pgas, double* res,
	int* ierr);

void eos_DT(Opacity* op, double Rho, double log10Rho, double T, double log10T, double* res,
            double* d_dlnRho_const_T, double* d_dlnT_const_Rho,
            double* Pgas, double* Prad, double* energy, double* entropy, int* ierr);

void kap_DT(Opacity* op, double Rho, double T, double lnfree_e,
	double* kappa, double* dlnkap_dlnRho, double* dlnkap_dlnT, int* ierr);

int nuclide_index(char* nuclei);

#endif  // OPACITY_H
