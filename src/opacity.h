#ifndef OPACITY_H
#define OPACITY_H

extern int NUM_EOS_RESULTS;
extern int NUM_CHEM_ISOS_POINTER;
extern int SOLSIZE;
extern int EOS_NAME_LENGTH;

typedef struct {
	int EOS_HANDLER, KAP_HANDLER;
	int SPECIES;
	double X, Y, Z, XC, XN, XO, XNe, ABAR, ZBAR, Z2BAR, Z53BAR, YE;
	int* NET_ISO;
	int* CHEM_ID;
	double* XA;
} Opacity;

extern void get_sol_x(double*);
extern void get_sol_chem_id(int*);
extern int get_num_chem_isos();

extern void get_eosDT_result_name(int index, char* name);

extern void init_mesa();
extern void shutdown_mesa();

extern void init_Opacity(Opacity*);
extern void shutdown_Opacity(Opacity*);
extern void eos_PT(Opacity* op, double Pgas, double T,
	double* Rho, double* dlnRho_dlnPgas_const_T, double* dlnRho_dlnT_const_Pgas,
	double* res, double* d_dlnRho_const_T, double *d_dlnT_const_Rho,
	int* ierr);
extern void kap_DT(Opacity* op, double Rho, double T,
    double lnfree_e, double d_lnfree_e_dlnRho, double d_lnfree_e_dlnT,
    double eta, double d_eta_dlnRho, double d_eta_dlnT,
	double* kappa, double* dlnkap_dlnRho, double* dlnkap_dlnT,
	int* ierr);

extern int nuclide_index(char* nuclei);

#endif  // OPACITY_H
