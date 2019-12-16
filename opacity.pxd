# cython: language_level=3

cdef extern from 'opacity.h' nogil:

    cdef int NUM_EOS_RESULTS
    cdef int NUM_CHEM_ISOS_POINTER
    cdef int SOLSIZE

    ctypedef struct Opacity:
        int EOS_HANDLER
        int KAP_HANDLER
        int SPECIES
        double X
        double Y
        double Z
        double XC
        double XN
        double XO
        double XNe
        double ABAR
        double ZBAR
        double Z2BAR
        double YE
        int* NET_ISO
        int* CHEM_ID
        double* XA

    cdef void get_sol_x(double*)
    cdef void get_sol_chem_id(int*)
    cdef void init_mesa()
    cdef int get_num_chem_isos()
    cdef void init_Opacity(Opacity*)
    cdef void shutdown_Opacity(Opacity*)
    cdef void eos_PT(Opacity*, double, double,
                     double*, double*, double*, double*,
                     double*,
                     int*)
    cdef void kap_DT(Opacity*, double, double, double,
                     double*, double*, double*, int*)

    cdef int nuclide_index(char*)