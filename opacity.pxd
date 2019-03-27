# cython: language_level=3

ctypedef double[7] species_double_array

cdef extern from 'opacity.h' nogil:
    cdef int _SPECIES

    cdef int NUM_EOS_RESULTS

    ctypedef struct Opacity:
        int EOS_HANDLER
        int KAP_HANDLER
        species_double_array XA
        double X
        double Y
        double Z
        double ABAR
        double ZBAR
        double Z2BAR
        double YE
        int* NET_ISO
        int* CHEM_ID

    cdef void init_Opacity(Opacity*)
    cdef void shutdown_Opacity(Opacity*)
    cdef void eos_PT(Opacity*, double, double,
                     double*, double*, double*, double*,
                     double*,
                     int*)
    cdef void kap_DT(Opacity*, double, double, double,
                     double*, double*, double*, int*)

    cdef int nuclide_index(char*)