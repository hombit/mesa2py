import os

import numpy as np 
cimport numpy as cnp

from opacity cimport *


cdef class Opac:
    cdef Opacity fort_opacity

    def __cinit__(self, mesa_dir='mesa'):
        os.environ['MESA_DIR'] = mesa_dir
        self.fort_opacity = init_Opacity()

    def __dealloc__(self):
        shutdown_Opacity(&self.fort_opacity)

    @property
    def X(self):
        return self.fort_opacity.X

    def PT_eos(self, pres, temp, return_grad=False):
        cdef tuple base_shape = cnp.broadcast(pres, temp).shape
        rho = np.empty(base_shape, np.double)
        log10Rho = np.empty(base_shape, np.double)
        dlnRho_dlnPgas_const_T = np.empty(base_shape, np.double)
        dlnRho_dlnT_const_Pgas = np.empty(base_shape, np.double)
        cdef species_double_array res_
        cdef species_double_array d_dlnRho_const_T_
        cdef species_double_array d_dlnT_const_Rho_
        cdef species_double_array d_dabar_const_TRho_
        cdef species_double_array d_dzbar_const_TRho_
        ierr = np.zeros(base_shape, dtype=np.int)
        
        cdef cnp.broadcast it = cnp.broadcast(
            pres, temp, rho, log10Rho,
            dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,
            ierr
        )
        cdef int i
        while cnp.PyArray_MultiIter_NOTDONE(it):
            eos_PT(&self.fort_opacity,
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 0))[0],
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 1))[0],
                   <double*> cnp.PyArray_MultiIter_DATA(it, 2),
                   <double*> cnp.PyArray_MultiIter_DATA(it, 3),
                   <double*> cnp.PyArray_MultiIter_DATA(it, 4),
                   <double*> cnp.PyArray_MultiIter_DATA(it, 5),
                   res_, d_dlnRho_const_T_, d_dlnT_const_Rho_,
                   d_dabar_const_TRho_, d_dzbar_const_TRho_,
                   <int*> cnp.PyArray_MultiIter_DATA(it, 6),
            )
            if (<int*> cnp.PyArray_MultiIter_DATA(it, 6))[0] != 0:
                (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0] = np.nan
                (<double*> cnp.PyArray_MultiIter_DATA(it, 3))[0] = np.nan
                (<double*> cnp.PyArray_MultiIter_DATA(it, 4))[0] = np.nan
                (<double*> cnp.PyArray_MultiIter_DATA(it, 5))[0] = np.nan
            cnp.PyArray_MultiIter_NEXT(it)
        if return_grad:
            return rho, dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas
        return rho

