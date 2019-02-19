import os
from collections import namedtuple

from libc.math cimport NAN

import numpy as np 
cimport numpy as cnp

from opacity cimport *


EOSGradients = namedtuple(
    'EOSGradients', 
    ('dlnRho_dlnPgas_const_T', 'dlnRho_dlnT_const_Pgas',
     'gamma1', 'gamma3'),
)


cdef class Opac:
    cdef Opacity fort_opacity

    def __cinit__(self, mesa_dir=None):
        if mesa_dir is not None:
            os.environ['MESA_DIR'] = mesa_dir
        self.fort_opacity = init_Opacity()

    def __dealloc__(self):
        shutdown_Opacity(&self.fort_opacity)

    @property
    def X(self):
        return self.fort_opacity.X

    def rho(self, pres, temp, return_grad=False):
        cdef tuple base_shape = cnp.broadcast(pres, temp).shape
        rho = np.empty(base_shape, np.double)
        log10Rho = np.empty(base_shape, np.double)
        dlnRho_dlnPgas_const_T = np.empty(base_shape, np.double)
        dlnRho_dlnT_const_Pgas = np.empty(base_shape, np.double)
        gamma1 = np.empty(base_shape, np.double)
        gamma3 = np.empty(base_shape, np.double)
        ierr = np.zeros(base_shape, dtype=np.int)
        
        cdef cnp.broadcast it = cnp.broadcast(
            pres, temp, rho, log10Rho,
            dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,
            gamma1, gamma3,
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
                   <double*> cnp.PyArray_MultiIter_DATA(it, 6),
                   <double*> cnp.PyArray_MultiIter_DATA(it, 7),
                   <int*> cnp.PyArray_MultiIter_DATA(it, 8))
            if (<int*> cnp.PyArray_MultiIter_DATA(it, 8))[0] != 0:
                (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 3))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 4))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 5))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 6))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 7))[0] = NAN
            cnp.PyArray_MultiIter_NEXT(it)
        if return_grad:
            grads = EOSGradients(dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,
                                 gamma1, gamma3)
            return rho, grads
        return rho

    def kappa(self, rho, temp, return_grad=False):
        cdef tuple base_shape = cnp.broadcast(rho, temp).shape
        kappa = np.empty(base_shape, np.double)
        dlnkap_dlnRho = np.empty(base_shape, np.double)
        dlnkap_dlnT = np.empty(base_shape, np.double)
        ierr = np.zeros(base_shape, np.int)
        
        cdef cnp.broadcast it = cnp.broadcast(
            rho, temp,
            kappa, dlnkap_dlnRho, dlnkap_dlnT,
            ierr
        )
        cdef int i
        while cnp.PyArray_MultiIter_NOTDONE(it):
            kap_DT(&self.fort_opacity,
                        (<double*> cnp.PyArray_MultiIter_DATA(it, 0))[0],
                        (<double*> cnp.PyArray_MultiIter_DATA(it, 1))[0],
                        <double*> cnp.PyArray_MultiIter_DATA(it, 2),
                        <double*> cnp.PyArray_MultiIter_DATA(it, 3),
                        <double*> cnp.PyArray_MultiIter_DATA(it, 4),
                        <int*> cnp.PyArray_MultiIter_DATA(it, 5))
            if (<int*> cnp.PyArray_MultiIter_DATA(it, 5))[0] != 0:
                (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 3))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 4))[0] = NAN
            cnp.PyArray_MultiIter_NEXT(it)
        if return_grad:
            return kappa, dlnkap_dlnRho, dlnkap_dlnT
        return kappa
