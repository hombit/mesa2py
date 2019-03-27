import os
from collections import namedtuple

from cython cimport view
from libc.math cimport NAN

import numpy as np 
cimport numpy as cnp

from opacity cimport *


EOSResults = namedtuple(
    'EOSResults', 
    ('dlnRho_dlnPgas_const_T', 'dlnRho_dlnT_const_Pgas',
     'mu', 'lnfree_e', 'grad_ad',),
)


cdef class Opac:
    cdef Opacity fort_opacity
    default_lnfree_e = <double> 0.

    cdef cnp.ndarray net_iso
    cdef cnp.ndarray chem_id
    cdef cnp.ndarray xa

    def __cinit__(self, composition, mesa_dir=None):
        if mesa_dir is not None:
            os.environ['MESA_DIR'] = mesa_dir

        init_mesa()

        self.fort_opacity.SPECIES = len(composition)
        self.net_iso = np.zeros(get_num_chem_isos(), dtype=np.intc)
        cdef int[:] view_net_iso = self.net_iso
        self.fort_opacity.NET_ISO = &view_net_iso[0]
        self.chem_id = np.empty(self.fort_opacity.SPECIES, dtype=np.intc)
        cdef int[:] view_chem_id = self.chem_id
        self.fort_opacity.CHEM_ID = &view_chem_id[0]
        self.xa = np.zeros(self.fort_opacity.SPECIES, dtype=np.double)
        cdef double[:] view_xa = self.xa
        self.fort_opacity.XA = &view_xa[0]

        norm = sum(composition.values())
        composition = {isotope: num_dens / norm for isotope, num_dens in composition.items()}

        for i, (isotope, num_dens) in enumerate(composition.items()):
            index = nuclide_index(isotope)
            if index < 0: # nuclide_not_found
                raise ValueError('Invalid isotope name {}'.format(isotope))
            self.chem_id[i] = index
            self.net_iso[index - 1] = i + 1
            self.xa[i] = num_dens

        init_Opacity(&self.fort_opacity)

    def __dealloc__(self):
        shutdown_Opacity(&self.fort_opacity)

    @property
    def X(self):
        return self.fort_opacity.X

    def rho(self, pres, temp, full_output=False):
        cdef tuple base_shape = cnp.broadcast(pres, temp).shape
        rho = np.empty(base_shape, np.double)
        log10Rho = np.empty(base_shape, np.double)
        dlnRho_dlnPgas_const_T = np.empty(base_shape, np.double)
        dlnRho_dlnT_const_Pgas = np.empty(base_shape, np.double)
        ierr = np.zeros(base_shape, dtype=np.int)

        mu = np.empty(base_shape, np.double)
        lnfree_e = np.empty(base_shape, np.double)
        grad_ad = np.empty(base_shape, np.double)

        res = view.array(shape=(NUM_EOS_RESULTS,), itemsize=sizeof(double), format='d')
        cdef double[:] res_view = res
        
        cdef cnp.broadcast it = cnp.broadcast(
            pres, temp, rho, log10Rho,
            dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,
            ierr,
            mu, lnfree_e, grad_ad,
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
                   &res_view[0],
                   <int*> cnp.PyArray_MultiIter_DATA(it, 6))
            if full_output:
                # These indexes are grabbed from eos_def
                (<double*> cnp.PyArray_MultiIter_DATA(it, 7))[0] = res[3]
                (<double*> cnp.PyArray_MultiIter_DATA(it, 8))[0] = res[4]
                (<double*> cnp.PyArray_MultiIter_DATA(it, 9))[0] = res[6]
                if (<int*> cnp.PyArray_MultiIter_DATA(it, 6))[0] != 0:
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 3))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 4))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 5))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 7))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 8))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 9))[0] = NAN
            elif (<int*> cnp.PyArray_MultiIter_DATA(it, 6))[0] != 0:
                (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0] = NAN
            cnp.PyArray_MultiIter_NEXT(it)
        if full_output:
            return rho, EOSResults(dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,
                                   mu, lnfree_e, grad_ad,)
        return rho

    def kappa(self, rho, temp, lnfree_e=default_lnfree_e, return_grad=False):
        cdef tuple base_shape = cnp.broadcast(rho, temp).shape
        kappa = np.empty(base_shape, np.double)
        dlnkap_dlnRho = np.empty(base_shape, np.double)
        dlnkap_dlnT = np.empty(base_shape, np.double)
        ierr = np.zeros(base_shape, np.int)
        
        cdef cnp.broadcast it = cnp.broadcast(
            rho, temp, lnfree_e,
            kappa, dlnkap_dlnRho, dlnkap_dlnT,
            ierr
        )
        cdef int i
        while cnp.PyArray_MultiIter_NOTDONE(it):
            kap_DT(&self.fort_opacity,
                        (<double*> cnp.PyArray_MultiIter_DATA(it, 0))[0],
                        (<double*> cnp.PyArray_MultiIter_DATA(it, 1))[0],
                        (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0],
                        <double*> cnp.PyArray_MultiIter_DATA(it, 3),
                        <double*> cnp.PyArray_MultiIter_DATA(it, 4),
                        <double*> cnp.PyArray_MultiIter_DATA(it, 5),
                        <int*> cnp.PyArray_MultiIter_DATA(it, 6))
            if (<int*> cnp.PyArray_MultiIter_DATA(it, 6))[0] != 0:
                (<double*> cnp.PyArray_MultiIter_DATA(it, 3))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 4))[0] = NAN
                (<double*> cnp.PyArray_MultiIter_DATA(it, 5))[0] = NAN
            cnp.PyArray_MultiIter_NEXT(it)
        if return_grad:
            return kappa, dlnkap_dlnRho, dlnkap_dlnT
        return kappa

