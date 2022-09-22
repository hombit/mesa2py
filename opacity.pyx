import math as m
import os
import weakref
from collections import OrderedDict, namedtuple

from cython cimport view
from libc.math cimport NAN

import numpy as np 
cimport numpy as cnp

from opacity cimport *


cdef class _Mesa:
    def __cinit__(self, mesa_dir=None):
        if mesa_dir is not None:
            os.environ['MESA_DIR'] = mesa_dir
        init_mesa()

    def __dealloc__(self):
        shutdown_mesa()


class Mesa(_Mesa):
    obj_ref = None

    def __new__(cls, mesa_dir=None):
        if cls.obj_ref is None or cls.obj_ref() is None:
            obj = super().__new__(cls, mesa_dir)
            cls.obj_ref = weakref.ref(obj)
        return cls.obj_ref()


EOSResults = namedtuple(
    'EOSResults', 
    ('dlnRho_dlnPgas_const_T', 'dlnRho_dlnT_const_Pgas',
     'mu', 'lnfree_e', 'grad_ad', 'c_p',),
)


cdef class _Opac:
    cdef object mesa

    cdef Opacity fort_opacity
    default_lnfree_e = <double> 0.

    cdef cnp.ndarray _net_iso
    cdef cnp.ndarray _chem_id
    cdef cnp.ndarray _xa

    def __cinit__(self, composition):
        self.mesa = Mesa()
        
        self._net_iso = np.zeros(get_num_chem_isos(), dtype=np.intc)
        cdef int[:] view_net_iso = self._net_iso
        self.fort_opacity.NET_ISO = &view_net_iso[0]

        self.fort_opacity.SPECIES = SOLSIZE if composition == 'solar' else len(composition)
        self._chem_id = np.empty(self.fort_opacity.SPECIES, dtype=np.intc)
        self._xa = np.zeros(self.fort_opacity.SPECIES, dtype=np.double)
        cdef int[:] view_chem_id = self._chem_id
        cdef double[:] view_xa = self._xa
        self.fort_opacity.CHEM_ID = &view_chem_id[0]
        self.fort_opacity.XA = &view_xa[0]

        if composition == 'solar':
            get_sol_x(&view_xa[0])
            get_sol_chem_id(&view_chem_id[0])

            for i, index in enumerate(self._chem_id):
                self._net_iso[index - 1] = i + 1
        else:
            try:
                composition = dict(composition)
            except ValueError as e:
                raise ValueError('Composition should be "solar" or convertible to dict') from e

            for i, (isotope, num_dens) in enumerate(composition.items()):
                index = nuclide_index(isotope)
                if index < 0: # nuclide_not_found
                    raise ValueError('Invalid isotope name {}'.format(isotope))
                self._chem_id[i] = index
                self._net_iso[index - 1] = i + 1
                self._xa[i] = num_dens

        init_Opacity(&self.fort_opacity)

    def __dealloc__(self):
        shutdown_Opacity(&self.fort_opacity)

    @property
    def X(self):
        return self.fort_opacity.X

    @property
    def Y(self):
        return self.fort_opacity.Y

    @property
    def Z(self):
        return self.fort_opacity.Z

    @property
    def net_iso(self):
        return self._net_iso.copy()

    @property
    def chem_id(self):
        return self._chem_id.copy()

    @property
    def xa(self):
        return self._xa.copy()

    def rho(self, pres, temp, full_output=False):
        """Density from gas pressure and temperature"""
        cdef tuple base_shape = cnp.broadcast(pres, temp).shape
        rho = np.empty(base_shape, np.double)
        log10Rho = np.empty(base_shape, np.double)
        dlnRho_dlnPgas_const_T = np.empty(base_shape, np.double)
        dlnRho_dlnT_const_Pgas = np.empty(base_shape, np.double)
        ierr = np.zeros(base_shape, dtype=np.intc)

        mu = np.empty(base_shape, np.double)
        lnfree_e = np.empty(base_shape, np.double)
        grad_ad = np.empty(base_shape, np.double)
        c_p = np.empty(base_shape, np.double)

        res = view.array(shape=(NUM_EOS_RESULTS,), itemsize=sizeof(double), format='d')
        cdef double[:] res_view = res
        
        cdef cnp.broadcast it = cnp.broadcast(
            pres, temp, rho, log10Rho,
            dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,
            ierr,
            mu, lnfree_e, grad_ad, c_p
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
                (<double*> cnp.PyArray_MultiIter_DATA(it, 10))[0] = res[9]
                if (<int*> cnp.PyArray_MultiIter_DATA(it, 6))[0] != 0:
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 3))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 4))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 5))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 7))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 8))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 9))[0] = NAN
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 10))[0] = NAN
            elif (<int*> cnp.PyArray_MultiIter_DATA(it, 6))[0] != 0:
                (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0] = NAN
            cnp.PyArray_MultiIter_NEXT(it)
        if full_output:
            return rho, EOSResults(dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,
                                   mu, lnfree_e, grad_ad, c_p)
        return rho


    def kappa(self, rho, temp, lnfree_e=default_lnfree_e, return_grad=False):
        cdef tuple base_shape = cnp.broadcast(rho, temp).shape
        kappa = np.empty(base_shape, np.double)
        dlnkap_dlnRho = np.empty(base_shape, np.double)
        dlnkap_dlnT = np.empty(base_shape, np.double)
        ierr = np.zeros(base_shape, np.intc)
        
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


def str2bytes(s):
    if isinstance(s, bytes):
        return s
    if isinstance(s, str):
        return s.encode('ascii')
    raise ValueError('Value ({}) must be str or bytes, not {}'.format(s, type(s)))


def signif(x, p):
    """Inspired by https://stackoverflow.com/a/59888924/5437597"""
    # We know that x is non-negative and finite
    if x == 0:
        return 0
    mags = 10 ** (p - m.floor(m.log10(x)))
    return round(x * mags) / mags


class Opac(_Opac):
    obj_refs = weakref.WeakValueDictionary()

    @staticmethod
    def normalize_composition(c):
        if isinstance(c, str):
            return c
        try:
            d = OrderedDict(c)
        except ValueError as e:
            raise ValueError('Composition must be str or convertible to dict') from e

        for isotope, num_dens in d.items():
            if num_dens >= 0 and not m.isinf(num_dens):
                continue
            raise ValueError(
                f'All composition values must be finite and non-negative, but {isotope} has number density {num_dens}'
            )
        
        bytes_isotopes = map(str2bytes, d)
        
        norm = sum(d.values())
        if norm == 0:
            raise ValueError('At least one composition value must be positive')
        normalized_dens = (num_dens / norm for num_dens in d.values())
        rounded_dens = (signif(x, 12) for x in normalized_dens)
        
        t = tuple((isotope, num_dens) for isotope, num_dens in zip(bytes_isotopes, rounded_dens) if num_dens > 0)
        return t

    def __new__(cls, composition):
        composition = cls.normalize_composition(composition)

        if composition not in cls.obj_refs:
            obj = super().__new__(cls, composition)
            cls.obj_refs[composition] = obj
        
        return cls.obj_refs[composition]

