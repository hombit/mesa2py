import math as m
import os
import re
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

        self.eos_result_dtype = np.dtype([(name, np.double) for name in self.eos_result_names])

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
    def eos_result_names(self):
        name_array = view.array(shape=(EOS_NAME_LENGTH,), itemsize=sizeof(char), format='c')
        cdef char[:] name_view = name_array

        cdef int i
        cdef bytes name_bytes

        # Names for res array components
        res_names = []
        for i in range(NUM_EOS_RESULTS):
            get_eosDT_result_name(i, &name_view[0])
            untrimmed_name = bytes(name_array).decode()
            # Trim whitespaces and remove stuff like "/"
            name = re.sub(r'\W', '', untrimmed_name)
            res_names.append(name)

        names = (
            res_names
            + [f'd_{name}_dlnRho_const_T' for name in res_names]
            + [f'd_{name}_dlnT_const_Rho' for name in res_names]
        )

        return names

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
        dlnRho_dlnPgas_const_T = np.empty(base_shape, np.double)
        dlnRho_dlnT_const_Pgas = np.empty(base_shape, np.double)

        # Various thermodynamic quantities and there derivatives
        eos_result = np.recarray(base_shape, dtype=self.eos_result_dtype)

        cdef int ierr = 0

        cdef int i_eos_result

        cdef cnp.broadcast it = cnp.broadcast(
            pres, temp, rho,
            dlnRho_dlnPgas_const_T, dlnRho_dlnT_const_Pgas,
            eos_result,
        )
        while cnp.PyArray_MultiIter_NOTDONE(it):
            eos_PT(&self.fort_opacity,
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 0))[0],                     # Pgas
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 1))[0],                     # T
                   <double*> cnp.PyArray_MultiIter_DATA(it, 2),                          # Rho
                   <double*> cnp.PyArray_MultiIter_DATA(it, 3),                          # dlnRho_dlnPgas_const_T
                   <double*> cnp.PyArray_MultiIter_DATA(it, 4),                          # dlnRho_dlnT_const_Pgas
                   <double*> cnp.PyArray_MultiIter_DATA(it, 5),                          # res
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 5)) + NUM_EOS_RESULTS,      # d_dlnRho_const_T
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 5)) + 2 * NUM_EOS_RESULTS,  # d_dlnT_const_Rho
                   &ierr)
            if ierr != 0:
                (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0] = NAN                 # Rho
                (<double*> cnp.PyArray_MultiIter_DATA(it, 3))[0] = NAN                 # dlnRho_dlnPgas_const_T
                (<double*> cnp.PyArray_MultiIter_DATA(it, 4))[0] = NAN                 # dlnRho_dlnT_const_Pgas
                for i_eos_result in range(3 * NUM_EOS_RESULTS):
                    (<double*> cnp.PyArray_MultiIter_DATA(it, 5))[i_eos_result] = NAN  # res and derivatives
            ierr = 0
            cnp.PyArray_MultiIter_NEXT(it)
        if full_output:
            return rho, eos_result
        return rho

    def kappa(self, rho, temp, eos=None, return_grad=False):
        """Opacity from density and temperature"""
        cdef tuple base_shape
        cdef cnp.ndarray lnfree_e, d_lnfree_e_d_lnRho, d_lnfree_e_d_lnT, eta, d_eta_d_lnRho, d_eta_d_lnT
        if eos is None:
            base_shape = cnp.broadcast(eos).shape

            lnfree_e = np.zeros(base_shape, np.double)
            d_lnfree_e_d_lnRho = np.zeros(base_shape, np.double)
            d_lnfree_e_d_lnT = np.zeros(base_shape, np.double)
            eta = np.zeros(base_shape, np.double)
            d_eta_d_lnRho = np.zeros(base_shape, np.double)
            d_eta_d_lnT = np.zeros(base_shape, np.double)
        else:
            lnfree_e = eos.lnfree_e
            d_lnfree_e_d_lnRho = eos.d_lnfree_e_dlnRho_const_T
            d_lnfree_e_d_lnT = eos.d_lnfree_e_dlnT_const_Rho
            eta = eos.eta
            d_eta_d_lnRho = eos.d_eta_dlnRho_const_T
            d_eta_d_lnT = eos.d_eta_dlnT_const_Rho

            base_shape = cnp.broadcast(rho, temp, lnfree_e).shape

        kappa = np.empty(base_shape, np.double)
        dlnkap_dlnRho = np.empty(base_shape, np.double)
        dlnkap_dlnT = np.empty(base_shape, np.double)

        cdef int ierr = 0

        cdef cnp.broadcast it = cnp.broadcast(
            rho, temp,
            lnfree_e, d_lnfree_e_d_lnRho, d_lnfree_e_d_lnT,
            eta, d_eta_d_lnRho, d_eta_d_lnT,
            kappa, dlnkap_dlnRho, dlnkap_dlnT,
        )
        while cnp.PyArray_MultiIter_NOTDONE(it):
            kap_DT(&self.fort_opacity,
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 0))[0],  # Rho
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 1))[0],  # T
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 2))[0],  # lnfree_e
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 3))[0],  # d_lnfree_e_d_lnRho
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 4))[0],  # d_lnfree_e_d_lnT
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 5))[0],  # eta
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 6))[0],  # d_eta_d_lnRho
                   (<double*> cnp.PyArray_MultiIter_DATA(it, 7))[0],  # d_eta_d_lnT
                   <double*> cnp.PyArray_MultiIter_DATA(it, 8),       # kappa
                   <double*> cnp.PyArray_MultiIter_DATA(it, 9),       # d_lnkap_d_lnRho
                   <double*> cnp.PyArray_MultiIter_DATA(it, 10),      # d_lnkap_d_lnT
                   &ierr)
            if ierr != 0:
                (<double*> cnp.PyArray_MultiIter_DATA(it, 8))[0] = NAN   # kappa
                (<double*> cnp.PyArray_MultiIter_DATA(it, 9))[0] = NAN   # d_lnkap_d_lnRho
                (<double*> cnp.PyArray_MultiIter_DATA(it, 10))[0] = NAN  # d_lnkap_d_lnT
            ierr = 0
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
