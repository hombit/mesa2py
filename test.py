#!/usr/bin/env python3

import os
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from parameterized import parameterized, param

from opacity import EOSResults, Mesa, Opac


class MesaUnitTestCase(unittest.TestCase):
    def test_singlton(self):
        m1 = Mesa()
        m2 = Mesa()
        assert m1 is m2


class OpacUnitTestCase(unittest.TestCase):
    def opac_equal(self, op1, op2):
        assert_array_equal(op1.net_iso, op2.net_iso)
        assert_array_equal(op1.chem_id, op2.chem_id)
        assert_allclose(op1.xa, op2.xa)

    def test_bytes_str_composition_keys(self):
        op_str = Opac({'h1': 0.7, 'he4': 0.3})
        op_bytes = Opac({b'h1': 0.7, b'he4': 0.3})
        assert op_str is op_bytes

    def test_composition_norm(self):
        op_norm = Opac({'h1': 0.8, 'he4': 0.2})
        op_unnorm = Opac({'h1': 1.0, 'he4': 0.25})
        assert op_norm is op_unnorm
    
    def test_hydrogen(self):
        op = Opac({'h1': 1.0})
        with self.subTest('X'):
            self.assertEqual(op.X, 1.0)
        
        p = 1e3
        t = 1e7
        rho, eos = op.rho(p, t, True)
        with self.subTest('rho'):
            assert_allclose(rho, p * eos.mu / (8.31e7 * t), rtol=1e-2)
        with self.subTest('mu'):
            assert_allclose(eos.mu, 0.5, rtol=1e-3)

    def test_negative_dens(self):
        with self.assertRaises(ValueError):
            Opac({'h1': 1.0, 'he4': -0.1})

    def test_all_zero_dens(self):
        with self.assertRaises(ValueError):
            Opac({'h1': 0.0, 'he4': 0.0})

    def test_inf_zero_dens(self):
        with self.assertRaises(ValueError):
            Opac({'h1': 1.0, 'he4': np.inf})

    def test_nan_zero_dens(self):
        with self.assertRaises(ValueError):
            Opac({'h1': 1.0, 'he4': np.nan})

    def test_cache(self):
        n = 1024
        ops = [Opac('solar') for _ in range(n)]
        assert all(ops[0] is op for op in ops)

    @parameterized.expand(
        [
            param(
                composition={'h1': 1.0},
            ),
            param(
                composition={'he4': 1.0},
            ),
            param(
                composition={'h1': 1.0, 'he4': 1.0, 'c12': 0.1},
            ),
            param(
                composition='solar',
            ),
        ]
    )
    def test_XYZ(self, composition):
        opac = Opac(composition)
        assert_allclose(opac.X + opac.Y + opac.Z, 1.0, atol=1e-12, rtol=0)

    @parameterized.expand(
        [
            param(
                composition={'h1': 1.0},
                p=1e3,
                t=1e4,
                rho=7.75826e-10,
                kappa=23.521732,
                eos=EOSResults(1.122646, -3.244283, 0.644622, -0.57993, 0.077202, 5.751599e+09),
            ),
            param(
                composition={'he4': 1.0},
                p=1e3,
                t=3e4,
                rho=7.76229824e-10,
                kappa=0.22766244,
                eos=EOSResults(1.02917524, -1.69374361, 1.93590453, -1.31933614, 0.22034105, 5.9938266e+09),
            ),
            param(
                composition='solar',
                p=1e3,
                t=1e4,
                rho=1.024002e-09,
                kappa=18.345603,
                eos=EOSResults(1.114417, -3.092285, 0.847222, -0.902434, 0.078645, 4.060632e+09),
            ),
            param(
                composition='solar',
                p=1e3,
                t=10**3.84,  # inside lowT and normal opacity tables blending interval
                rho=2.22718465e-09,
                kappa=0.26554388,
                eos=EOSResults(1.00607563, -1.15243794, 1.2775947, -4.65811899, 0.18233674, 4.20004691e+08),
            ),
        ]
    )
    def test_regression(self, composition, p, t, rho, kappa, eos):
        op = Opac(composition)
        with self.subTest('rho'):
            assert_allclose(op.rho(p, t), rho, rtol=1e-5)
        with self.subTest('kappa'):
            assert_allclose(op.kappa(rho, t), kappa, rtol=1e-5)
        _, eos_actual = op.rho(p, t, True)
        for field in eos._fields:
            with self.subTest('eos-{}'.format(field)):
                assert_allclose(getattr(eos_actual, field), getattr(eos, field), rtol=1e-5)

