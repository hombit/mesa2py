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

    def test_cache(self):
        n = 1024
        ops = [Opac('solar') for _ in range(n)]
        assert all(ops[0] is op for op in ops)

    @parameterized.expand(
        [
            param(
                composition={'h1': 1.0},
                p=1e3,
                t=1e4,
                rho=7.75826e-10,
                kappa=23.51589,
                eos=EOSResults(1.122646, -3.244283, 0.644622, -0.57993, 0.077202, 5.751599e+09),
            ),
            param(
                composition={'he4': 1.0},
                p=1e3,
                t=1e4,
                rho=4.802586e-09,
                kappa=0.0007247588733323767,
                eos=EOSResults(1.001195, -1.036848, 3.993066, -7.424327, 0.28037, 86625632.0),
            ),
            param(
                composition='solar',
                p=1e3,
                t=1e4,
                rho=1.024002e-09,
                kappa=18.259231,
                eos=EOSResults(1.114417, -3.092285, 0.847222, -0.902434, 0.078645, 4.060632e+09),
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

