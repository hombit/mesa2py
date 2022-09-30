#!/usr/bin/env python3

import os
import unittest

import numpy as np
from numpy.testing import assert_allclose, assert_array_equal
from parameterized import parameterized, param

from opacity import Mesa, Opac


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
        t = 1e6
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
                rho=7.749773e-10,
                kappa=23.521732,
                eos={'mu': 0.643923, 'Cp': 5.7465e+09, 'lnfree_e': -0.57849, 'grad_ad': 0.07727},
            ),
            param(
                composition={'he4': 1.0},
                p=1e3,
                t=3e4,
                rho=7.761565e-10,
                kappa=0.227659,
                eos={'mu': 1.935747, 'Cp': 6.01333e+09, 'lnfree_e': -1.32141, 'grad_ad': 0.219777},
            ),
            param(
                composition='solar',
                p=1e3,
                t=1e4,
                rho=1.02442e-09,
                kappa=18.323045,
                eos={'mu': 0.85119, 'Cp': 4.055806e+09, 'lnfree_e': -0.90484, 'grad_ad': 0.078819},
            ),
            param(
                composition='solar',
                p=1e3,
                t=10**3.84,  # inside lowT and normal opacity tables blending interval
                rho=2.229603e-09,
                kappa=0.265253,
                eos={'mu': 1.282465, 'Cp': 4.215608e+08, 'lnfree_e': -4.651135, 'grad_ad': 0.181976},
            ),
        ]
    )
    def test_regression(self, composition, p, t, rho, kappa, eos):
        op = Opac(composition)

        with self.subTest('rho'):
            assert_allclose(op.rho(p, t), rho, rtol=1e-5, err_msg='rho')

        with self.subTest('kappa'):
            assert_allclose(op.kappa(rho, t), kappa, rtol=1e-5, err_msg='kappa')

        with self.subTest('kappa w/ vs w/o eos'):
            rho_actual, eos_actual = op.rho(p, t, full_output=True)
            assert_allclose(op.kappa(rho_actual, t), op.kappa(rho_actual, t, eos=eos_actual), rtol=1e-5,
                            err_msg='kappa w/ vs w/o eos')

        _, eos_actual = op.rho(p, t, full_output=True)
        for field, value in eos.items():
            subtest_name = f'eos.{field}'
            with self.subTest(subtest_name):
                assert_allclose(getattr(eos_actual, field), value, equal_nan=False, rtol=1e-5, err_msg=subtest_name)
