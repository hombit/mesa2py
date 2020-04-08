#!/usr/bin/env python3

import os
from time import sleep

import numpy as np

from opacity import Mesa, Opac

ops = [
    Opac({b'h1': 0.7, b'he4': 0.3}),
    Opac('solar'),
]

for opacity in ops:
    print(f'X: {opacity.X}')
    p = 1e4
    temp = 1e5
    rho = opacity.rho(p, temp)
    print(f'rho = {rho}')
    rho, eos = opacity.rho(p, temp, True)
    print(f'eos = {eos}')
    kappa = opacity.kappa(rho, temp)
    print(f'kappa = {kappa}')
    kappa_with_free_e = opacity.kappa(rho, temp, eos.lnfree_e)
    print(f'kappa with free e = {kappa_with_free_e}')

del ops[1]
del ops[0]
