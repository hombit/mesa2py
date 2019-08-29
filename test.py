#!/usr/bin/env python3

import os

import numpy as np

from opacity import Opac

# opacity = Opac({b'h1': 0.7, b'he4': 0.3})

opacity = Opac('solar')

print(opacity.X)
p = 1e4
temp = 1e4
rho = opacity.rho(p, temp)
print(rho)
rho, eos = opacity.rho(p, temp, True)
kappa = opacity.kappa(rho, temp)
print(kappa)
kappa_with_free_e = opacity.kappa(rho, temp, eos.lnfree_e)
print(kappa_with_free_e)

del opacity
