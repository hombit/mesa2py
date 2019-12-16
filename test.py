#!/usr/bin/env python3

import os

import numpy as np

from opacity import Opac

opacity = Opac({b'h1': 1.0})

# opacity = Opac('solar')

print(opacity.X)
p = 1e5
temp = 1e5
rho = opacity.rho(p, temp)

# energy, Pgas, Prad, entropy = opacity.energy(rho, temp)
# print('Energy = ', energy, Pgas, Prad, entropy)
print(rho)
rho, eos = opacity.rho(p, temp, True)
print(eos)
kappa = opacity.kappa(rho, temp)
print(kappa)
kappa_with_free_e = opacity.kappa(rho, temp, eos.lnfree_e)
print(kappa_with_free_e)

del opacity
