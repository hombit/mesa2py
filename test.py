#!/usr/bin/env python3

import os

import numpy as np

from opacity import Opac

opacity = Opac(mesa_dir='/mesa')

print(opacity.X)
p = 1e4
temp = 1e4
rho = opacity.rho(p, temp)
print(rho)
print(opacity.rho(p, temp, True))
kappa = opacity.kappa(rho, temp)
print(kappa)

del opacity
