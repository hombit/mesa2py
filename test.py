#!/usr/bin/env python3

import os

import numpy as np

from opacity import opacity

os.environ['MESA_DIR'] = os.path.abspath('./mesa')

opacity.init_opacity()
print(opacity.get_hbar())
print(opacity.xa)
opacity.shutdown_opacity()
