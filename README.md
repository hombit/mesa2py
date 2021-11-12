# Yet another Python3 binding to MESA

It covers only `eos` and `kappa` modules from [MESA code](http://mesa.sourceforge.net). 

### Installation using Docker

You can use the latest pre-build Docker image:

```shell
$ docker run --rm -ti ghcr.io/hombit/mesa2py
```

Or build a docker image by yourself

``` shell
$ git clone https://github.com/hombit/mesa2py
$ cd mesa2py
$ docker build -t mesa2py .
```

### Usage

`mesa2py` contains `opacity` module. Run docker 'mesa2py' image as a container and run Python3 in it

``` shell
$ docker run --rm -ti mesa2py python3
```

Then import `Opac` class from the `opacity` module

``` python3
from opacity import Opac
composition = 'solar'  # solar chemical composition
mesaop = Opac(composition)
```

`Opac` class constructor requires one argument - chemical composition of matter. 
It should be a dictionary with format {'isotope_name': abundance}, e.g. {'h1': 0.7, 'he4': 0.3}, look for full list of available isotopes in the Mesa source code. Also you can use 'solar' string to set the solar composition.

`Opac` object has two methods: `rho` and `kappa`. `rho` calculates bulk density from pressure and temperature, using MESA EOS tables. Usage of `rho`:

``` python3
>>> P = 2e4  # pressure in cgs units
>>> T = 3e5  # temperature in K
>>> rho = mesaop.rho(P, T)  # bulk density in g/cm^3
>>> print(rho)
4.92741807e-10
```

You can also use multidimensional `numpy` arrays as an input, usual `numpy` broadcasting rules will apply, for example
```python3
>>> import numpy as np

>>> P = np.logspace(3, 5, 11).reshape(-1, 1)
>>> T = np.logspace(4, 6, 11).reshape(1, -1)
>>> rho = mesaop.rho(P, T)
>>> print(rho.shape)
(11, 11)
```

If a keyword argument 'full_output' is True, `rho` also returns 'eos' namedtuple, that contains several thermodynamic functions.

``` python3
>>> rho, eos = mesaop.rho(P, T, full_output=True)
>>> print(rho)
4.927418070544516e-10
>>> print(eos)
EOSResults(dlnRho_dlnPgas_const_T=array(1.01938468), dlnRho_dlnT_const_Pgas=array(-0.98225163), mu=array(0.61375111), lnfree_e=array(-0.16125285), grad_ad=array(0.25), c_p=array(2.25940722e+15))
```

`kappa` method calculates opacity coefficient from bulk density and temperature, using MESA opacity tables. Usage of `kappa`:

``` python3
>>> rho = 5e-10  # bulk density in g/cm^3
>>> T = 3e5  # temperature in K
>>> kappa = mesaop.kappa(rho, T)  # opacity coefficient in cm^2/g, opacity tables from MESA code are used
>>> print(kappa)
0.35719242
```
