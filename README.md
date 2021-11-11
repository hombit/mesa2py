# Yet another Python3 binding to MESA

It covers only `eos` and `kappa` modules from MESA code. 

### Installation using Docker

Just build docker image

``` shell
$ cd ~/mesa2py
$ docker build -t mesa2py .
```

### Usage

mesa2py contains `opacity` module. Run docker 'mesa2py' image as a container and run Python3 in it

``` shell
$ docker run --rm -ti mesa2py python3
```

Then import `Opac` object from the `opacity` module

``` python3
from opacity import Opac
composition = 'solar'  # solar chemical composition
mesaop = Opac(composition)
```

`Opac` object requires one argument - chemical composition of matter. 

It should be dictionary with format {'isotope_name': abundance}, e.g. {'h1': 0.7, 'he4': 0.3}. Use 'solar' string in case of solar composition.

`Opac` object has two methods - `rho` and `kappa`. `rho` calculates bulk density from pressure and temperature, using MESA EOS tables. Usage of `rho`:

``` python3
P = 2e4  # pressure in cgs units
T = 3e5  # temperature in K
rho = mesaop.rho(P, T)  # bulk density in g/cm^3
print(rho)
```
	4.92741807e-10
If keyword argument 'full_output' is True, `rho` also returns 'eos' namedtuple, that contains several thermodynamic functions.

``` python3
rho, eos = mesaop.rho(P, T, full_output=True)
print(rho)
print(eos)
```

	4.927418070544516e-10
	EOSResults(dlnRho_dlnPgas_const_T=array(1.01938468), dlnRho_dlnT_const_Pgas=array(-0.98225163), mu=array(0.61375111), lnfree_e=array(-0.16125285), grad_ad=array(0.25), c_p=array(2.25940722e+15))


`kappa` calculates opacity coefficient from bulk density and temperature, using MESA opacity tables. Usage of `kappa`:

``` python3
rho = 5e-10  # bulk density in g/cm^3
T = 3e5  # temperature in K
kappa = mesaop.kappa(rho, T)  # opacity coefficient in cm^2/g, opacity tables from MESA code are used
print(kappa)
```
	0.35719242
