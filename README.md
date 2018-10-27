![Travis-CI staus](https://api.travis-ci.org/ESA-VirES/MagneticModel.svg?branch=master)


## Magnetic Model

This repository contains various utilities related to Earth magnetic field
modelling and spherical harmonics.

The repository contains following directories:

- `eoxmagmod` - Collection models of the Earth magnetic field - python module
- `qdipole` - Quasi-Dipole apex coordinates evaluation - Fortran code compiled
  as a shared library (dependency of the `eoxmagmod` package)
- `libcdf` - [CDF library](https://cdf.gsfc.nasa.gov/) source installation
  (dependency of the `eoxmagmod` package)

### Installation from Sources

#### CDF

```
$ cd libcdf/
$ make build
$ sudo make install
```

#### QDIPOLE

```
$ cd qdipole/
$ ./configure
$ make build
$ sudo make install
```

#### EOxMagMod
Requires QDIPOLE, CDF libraries + NumPy and SpacePy Python packages
to be installed.
NumPy and SpacePy can be installed using `pip`.

```
$ cd eoxmagmod/
$ python ./setup.py build
$ sudo python ./setup.py install
```
