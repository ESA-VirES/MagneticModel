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

### Conda installation

The package contains the `conda-buidl` scripts allowing local conda build and
installation following this procedure:

```
conda activate <target-environment>
conda install numpy scipy matplotlib h5py networkx cython conda-build
conda-build ./qdipole
conda-build ./cdf
conda-build purge
conda install --use-local qdipole cdf
pip install spacepy
pip install eoxmagmod
```
