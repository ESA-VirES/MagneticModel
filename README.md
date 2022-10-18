

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

By default, the library gets installed in `/usr/local/cdf` directory.
To install it to a different path override the `INSTALLDIR`
variable:
```
$ make install INSTALLDIR=<install directory>
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

The package contains the `conda-build` scripts allowing local conda build and
installation following this procedure:

1) build the binary dependencies:
```
conda install conda-build
conda build ./qdipole
conda build ./libcdf
conda build purge
```
Tested on GNU/Linux. Possibly works on other POSIX systems. Does not work on MS
Windows (primarily because of a missing Fortran compiler).

2) install the `eoxmagmod` in your conda environment:
```
conda activate <target-environment>
conda install --use-local qdipole cdf
conda install numpy scipy matplotlib h5py networkx
conda install gcc_linux-64           # spacepy and eoxmagmod require C compiler
conda install gfortran_linux-64      # spacepy requires Fortran compiler
conda deactivate
conda activate <target-environment>  # re-activation is required to update the environment variables
pip install spacepy
pip install ./eoxmagmod
```

The `gfortran_linux-64` and`gcc_linux-64` compilers work on a x86_64 GNU/Linux system.
Other platforms might provide different compilers.
