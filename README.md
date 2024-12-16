## VirES for Swarm - Magnetic Models, Coordinates and Other Utilities

This repository contains various utilities for calculation of Earth magnetic
field models, magnetic coordinates and other auxiliary variables.

The repository contains following directories:

- `eoxmagmod` - Python package for calculation of the Earth magnetic field,
  magnetic coordinates, time conversions and other auxiliary variables.
- `qdipole` - Quasi-Dipole apex coordinates evaluation - Fortran code compiled
  as a shared library (dependency of the `eoxmagmod` package)
- `libcdf` - [CDF library](https://cdf.gsfc.nasa.gov/) source installation
  (dependency of the `eoxmagmod` package)

For more details read the [on-line documentation](https://esa-vires.github.io/MagneticModel/)

### Installation

#### Conda Installation

Step in the MagneticModel directory and follow these steps:

1) Build the binary dependencies:
```
conda install conda-build
conda build ./qdipole
conda build ./libcdf
conda build purge
```
Tested on GNU/Linux. Possibly works on other POSIX systems.

2) Create an new Conda environment and install `eoxmagmod` dependencies

```
conda env create -n <environment-name> -f conda_env.yaml
```

3) Install `eoxmagmod` in your new Conda environment:
```
conda activate <environment-name>
pip install ./eoxmagmod
```

#### Installation from Sources

The installation requires:
- GCC C compiler
- GFortran Fortran compiler
- `make`

1) Build and install `cdf` library

The CDF build script download sources and builds the [NASA CDF library](https://cdf.gsfc.nasa.gov/)

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

2) Build and install `qdipole` library

```
$ cd qdipole/
$ ./configure
$ make build
$ sudo make install
```

By default the ``qdipole`` library is installed in the ``/usr`` system directory.
To install the library in a custom location use configuration with a custom prefix::

```
./configure --prefix=<install prefix>
```

3) Build and install `eoxmagmod` pip dependencies

```
pip3 install 'numpy<2'
pip3 install scipy
pip3 install 'SpacePy>0.5'
```

4) Finally, install `eoxmagmod` package

```
pip3 install ./eoxmagmod
```

### Testing

To test the fresh installation, leave the `MagneticModel` directory and run
the following command:

```
python3 -m unittest discover -p '[a-z]*.py' -v eoxmagmod
```
