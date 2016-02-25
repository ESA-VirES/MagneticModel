## Magnetic Model

This repository contains various utilities related to Earth magnetic field
modelling and spherical harmonics.

The repository contains following directories:

- `eoxmagmod` - Collection models of the Earth magnetic field - python module
- `wmm` - World Magnetic Model 2015 - re-packed to be compiled as a shared
  library (dependency of the `eoxmagmod` package.) 

### Installation from Sources

#### WMM

```
$ cd wmm/ 
$ ./configure --prefix=/usr/ --libdir='${exec_prefix}/lib64'
$ make build
$ sudo make install
```

#### EOxMagMod

```
$ cd eoxmagmod/
$ cd eoxmagmod/qdipolelib
$ make clean build install
$ cd -
$ python ./setup.py build
$ sudo python ./setup.py install
```
