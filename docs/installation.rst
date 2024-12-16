Installation
============

The ``eoxmagmod`` packages is installed from sources. The sources can be
obtained either from the `releases`_ (``eoxmagmod-<version>.tar.gz``) or by
cloning the `git repository`_ (latest development version).

The sources contain these directories:

- ``eoxmagmod`` -  containing the sources of the ``eoxmagmod`` Python package
- ``qdipole`` - containing sources of the shared library calculating
  Quasi-Dipole apex coordinates (QD-latitude, QD-longitude, and Magnetic Local
  Time). This shared library is a dependency of the ``eoxmagmod`` package.
- ``libcdf`` - containing installator of the `NASA CDF library`_.
  The installator script downloads its sources, builds and installs the shared
  library.  This shared library is a dependency of the ``eoxmagmod`` package.

The installation requires C and Fortran compilers and Make build tool.

The package is tested to work on x86_64 GNU/Linux using GCC and GFortran
compilers. Porting to other platforms should be technically possible, although
formally not supported.

.. _releases: https://github.com/ESA-VirES/MagneticModel/releases
.. _git repository: https://github.com/ESA-VirES/MagneticModel
.. _NASA CDF library: https://cdf.gsfc.nasa.gov/


Conda Installation
------------------

`Conda` is package management system used for distribution of Python packages
and their dependencies. Conda is part of the `Anaconda`_ and `Miniconda`_
installers.

Due to the automatic handling of dependencies the Conda installation is simpler
than the generic source installation.

.. _Anaconda: https://docs.anaconda.com/anaconda/
.. _Miniconda: https://docs.anaconda.com/miniconda/

To install  ``eoxmagmod`` in a Conda environment follow these steps:

Before the start, change to the root source directory ``MagneticModel/``.

1) Build the shared libraries::

    conda install conda-build  # install build tools
    conda build ./qdipole      # build qdipole local binary installation package
    conda build ./libcdf       # build cdf local binary installation package
    conda build purge          # remove build tools and artefacts

The built packages can be listed with::

    conda search --use-local --offline


The output may look like::

    Loading channels: done
    # Name                       Version           Build  Channel
    cdf                            3.9.0      h3218e01_0  .../conda/conda-bld
    qdipole                        0.6.0      h3218e01_0  .../conda/conda-bld


2) Create a new Conda environment with the ``eoxmagmod`` dependencies:

The sources contain Conda environment file ``conda_env.yaml``  defining all
required dependencies. Its content may look like:

.. literalinclude:: ../conda_env.yaml
  :language: YAML


Create a new Conda environment from the provided ``conda_env.yaml`` file::

    conda env create -n <environment-name> -f conda_env.yaml


3) Activate the new Conda environment and install the ``eoxmagmod`` package::

    conda activate <environment-name>
    pip install ./eoxmagmod


Generic Source Installation
---------------------------

The installation requires the C (GCC) and Fortran (GFortran) compilers and
the `make` build command.

Before the start, change to the root source directory ``MagneticModel/``.

1) Build and install CDF library::

    cd libcdf/
    make build
    make test
    make install
    make clean

By default, the library is installed in ``/usr/local/cdf`` directory.
To install it to a different path override the ``INSTALLDIR`` variable::

    make install INSTALLDIR=<install directory>

2) Build and install ``qdipole`` library::

    cd qdipole/
    ./configure
    make build
    make install
    make clean
    make test

By default the ``qdipole`` library is installed in the ``/usr`` system directory.
To install the library in a custom location use configuration with a custom prefix::

   ./configure --prefix=<install prefix>

3) Install other ``eoxmagmod`` dependencies (assuming pure ``pip`` installation)::

    pip3 install 'numpy<2'
    pip3 install scipy
    pip3 install 'SpacePy>0.5'

3) Finally install ``eoxmagmod`` package::

    pip3 install ./eoxmagmod


Test Installation
-----------------
To test the installation, leave the ``MagneticModel`` directory and run the
test command::

    cd ..      # make sure you leave the MagneticModel directory!
    python3 -m unittest discover -p '[a-z]*.py' -v eoxmagmod
