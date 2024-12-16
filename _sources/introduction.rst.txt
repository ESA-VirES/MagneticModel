Introduction
============

This is the documentation of the ``eoxmagmod`` utility Python package used by
the *VirES for Swarm* (https://vires.services) service for the calculation
of:

- geomagnetic models (Spherical Harmonic expansion)
- magnetic coordinates (quasi-dipole and dipole latitude and longitude, magnetic local time)
- coordinate system conversions (geodetic, geocentric spherical, geocentric Cartesian coordinates)
- time conversions (modified Julian date 2000, decimal year)
- and other auxiliary calculations.

This package has been developed and maintained as a collection of various
performance-optimized support algorithms serving the VirES for Swarm service,
usability and portability being our prime objective.

Nonetheless, this packages if freely available to anybody under the terms of
the :doc:`MIT license<license>`. This documentation then provides information
about the installation and use of the module.

The ``eoxmagmod`` is pre-installed in the *Swarm Virtual Research Environment*
(https://vre.vires.services), a JupyterHub instance dedicated to exploration
of data provided by the VirES for Swarm service.
See an `example Jupyter notebook`_.

The ``eoxmagmod`` package may include 3rd party datasets. The package shall not
be used as the primary source of these data. These datasets are added for
for testing and may not be up-to-date.
The package itself currently does not provide any means for the synchronization
of these datasets.

To report an issue create a ticket in the `Git repository`_.

.. _example Jupyter notebook: https://notebooks.vires.services/notebooks/04b1_geomag-models-eoxmagmod
.. _Git repository: https://github.com/ESA-VirES/MagneticModel
