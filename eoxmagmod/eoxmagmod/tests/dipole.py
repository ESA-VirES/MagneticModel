#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Expansion - Geomagnetic Model - tests
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2018 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring,no-name-in-module,too-few-public-methods,line-too-long

from unittest import TestCase, main
from itertools import product
from random import random
from numpy import asarray, stack
from numpy.testing import assert_allclose
from eoxmagmod._pymm import (
    POTENTIAL, GRADIENT, POTENTIAL_AND_GRADIENT,
    GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
    convert, relradpow, loncossin, legendre, spharpot, sphargrd,
)
from eoxmagmod.dipole_coords import convert_to_dipole, vrot_from_dipole
from eoxmagmod.dipole import sheval_dipole
from eoxmagmod.tests.data import mma_external, mma_internal


class DipoleSphericalHarmonicsMixIn:
    options = {}
    scale_potential = 1.0
    scale_gradient = [1.0, 1.0, 1.0]
    source_coordinate_system = None
    target_coordinate_system = None
    is_internal = True
    degree = None
    coeff = None
    lat_ngp = None
    lon_ngp = None

    @classmethod
    def eval_sheval(cls, coords, mode):
        return sheval_dipole(
            coords, mode=mode, is_internal=cls.is_internal,
            degree=cls.degree, coef=cls.coeff,
            lat_ngp=cls.lat_ngp, lon_ngp=cls.lon_ngp,
            coord_type_in=cls.source_coordinate_system,
            coord_type_out=cls.target_coordinate_system,
            **cls.options
        )

    @classmethod
    def reference_sheval(cls, coords):

        coords_dipole = convert_to_dipole(
            coords, cls.lat_ngp, cls.lon_ngp, cls.source_coordinate_system
        )

        potential, gradient = cls._spherical_harmonics(
            coords_dipole[..., 0],
            coords_dipole[..., 1],
            coords_dipole[..., 2],
        )

        gradient = cls._rotate_gradient(gradient, coords_dipole, coords)

        potential *= cls.scale_potential
        gradient *= cls.scale_gradient

        return potential, gradient

    @classmethod
    def get_series(cls, degree, latitude, longitude, radius):
        rad_series = relradpow(radius, degree, is_internal=cls.is_internal)
        cos_sin_series = loncossin(longitude, degree)
        p_series, dp_series = legendre(latitude, degree)
        return rad_series, cos_sin_series, p_series, dp_series

    @classmethod
    def _spherical_harmonics(cls, latitude, longitude, radius):
        rad_series, cos_sin_series, p_series, dp_series = cls.get_series(
            cls.degree, latitude, longitude, radius
        )
        potential = spharpot(
            radius, cls.coeff, p_series, rad_series,
            cos_sin_series, degree=cls.degree
        )
        gradient = sphargrd(
            latitude, cls.coeff, p_series, dp_series, rad_series,
            cos_sin_series, is_internal=cls.is_internal, degree=cls.degree
        )
        return potential, gradient

    @classmethod
    def _rotate_gradient(cls, vectors, coords_dipole, coords_in):
        lat_dipole, lon_dipole = coords_dipole[..., 0], coords_dipole[..., 1]
        if cls.target_coordinate_system == GEOCENTRIC_CARTESIAN:
            lat_out, lon_out = None, None
        else:
            coords_out = convert(
                coords_in, cls.source_coordinate_system,
                cls.target_coordinate_system
            )
            lat_out, lon_out = coords_out[..., 0], coords_out[..., 1]

        return vrot_from_dipole(
            vectors, cls.lat_ngp, cls.lon_ngp, lat_dipole, lon_dipole,
            lat_out, lon_out, cls.target_coordinate_system
        )

    def test_sheval_potential_and_gradient(self):
        coords = self.coordinates
        potential_ref, gradient_ref = self.reference_sheval(coords)
        potential, gradient = self.eval_sheval(coords, POTENTIAL_AND_GRADIENT)
        assert_allclose(potential, potential_ref, atol=1e-6)
        assert_allclose(gradient, gradient_ref, atol=1e-6)

    def test_sheval_potential(self):
        coords = self.coordinates
        potential_ref, _ = self.reference_sheval(coords)
        potential = self.eval_sheval(coords, POTENTIAL)
        assert_allclose(potential, potential_ref, atol=1e-6)

    def test_sheval_gradient(self):
        coords = self.coordinates
        _, gradient_ref = self.reference_sheval(coords)
        gradient = self.eval_sheval(coords, GRADIENT)
        assert_allclose(gradient, gradient_ref, atol=1e-6)

#-------------------------------------------------------------------------------
# sources

class SourceSpherical:
    source_coordinate_system = GEOCENTRIC_SPHERICAL

    @property
    def coordinates(self):
        return asarray([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ])


class SourceGeodetic:
    source_coordinate_system = GEODETIC_ABOVE_WGS84

    @property
    def coordinates(self):
        return asarray([
            (lat, lon, -50. + 200.0 * random()) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ])


class SourceCartesian:
    source_coordinate_system = GEOCENTRIC_CARTESIAN

    @property
    def coordinates(self):
        return convert(asarray([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]), GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN)

#-------------------------------------------------------------------------------

class SHTypeInternal:
    lat_ngp = 80.08
    lon_ngp = -72.22
    is_internal = True
    degree = mma_internal.DEGREE
    coeff = stack((mma_internal.COEF_G, mma_internal.COEF_H), axis=-1)


class SHTypeExternal:
    lat_ngp = 80.08
    lon_ngp = -72.22
    is_internal = False
    degree = mma_external.DEGREE
    coeff = stack((mma_external.COEF_Q, mma_external.COEF_S), axis=-1)

#-------------------------------------------------------------------------------

class TestDipoleSHEvalCartesian2CartesianInternal(TestCase, SourceCartesian, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestDipoleSHEvalCartesian2CartesianExternal(TestCase, SourceCartesian, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestDipoleSHEvalCartesian2SphericalInternal(TestCase, SourceCartesian, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestDipoleSHEvalCartesian2SphericalExternal(TestCase, SourceCartesian, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestDipoleSHEvalCartesian2WGS84Internal(TestCase, SourceCartesian, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestDipoleSHEvalCartesian2WGS84External(TestCase, SourceCartesian, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestDipoleSHEvalSpherical2CartesianInternal(TestCase, SourceSpherical, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestDipoleSHEvalSpherical2CartesianExternal(TestCase, SourceSpherical, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestDipoleSHEvalSpherical2SphericalInternal(TestCase, SourceSpherical, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestDipoleSHEvalSpherical2SphericalExternal(TestCase, SourceSpherical, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestDipoleSHEvalSpherical2WGS84Internal(TestCase, SourceSpherical, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestDipoleSHEvalSpherical2WGS84External(TestCase, SourceSpherical, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestDipoleSHEvalWGS842CartesianInternal(TestCase, SourceGeodetic, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestDipoleSHEvalWGS842CartesianExternal(TestCase, SourceGeodetic, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestDipoleSHEvalWGS842SphericalInternal(TestCase, SourceGeodetic, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestDipoleSHEvalWGS842SphericalExternal(TestCase, SourceGeodetic, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestDipoleSHEvalWGS842WGS84Internal(TestCase, SourceGeodetic, SHTypeInternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestDipoleSHEvalWGS842WGS84External(TestCase, SourceGeodetic, SHTypeExternal, DipoleSphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestDipoleSHEvalCart2CartScaled(TestDipoleSHEvalCartesian2CartesianInternal):
    options = {"scale_potential": 2.0, "scale_gradient": [0.5, 1.0, -1.0]}
    scale_potential = 2.0
    scale_gradient = [0.5, 1.0, -1.0]

class TestDipoleSHEvalSph2SphScaled(TestDipoleSHEvalSpherical2SphericalInternal):
    options = {"scale_gradient": -1.0}
    scale_gradient = [-1.0, -1.0, -1.0]

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
