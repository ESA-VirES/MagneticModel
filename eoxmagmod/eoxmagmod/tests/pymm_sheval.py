#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Expansion - Geomagnetic Model - tests
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2018-2022 EOX IT Services GmbH
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
    convert, vrot_sph2geod, vrot_sph2cart,
    relradpow, loncossin, legendre,
    spharpot, sphargrd, sheval,
)
from eoxmagmod.tests.data import sifm, mma_external


class SphericalHarmonicsMixIn:
    options = {}
    scale_potential = 1.0
    scale_gradient = [1.0, 1.0, 1.0]
    source_coordinate_system = None
    target_coordinate_system = None
    is_internal = True
    degree = None
    coeff = None

    @classmethod
    def eval_sheval(cls, coords, mode):
        return sheval(
            coords, mode=mode, is_internal=cls.is_internal,
            degree=cls.degree, coef=cls.coeff,
            coord_type_in=cls.source_coordinate_system,
            coord_type_out=cls.target_coordinate_system,
            **cls.options
        )

    @classmethod
    def reference_sheval(cls, coords):

        coords_spherical = convert(
            coords, cls.source_coordinate_system, GEOCENTRIC_SPHERICAL
        )

        potential, gradient = cls._spherical_harmonics(
            coords_spherical[..., 0],
            coords_spherical[..., 1],
            coords_spherical[..., 2],
        )

        gradient = cls._rotate_gradient(gradient, coords_spherical)

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
    def _rotate_gradient(cls, vectors, coords):
        if cls.target_coordinate_system == GEOCENTRIC_SPHERICAL:
            return vectors
        if cls.target_coordinate_system == GEOCENTRIC_CARTESIAN:
            latd = coords[..., 0]
            lond = coords[..., 1]
            return vrot_sph2cart(vectors, latd, lond)
        if cls.target_coordinate_system == GEODETIC_ABOVE_WGS84:
            dlatd = convert(
                coords, GEOCENTRIC_SPHERICAL, cls.target_coordinate_system
            )[..., 0] - coords[..., 0]
            return vrot_sph2geod(vectors, dlatd)
        return None

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
    is_internal = True
    degree = sifm.DEGREE
    coeff = stack((sifm.COEF_G, sifm.COEF_H), axis=-1)


class SHTypeExternal:
    is_internal = False
    degree = mma_external.DEGREE
    coeff = stack((mma_external.COEF_Q, mma_external.COEF_S), axis=-1)

#-------------------------------------------------------------------------------

class TestSHEvalCartesian2CartesianInternal(TestCase, SourceCartesian, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalCartesian2CartesianExternal(TestCase, SourceCartesian, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalCartesian2SphericalInternal(TestCase, SourceCartesian, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalCartesian2SphericalExternal(TestCase, SourceCartesian, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalCartesian2WGS84Internal(TestCase, SourceCartesian, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalCartesian2WGS84External(TestCase, SourceCartesian, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalSpherical2CartesianInternal(TestCase, SourceSpherical, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalSpherical2CartesianExternal(TestCase, SourceSpherical, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalSpherical2SphericalInternal(TestCase, SourceSpherical, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalSpherical2SphericalExternal(TestCase, SourceSpherical, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalSpherical2WGS84Internal(TestCase, SourceSpherical, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalSpherical2WGS84External(TestCase, SourceSpherical, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalWGS842CartesianInternal(TestCase, SourceGeodetic, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalWGS842CartesianExternal(TestCase, SourceGeodetic, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalWGS842SphericalInternal(TestCase, SourceGeodetic, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalWGS842SphericalExternal(TestCase, SourceGeodetic, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalWGS842WGS84Internal(TestCase, SourceGeodetic, SHTypeInternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalWGS842WGS84External(TestCase, SourceGeodetic, SHTypeExternal, SphericalHarmonicsMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalCart2CartScaled(TestSHEvalCartesian2CartesianInternal):
    options = {"scale_potential": 2.0, "scale_gradient": [0.5, 1.0, -1.0]}
    scale_potential = 2.0
    scale_gradient = [0.5, 1.0, -1.0]

class TestSHEvalSph2SphScaled(TestSHEvalSpherical2SphericalInternal):
    options = {"scale_gradient": -1.0}
    scale_gradient = [-1.0, -1.0, -1.0]

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
