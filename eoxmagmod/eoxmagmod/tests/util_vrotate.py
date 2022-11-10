#-------------------------------------------------------------------------------
#
#  Testing vector rotation.
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
# pylint: disable=missing-docstring

from unittest import TestCase, main
from numpy import zeros, meshgrid, linspace
from numpy.random import random
from numpy.testing import assert_allclose
from eoxmagmod.util import vrotate
from eoxmagmod import (
    GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
    convert, vrot_sph2geod, vrot_sph2cart, vrot_cart2sph,
)


class VectorRotationMixIn:
    source_coordinate_system = None
    target_coordinate_system = None
    shape = (19, 19)

    @property
    def vectors(self):
        return 2.0*random(self.shape + (3,)) - 1.0

    @property
    def coords(self):
        coords = zeros(self.shape + (3,))
        coords[..., 1], coords[..., 0] = meshgrid(
            linspace(-180, 180, self.shape[1]),
            linspace(-90, 90, self.shape[0])
        )
        return coords

    def test_vrotate(self):
        vectors = self.vectors
        coords = self.coords
        assert_allclose(
            self.eval_vrotate(vectors, coords),
            self.reference_vrotate(vectors, coords)
        )

    @classmethod
    def eval_vrotate(cls, vectors, coords):
        coord_in, coord_out = cls.get_input_and_output_coordinates(coords)
        return vrotate(
            vectors, coord_in=coord_in, coord_out=coord_out,
            coord_type_in=cls.source_coordinate_system,
            coord_type_out=cls.target_coordinate_system,
        )

    @classmethod
    def get_input_and_output_coordinates(cls, coords):
        raise NotImplementedError

    @classmethod
    def reference_vrotate(cls, vectors, coords):
        raise NotImplementedError


class VectorRotationIdentityMixIn(VectorRotationMixIn):

    @classmethod
    def get_input_and_output_coordinates(cls, coords):
        return None, None

    @classmethod
    def reference_vrotate(cls, vectors, coords):
        return vectors


class VectorRotationSphericalToCartesianMixIn(VectorRotationMixIn):

    @classmethod
    def get_input_and_output_coordinates(cls, coords):
        return coords, None

    @classmethod
    def reference_vrotate(cls, vectors, coords):
        return vrot_sph2cart(vectors, coords[..., 0], coords[..., 1])


class VectorRotationCartesianToSphericalMixIn(VectorRotationMixIn):

    @classmethod
    def get_input_and_output_coordinates(cls, coords):
        return None, coords

    @classmethod
    def reference_vrotate(cls, vectors, coords):
        return vrot_cart2sph(vectors, coords[..., 0], coords[..., 1])


class VectorRotationSphericalToGeodeticMixIn(VectorRotationMixIn):

    @classmethod
    def get_input_and_output_coordinates(cls, coords):
        return coords, convert(
            coords, GEOCENTRIC_SPHERICAL, cls.target_coordinate_system
        )

    @classmethod
    def reference_vrotate(cls, vectors, coords):
        coords_in, coords_out = cls.get_input_and_output_coordinates(coords)
        return vrot_sph2geod(vectors, coords_out[..., 0] - coords_in[..., 0])


class VectorRotationGeodeticToSphericalMixIn(VectorRotationSphericalToGeodeticMixIn):

    @classmethod
    def get_input_and_output_coordinates(cls, coords):
        return convert(
            coords, cls.source_coordinate_system, GEOCENTRIC_SPHERICAL,
        ), coords

#-------------------------------------------------------------------------------

class TestVectorRotateWGS84ToWGS84(TestCase, VectorRotationIdentityMixIn):
    source_coordinate_system = GEODETIC_ABOVE_WGS84
    target_coordinate_system = GEODETIC_ABOVE_WGS84


class TestVectorRotateWGS84ToSpherical(TestCase, VectorRotationGeodeticToSphericalMixIn):
    source_coordinate_system = GEODETIC_ABOVE_WGS84
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestVectorRotateWGS84ToCartesian(TestCase, VectorRotationSphericalToCartesianMixIn):
    source_coordinate_system = GEODETIC_ABOVE_WGS84
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestVectorRotateSphericalToWGS84(TestCase, VectorRotationGeodeticToSphericalMixIn):
    source_coordinate_system = GEOCENTRIC_SPHERICAL
    target_coordinate_system = GEODETIC_ABOVE_WGS84


class TestVectorRotateSphericalToSpherical(TestCase, VectorRotationIdentityMixIn):
    source_coordinate_system = GEOCENTRIC_SPHERICAL
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestVectorRotateSphericalToCartesian(TestCase, VectorRotationSphericalToCartesianMixIn):
    source_coordinate_system = GEOCENTRIC_SPHERICAL
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestVectorRotateCartesianToWGS84(TestCase, VectorRotationCartesianToSphericalMixIn):
    source_coordinate_system = GEOCENTRIC_CARTESIAN
    target_coordinate_system = GEODETIC_ABOVE_WGS84


class TestVectorRotateCartesianToSpherical(TestCase, VectorRotationCartesianToSphericalMixIn):
    source_coordinate_system = GEOCENTRIC_CARTESIAN
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestVectorRotateCartesianToCartesian(TestCase, VectorRotationIdentityMixIn):
    source_coordinate_system = GEOCENTRIC_CARTESIAN
    target_coordinate_system = GEOCENTRIC_CARTESIAN

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
