#-------------------------------------------------------------------------------
#
#  Magnetic Dipole Coordinates - tests
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
# pylint: disable=missing-docstring, invalid-name

from unittest import TestCase, main
from itertools import product
from math import pi
from numpy import array, asarray, zeros, linspace, meshgrid, sin, cos, dot
from numpy.random import random
from numpy.testing import assert_allclose
from eoxmagmod._pywmm import (
    convert, vrot_sph2cart, vrot_cart2sph,
    GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
)
from eoxmagmod.dipole_coords import (
    get_dipole_rotation_matrix,
    convert_to_dipole,
    vrot_dipole2spherical,
)

DEG2RAD = pi/180.0


class TestDipoleRotationMatrix(TestCase):

    @staticmethod
    def reference_rotation_matrix(latitude, longitude):
        sin_lat, cos_lat = sin(DEG2RAD*latitude), cos(DEG2RAD*latitude)
        sin_lon, cos_lon = sin(DEG2RAD*longitude), cos(DEG2RAD*longitude)
        matrix = dot(
            # rotate around azimuth axis by -longitude
            array([
                [cos_lon, -sin_lon, 0],
                [sin_lon, cos_lon, 0],
                [0, 0, 1],
            ]),
            # rotate around elevation axis by 90dg - latitude
            array([
                [sin_lat, 0, cos_lat],
                [0, 1, 0],
                [-cos_lat, 0, sin_lat],
            ])
        )
        return matrix

    @staticmethod
    def eval_rotation_matrix(latitude, longitude):
        return get_dipole_rotation_matrix(latitude, longitude)

    def test_rotation_matrix(self):
        coords = [
            (lat, lon) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]

        for lat, lon in coords:
            matrix = self.eval_rotation_matrix(lat, lon)
            assert_allclose(
                matrix,
                self.reference_rotation_matrix(lat, lon),
                atol=1e-14
            )
            assert_allclose(
                dot(matrix.transpose(), matrix),
                [(1, 0, 0), (0, 1, 0), (0, 0, 1)],
                atol=1e-14
            )


class TestConvertToDipoleCoordinates(TestCase):

    @staticmethod
    def reference_convert_to_dipole(coords, latitude, longitude):
        rotation_matrix = get_dipole_rotation_matrix(latitude, longitude)
        coords = convert(coords, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN)
        coords = dot(coords, rotation_matrix)
        return coords

    @staticmethod
    def eval_convert_to_dipole(coords, latitude, longitude):
        # to avoid pole longitude ambiguity compare Cartesian coordinates
        return convert(
            convert_to_dipole(coords, latitude, longitude),
            GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN
        )

    @property
    def coordinates(self):
        return array([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ])


    def test_convert_to_dipole(self):
        north_pole_coords = [
            (lat, lon) for lat, lon
            in product(range(-90, 91, 10), range(-180, 181, 20))
        ]
        for lat, lon in north_pole_coords:
            coords = self.coordinates
            assert_allclose(
                self.eval_convert_to_dipole(coords, lat, lon),
                self.reference_convert_to_dipole(coords, lat, lon),
                atol=1e-8
            )

    def test_convert_to_dipole_sanity_check(self):
        assert_allclose(
            self.eval_convert_to_dipole([
                (80, -170, 1.0),
                (-80, 10, 1.0),
                (-10, -170, 1.0),
                (10, 10, 1.0),
                (0, -80, 1.0),
                (0, 100, 1.0),
            ], 80, -170),
            [
                (0, 0, 1),
                (0, 0, -1),
                (1, 0, 0),
                (-1, 0, 0),
                (0, 1, 0),
                (0, -1, 0),
            ],
            atol=1e-12
        )


class TestVRotDipoleToSpherical(TestCase):
    shape = (37, 37)

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
        coords[..., 2] = 1.0
        return coords

    @staticmethod
    def reference_vrot_dipole2spherical(vectors, coords, latitude, longitude):
        coords = asarray(coords)
        rotation_matrix = get_dipole_rotation_matrix(latitude, longitude)
        coords_dipole = convert_to_dipole(coords, latitude, longitude)
        lat_spherical = coords[..., 0]
        lon_spherical = coords[..., 1]
        lat_dipole = coords_dipole[..., 0]
        lon_dipole = coords_dipole[..., 1]
        vectors = vrot_sph2cart(vectors, lat_dipole, lon_dipole)
        vectors = dot(vectors, rotation_matrix.transpose())
        vectors = vrot_cart2sph(vectors, lat_spherical, lon_spherical)
        return vectors

    @staticmethod
    def eval_vrot_dipole2spherical(vectors, coords, latitude, longitude):
        coords = asarray(coords)
        coords_dipole = convert_to_dipole(coords, latitude, longitude)
        lat_spherical = coords[..., 0]
        lon_spherical = coords[..., 1]
        lat_dipole = coords_dipole[..., 0]
        lon_dipole = coords_dipole[..., 1]
        return vrot_dipole2spherical(
            vectors, lat_dipole, lon_dipole, lat_spherical, lon_spherical,
            latitude, longitude
        )

    def test_vrot_dipole2spherical(self):
        north_pole_coords = [
            (lat, lon) for lat, lon
            in product(range(-90, 91, 10), range(-180, 181, 20))
        ]
        for lat, lon in north_pole_coords:
            coords = self.coords
            vects = self.vectors
            assert_allclose(
                self.eval_vrot_dipole2spherical(vects, coords, lat, lon),
                self.reference_vrot_dipole2spherical(vects, coords, lat, lon),
                atol=1e-12
            )

    def test_vrot_dipole2spherical_sanity_check(self):
        lat_ngp, lon_ngp = 80, -170
        vectors = array([
            (1, 0, 0),
            (0, 1, 0),
            (0, 0, 1),
        ])
        sin10, cos10 = sin(DEG2RAD*10), cos(DEG2RAD*10)
        input_output_pairs = [
            ((-10, -170), (0, 0), [(1, 0, 0), (0, 1, 0), (0, 0, 1),]),
            ((10, 10), (0, 180), [(1, 0, 0), (0, 1, 0), (0, 0, 1),]),
            (
                (0, -80), (0, 90),
                [(cos10, -sin10, 0), (sin10, cos10, 0), (0, 0, 1)]
            ),
            (
                (0, 100), (0, -90),
                [(cos10, sin10, 0), (-sin10, cos10, 0), (0, 0, 1)]
            ),
            ((80, -170), (90, 0), [(1, 0, 0), (0, 1, 0), (0, 0, 1),]),
            ((-80, 10), (-90, 0), [(-1, 0, 0), (0, -1, 0), (0, 0, 1),]),
        ]
        for (lat_sph, lon_sph), (lat_dip, lon_dip), expected in input_output_pairs:
            assert_allclose(
                vrot_dipole2spherical(
                    vectors, lat_dip, lon_dip, lat_sph, lon_sph, lat_ngp, lon_ngp
                ), expected, atol=1e-12
            )


if __name__ == "__main__":
    main()
