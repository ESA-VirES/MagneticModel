#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Expansion with 2D Fourier series coefficients
#  - Geomagnetic Model  - tests
#
# Author: Martin Paces <martin.paces@eox.at>
#-------------------------------------------------------------------------------
# Copyright (C) 2022 EOX IT Services GmbH
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
from collections import namedtuple
from numpy import pi, stack, empty, zeros, prod, asarray
from numpy.random import uniform
from numpy.testing import assert_allclose
from numpy.lib.stride_tricks import as_strided
from eoxmagmod._pymm import (
    POTENTIAL, GRADIENT, POTENTIAL_AND_GRADIENT,
    GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
    INTERP_C0, INTERP_C1,
    convert, fourier2d, sheval, sheval2dfs,
)
from eoxmagmod.tests.data import mio


Fourier2DCoeffSet = namedtuple(
    "Fourier2DCoeffSet", [
        "coef_ab", "coef_nm", "min_degree1", "min_degree2", "scale1", "scale2",
    ]
)


class SphericalHarmonicsWithFourier2DCoeffMixIn:
    options = {}
    scale_potential = 1.0
    scale_gradient = [1.0, 1.0, 1.0]
    source_coordinate_system = None
    target_coordinate_system = None
    is_internal = True
    coef_set = None
    time1_range = (0.0, 1.0)
    time2_range = (0.0, 24.0)

    def times(self, shape=1):
        return (
            uniform(*self.time1_range, shape),
            uniform(*self.time2_range, shape),
        )

    def coordinates(self, shape=1):
        raise NotImplementedError

    @classmethod
    def eval_sheval2dfs(cls, times1, times2, coords, mode):
        return sheval2dfs(
            times1, times2, coords, cls.coef_set,
            coord_type_in=cls.source_coordinate_system,
            coord_type_out=cls.target_coordinate_system,
            mode=mode,
            is_internal=cls.is_internal,
            **cls.options
        )

    @classmethod
    def eval_reference_sheval2dfs(cls, times1, times2, coords):

        # make sure the times are arrays
        times1 = asarray(times1)
        times2 = asarray(times2)
        assert times1.shape == times2.shape

        # allocate multi-time coefficient arrays
        degree = cls.coef_set.coef_nm[:, 0].max()
        coeff_size = (degree+2)*(degree+1)//2
        coeff = zeros((*times1.shape, coeff_size, 2))

        # evaluate coefficients from the 2D Fourier series

        def _eval_fourier2d_coeff(times1, times2, coef_set):

            # mapping raw coefficients to G/H arrays
            degree = coef_set.coef_nm[:, 0]
            order = coef_set.coef_nm[:, 1]
            idx = ((degree+1)*degree)//2 + abs(order)

            coef_g_sel = (order >= 0).nonzero()[0]
            coef_h_sel = (order < 0).nonzero()[0]
            coef_g_idx = idx[coef_g_sel]
            coef_h_idx = idx[coef_h_sel]

            coeff_sparse = fourier2d(
                times1, times2,
                coef_set.coef_ab,
                coef_set.min_degree1,
                coef_set.min_degree2,
                coef_set.scale1,
                coef_set.scale2
            )

            coeff[..., coef_g_idx, 0] = coeff_sparse[..., coef_g_sel]
            coeff[..., coef_h_idx, 1] = coeff_sparse[..., coef_h_sel]

        _eval_fourier2d_coeff(times1, times2, cls.coef_set)

        # evaluate the spherical harmonics
        return sheval(
            coords,
            mode=POTENTIAL_AND_GRADIENT,
            is_internal=cls.is_internal,
            degree=degree,
            coef=coeff,
            coord_type_in=cls.source_coordinate_system,
            coord_type_out=cls.target_coordinate_system,
            **cls.options
        )

    def _test_sheval2dfs_potential_and_gradient(self, times_shape, coords_shape):
        times1, times2 = self.times(times_shape)
        coords = self.coordinates(coords_shape)
        potential, gradient = self.eval_sheval2dfs(times1, times2, coords, POTENTIAL_AND_GRADIENT)
        potential_ref, gradient_ref = self.eval_reference_sheval2dfs(times1, times2, coords)
        assert_allclose(potential, potential_ref, atol=1e-6)
        assert_allclose(gradient, gradient_ref, atol=1e-6)

    def _test_sheval2dfs_potential(self, times_shape, coords_shape):
        times1, times2 = self.times(times_shape)
        coords = self.coordinates(coords_shape)
        potential = self.eval_sheval2dfs(times1, times2, coords, POTENTIAL)
        potential_ref, _ = self.eval_reference_sheval2dfs(times1, times2, coords)
        assert_allclose(potential, potential_ref, atol=1e-6)

    def _test_sheval2dfs_gradient(self, times_shape, coords_shape):
        times1, times2 = self.times(times_shape)
        coords = self.coordinates(coords_shape)
        gradient = self.eval_sheval2dfs(times1, times2, coords, GRADIENT)
        _, gradient_ref = self.eval_reference_sheval2dfs(times1, times2, coords)
        assert_allclose(gradient, gradient_ref, atol=1e-6)

    def test_sheval2dfs_potential_and_gradient_T0X0(self):
        self._test_sheval2dfs_potential_and_gradient((), ())

    def test_sheval2dfs_potential_T0X0(self):
        self._test_sheval2dfs_potential((), ())

    def test_sheval2dfs_gradient_T0X0(self):
        self._test_sheval2dfs_gradient((), ())

    def test_sheval2dfs_potential_and_gradient_T0X1(self):
        self._test_sheval2dfs_potential_and_gradient((), (10,))

    def test_sheval2dfs_potential_T0X1(self):
        self._test_sheval2dfs_potential((), (10,))

    def test_sheval2dfs_gradient_T0X1(self):
        self._test_sheval2dfs_gradient((), (10,))

    def test_sheval2dfs_potential_and_gradient_T0X2(self):
        self._test_sheval2dfs_potential_and_gradient((), (10, 5))

    def test_sheval2dfs_potential_T0X2(self):
        self._test_sheval2dfs_potential((), (10, 5))

    def test_sheval2dfs_gradient_T0X2(self):
        self._test_sheval2dfs_gradient((), (10, 5))

    def test_sheval2dfs_potential_and_gradient_T1X0(self):
        self._test_sheval2dfs_potential_and_gradient((10,), ())

    def test_sheval2dfs_potential_T1X0(self):
        self._test_sheval2dfs_potential((10,), ())

    def test_sheval2dfs_gradient_T1X0(self):
        self._test_sheval2dfs_gradient((10,), ())

    def test_sheval2dfs_potential_and_gradient_T1X1(self):
        self._test_sheval2dfs_potential_and_gradient((10,), (10,))

    def test_sheval2dfs_potential_T1X1(self):
        self._test_sheval2dfs_potential((10,), (10,))

    def test_sheval2dfs_gradient_T1X1(self):
        self._test_sheval2dfs_gradient((10,), (10,))

    def test_sheval2dfs_potential_and_gradient_T1X2(self):
        self._test_sheval2dfs_potential_and_gradient((10,), (10, 5))

    def test_sheval2dfs_potential_T1X2(self):
        self._test_sheval2dfs_potential((10,), (10, 5))

    def test_sheval2dfs_gradient_T1X2(self):
        self._test_sheval2dfs_gradient((10,), (10, 5))

    def test_sheval2dfs_potential_and_gradient_T2X1(self):
        self._test_sheval2dfs_potential_and_gradient((10, 5), (10,))

    def test_sheval2dfs_potential_T2X1(self):
        self._test_sheval2dfs_potential((10, 5), (10,))

    def test_sheval2dfs_gradient_T2X1(self):
        self._test_sheval2dfs_gradient((10, 5), (10,))

#-------------------------------------------------------------------------------
# type of SH expansion

class SHTypeInternalMIO:
    is_internal = True
    coef_set = Fourier2DCoeffSet(
        mio.COEFF_I, mio.NMMAP, mio.SMIN, mio.PMIN, 2*pi, pi/12
    )


class SHTypeExternalMIO:
    is_internal = False
    coef_set = Fourier2DCoeffSet(
        mio.COEFF_E, mio.NMMAP, mio.SMIN, mio.PMIN, 2*pi, pi/12
    )


#-------------------------------------------------------------------------------
# sources coordinate systems

class SourceSpherical:
    source_coordinate_system = GEOCENTRIC_SPHERICAL

    def coordinates(self, shape=1):
        return stack((
            uniform(-90, 90, shape),
            uniform(-180, 180, shape),
            uniform(6371.2, 6371.2*1.5, shape),
        ), axis=-1)


class SourceGeodetic:
    source_coordinate_system = GEODETIC_ABOVE_WGS84

    def coordinates(self, shape=1):
        return stack((
            uniform(-90, 90, shape),
            uniform(-180, 180, shape),
            uniform(-50, 150, shape),
        ), axis=-1)


class SourceCartesian:
    source_coordinate_system = GEOCENTRIC_CARTESIAN

    def coordinates(self, shape=1):
        return convert(
            stack((
                uniform(-90, 90, shape),
                uniform(-180, 180, shape),
                uniform(6371.2, 6371.2*1.5, shape),
            ), axis=-1),
            GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN
        )

#-------------------------------------------------------------------------------

class TestSHEval2DFourierCartesian2CartesianInternalMIO(TestCase, SourceCartesian, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEval2DFourierCartesian2CartesianExternalMIO(TestCase, SourceCartesian, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEval2DFourierCartesian2SphericalInternalMIO(TestCase, SourceCartesian, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEval2DFourierCartesian2SphericalExternalMIO(TestCase, SourceCartesian, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEval2DFourierCartesian2WGS84InternalMIO(TestCase, SourceCartesian, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEval2DFourierCartesian2WGS84ExternalMIO(TestCase, SourceCartesian, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEval2DFourierSpherical2CartesianInternalMIO(TestCase, SourceSpherical, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEval2DFourierSpherical2CartesianExternalMIO(TestCase, SourceSpherical, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEval2DFourierSpherical2SphericalInternalMIO(TestCase, SourceSpherical, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEval2DFourierSpherical2SphericalExternalMIO(TestCase, SourceSpherical, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEval2DFourierSpherical2WGS84InternalMIO(TestCase, SourceSpherical, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEval2DFourierSpherical2WGS84ExternalMIO(TestCase, SourceSpherical, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEval2DFourierWGS842CartesianInternalMIO(TestCase, SourceGeodetic, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEval2DFourierWGS842CartesianExternalMIO(TestCase, SourceGeodetic, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEval2DFourierWGS842SphericalInternalMIO(TestCase, SourceGeodetic, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEval2DFourierWGS842SphericalExternalMIO(TestCase, SourceGeodetic, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEval2DFourierWGS842WGS84InternalMIO(TestCase, SourceGeodetic, SHTypeInternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEval2DFourierWGS842WGS84ExternalMIO(TestCase, SourceGeodetic, SHTypeExternalMIO, SphericalHarmonicsWithFourier2DCoeffMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEval2DFourierCart2CartScaled(TestSHEval2DFourierCartesian2CartesianInternalMIO):
    options = {"scale_potential": 2.0, "scale_gradient": [0.5, 1.0, -1.0]}
    scale_potential = 2.0
    scale_gradient = [0.5, 1.0, -1.0]

class TestSHEval2DFourierSph2SphScaled(TestSHEval2DFourierSpherical2SphericalInternalMIO):
    options = {"scale_gradient": -1.0}
    scale_gradient = [-1.0, -1.0, -1.0]

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
