#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Expansion with temporal interpolation - Geomagnetic Model
#  - tests
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
from numpy import stack, empty, zeros, prod, asarray
from numpy.random import uniform
from numpy.testing import assert_allclose
from numpy.lib.stride_tricks import as_strided
from eoxmagmod._pymm import (
    POTENTIAL, GRADIENT, POTENTIAL_AND_GRADIENT,
    GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
    INTERP_C0, INTERP_C1,
    convert, interp, sheval, shevaltemp,
)
from eoxmagmod.tests.data import chaos_core, chaos_mma


# Coeficient set:
# coef_times - [Nt] time nodes of the SH coefficient time-series
# coef_coef  - [Nc,Nt] interpolated time-series of the SH coefficients
# coef_nm    - [Nc,2] degree/order mapping of the coefficients
# spline_order - currently only 1 and 2 are allowed

CoeffSet = namedtuple("CoeffSet", ["coef_times", "coef_coef", "coef_nm", "spline_order"])


class SphericalHarmonicsWithCoeffInterpolationMixIn:

    options = {}
    scale_potential = 1.0
    scale_gradient = [1.0, 1.0, 1.0]
    source_coordinate_system = None
    target_coordinate_system = None
    is_internal = True
    coef_sets = [] # list of coefficient sets

    def times(self, shape=1):

        start = max([cs.coef_times[0] for cs in self.coef_sets])
        stop = min([cs.coef_times[-1] for cs in self.coef_sets])

        return uniform(start, stop, shape)

    def coordinates(self, shape=1):
        raise NotImplementedError

    @classmethod
    def eval_shevaltemp(cls, times, coords, mode):
        return shevaltemp(
            times,
            coords,
            coef_set_list=cls.coef_sets,
            coord_type_in=cls.source_coordinate_system,
            coord_type_out=cls.target_coordinate_system,
            mode=mode,
            is_internal=cls.is_internal,
            **cls.options
        )

    @classmethod
    def eval_reference_shevaltemp(cls, times, coords):

        # make sure the time is an array
        times = asarray(times)

        # allocate multi-time coefficient arrays
        degree = max(coef_set.coef_nm[:, 0].max() for coef_set in cls.coef_sets)
        coeff_size = (degree+2)*(degree+1)//2
        coeff = zeros((*times.shape, coeff_size, 2))

        # interpolate coefficient sets

        def _interp_coef_set(times, coef_set):
            options = {
                1: {"kind":INTERP_C0,"extrapolate":True},
                2: {"kind":INTERP_C1},
            }[coef_set.spline_order]

            # mapping raw coefficients to G/H arrays
            degree = coef_set.coef_nm[:, 0]
            order = coef_set.coef_nm[:, 1]
            idx = ((degree+1)*degree)//2 + abs(order)

            coef_g_sel = (order >= 0).nonzero()[0]
            coef_h_sel = (order < 0).nonzero()[0]
            coef_g_idx = idx[coef_g_sel]
            coef_h_idx = idx[coef_h_sel]

            coeff_sparse = interp(
                times,
                coef_set.coef_times,
                coef_set.coef_coef,
                **options,
            )

            coeff[..., coef_g_idx, 0] = coeff_sparse[..., coef_g_sel]
            coeff[..., coef_h_idx, 1] = coeff_sparse[..., coef_h_sel]

        for coef_set in cls.coef_sets:
            _interp_coef_set(times, coef_set)

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

    def _test_shevaltemp_potential_and_gradient(self, times_shape, coords_shape):
        times = self.times(times_shape)
        coords = self.coordinates(coords_shape)
        potential, gradient = self.eval_shevaltemp(times, coords, POTENTIAL_AND_GRADIENT)
        potential_ref, gradient_ref = self.eval_reference_shevaltemp(times, coords)
        assert_allclose(potential, potential_ref, atol=1e-6)
        assert_allclose(gradient, gradient_ref, atol=1e-6)

    def _test_shevaltemp_potential(self, times_shape, coords_shape):
        times = self.times(times_shape)
        coords = self.coordinates(coords_shape)
        potential = self.eval_shevaltemp(times, coords, POTENTIAL)
        potential_ref, _ = self.eval_reference_shevaltemp(times, coords)
        assert_allclose(potential, potential_ref, atol=1e-6)

    def _test_shevaltemp_gradient(self, times_shape, coords_shape):
        times = self.times(times_shape)
        coords = self.coordinates(coords_shape)
        gradient = self.eval_shevaltemp(times, coords, GRADIENT)
        _, gradient_ref = self.eval_reference_shevaltemp(times, coords)
        assert_allclose(gradient, gradient_ref, atol=1e-6)

    def test_shevaltemp_potential_and_gradient_T0X0(self):
        self._test_shevaltemp_potential_and_gradient((), ())

    def test_shevaltemp_potential_T0X0(self):
        self._test_shevaltemp_potential((), ())

    def test_shevaltemp_gradient_T0X0(self):
        self._test_shevaltemp_gradient((), ())

    def test_shevaltemp_potential_and_gradient_T0X1(self):
        self._test_shevaltemp_potential_and_gradient((), (10,))

    def test_shevaltemp_potential_T0X1(self):
        self._test_shevaltemp_potential((), (10,))

    def test_shevaltemp_gradient_T0X1(self):
        self._test_shevaltemp_gradient((), (10,))

    def test_shevaltemp_potential_and_gradient_T0X2(self):
        self._test_shevaltemp_potential_and_gradient((), (10, 5))

    def test_shevaltemp_potential_T0X2(self):
        self._test_shevaltemp_potential((), (10, 5))

    def test_shevaltemp_gradient_T0X2(self):
        self._test_shevaltemp_gradient((), (10, 5))

    def test_shevaltemp_potential_and_gradient_T1X0(self):
        self._test_shevaltemp_potential_and_gradient((10,), ())

    def test_shevaltemp_potential_T1X0(self):
        self._test_shevaltemp_potential((10,), ())

    def test_shevaltemp_gradient_T1X0(self):
        self._test_shevaltemp_gradient((10,), ())

    def test_shevaltemp_potential_and_gradient_T1X1(self):
        self._test_shevaltemp_potential_and_gradient((10,), (10,))

    def test_shevaltemp_potential_T1X1(self):
        self._test_shevaltemp_potential((10,), (10,))

    def test_shevaltemp_gradient_T1X1(self):
        self._test_shevaltemp_gradient((10,), (10,))

    def test_shevaltemp_potential_and_gradient_T1X2(self):
        self._test_shevaltemp_potential_and_gradient((10,), (10, 5))

    def test_shevaltemp_potential_T1X2(self):
        self._test_shevaltemp_potential((10,), (10, 5))

    def test_shevaltemp_gradient_T1X2(self):
        self._test_shevaltemp_gradient((10,), (10, 5))

    def test_shevaltemp_potential_and_gradient_T2X1(self):
        self._test_shevaltemp_potential_and_gradient((10, 5), (10,))

    def test_shevaltemp_potential_T2X1(self):
        self._test_shevaltemp_potential((10, 5), (10,))

    def test_shevaltemp_gradient_T2X1(self):
        self._test_shevaltemp_gradient((10, 5), (10,))

#-------------------------------------------------------------------------------
# type of SH expansion

class SHTypeInternalStatic:
    is_internal = True
    coef_sets = [
        CoeffSet(chaos_core.TIMES[:1], chaos_core.COEFF[:,:1], chaos_core.NMMAP, 1),
    ]

class SHTypeInternalCore:
    is_internal = True
    coef_sets = [
        CoeffSet(chaos_core.TIMES, chaos_core.COEFF, chaos_core.NMMAP, 2),
    ]

class SHTypeInternalMMA:
    is_internal = True
    coef_sets = [
        CoeffSet(chaos_mma.TIMES_I_1, chaos_mma.COEFF_I_1, chaos_mma.NMMAP_I_1, 2),
        CoeffSet(chaos_mma.TIMES_I_2, chaos_mma.COEFF_I_2, chaos_mma.NMMAP_I_2, 2),
    ]

class SHTypeExternalMMA:
    is_internal = True
    coef_sets = [
        CoeffSet(chaos_mma.TIMES_E_1, chaos_mma.COEFF_E_1, chaos_mma.NMMAP_E_1, 2),
        CoeffSet(chaos_mma.TIMES_E_2, chaos_mma.COEFF_E_2, chaos_mma.NMMAP_E_2, 2),
    ]

#-------------------------------------------------------------------------------
# sources


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

class TestSHEvalTempCartesian2CartesianInternalCore(TestCase, SourceCartesian, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalTempCartesian2CartesianInternalMMA(TestCase, SourceCartesian, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalTempCartesian2CartesianExternalMMA(TestCase, SourceCartesian, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalTempCartesian2SphericalInternalCore(TestCase, SourceCartesian, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalTempCartesian2SphericalInternalMMA(TestCase, SourceCartesian, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalTempCartesian2SphericalExternalMMA(TestCase, SourceCartesian, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalTempCartesian2WGS84InternalCore(TestCase, SourceCartesian, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalTempCartesian2WGS84InternalMMA(TestCase, SourceCartesian, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalTempCartesian2WGS84ExternalMMA(TestCase, SourceCartesian, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalTempSpherical2CartesianInternalCore(TestCase, SourceSpherical, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalTempSpherical2CartesianInternalMMA(TestCase, SourceSpherical, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalTempSpherical2CartesianExternalMMA(TestCase, SourceSpherical, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalTempSpherical2SphericalInternalCore(TestCase, SourceSpherical, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalTempSpherical2SphericalInternalMMA(TestCase, SourceSpherical, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalTempSpherical2SphericalExternalMMA(TestCase, SourceSpherical, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalTempSpherical2WGS84InternalCore(TestCase, SourceSpherical, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalTempSpherical2WGS84InternalMMA(TestCase, SourceSpherical, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalTempSpherical2WGS84ExternalMMA(TestCase, SourceSpherical, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalTempWGS842CartesianInternalCore(TestCase, SourceGeodetic, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalTempWGS842CartesianInternalMMA(TestCase, SourceGeodetic, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN

class TestSHEvalTempWGS842CartesianExternalMMA(TestCase, SourceGeodetic, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_CARTESIAN


class TestSHEvalTempWGS842SphericalInternalCore(TestCase, SourceGeodetic, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalTempWGS842SphericalInternalMMA(TestCase, SourceGeodetic, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL

class TestSHEvalTempWGS842SphericalExternalMMA(TestCase, SourceGeodetic, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEOCENTRIC_SPHERICAL


class TestSHEvalTempWGS842WGS84InternalCore(TestCase, SourceGeodetic, SHTypeInternalCore, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalTempWGS842WGS84InternalMMA(TestCase, SourceGeodetic, SHTypeInternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

class TestSHEvalTempWGS842WGS84ExternalMMA(TestCase, SourceGeodetic, SHTypeExternalMMA, SphericalHarmonicsWithCoeffInterpolationMixIn):
    target_coordinate_system = GEODETIC_ABOVE_WGS84

#-------------------------------------------------------------------------------

class TestSHEvalTempCart2CartScaled(TestSHEvalTempCartesian2CartesianInternalCore):
    options = {"scale_potential": 2.0, "scale_gradient": [0.5, 1.0, -1.0]}
    scale_potential = 2.0
    scale_gradient = [0.5, 1.0, -1.0]

class TestSHEvalTempSph2SphScaled(TestSHEvalTempSpherical2SphericalInternalCore):
    options = {"scale_gradient": -1.0}
    scale_gradient = [-1.0, -1.0, -1.0]

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
