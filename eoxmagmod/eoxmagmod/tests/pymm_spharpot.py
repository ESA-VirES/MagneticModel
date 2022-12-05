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
# pylint: disable=missing-docstring

from unittest import TestCase, main
import operator
from itertools import chain, product
from functools import reduce
from random import random
from numpy import asarray, zeros, empty
from numpy.lib.stride_tricks import as_strided
from numpy.testing import assert_allclose
from eoxmagmod._pymm import legendre, lonsincos, relradpow, spharpot
from eoxmagmod.tests.data import sifm
from eoxmagmod.tests.data import mma_external


def _reshape_variable(shape, variable, preserved_dimensions=0):
    ndim = variable.ndim - preserved_dimensions
    if len(shape) > ndim:
        variable = as_strided(
            variable,
            shape=(*shape, *variable.shape[ndim:]),
            strides=(
                *variable.strides[:ndim],
                *((0,) * (len(shape) - ndim)),
                *variable.strides[ndim:],
            )
        )
    return variable


class SphericalHarmonicsCommonMixIn:
    is_internal = True
    degree = None
    coef_g = None
    coef_h = None

    @classmethod
    def reference_potential(cls, degree, coef_g, coef_h, latitude, longitude, radius):

        (
            shape, size, coef_g, coef_h, latitude, longitude, radius,
        ) = cls._ravel_inputs(degree, coef_g, coef_h, latitude, longitude, radius)

        result = empty(size)

        for i in range(size):
            result[i] = cls.reference_potential_scalar(
                degree, coef_g[i], coef_h[i],
                latitude[i], longitude[i], radius[i]
            )

        return result.reshape(shape)

    @classmethod
    def reference_potential_scalar(cls, degree, coef_g, coef_h, latitude, longitude, radius):
        rad_series, sin_series, cos_series, p_series, _ = cls.get_series(
            degree, latitude, longitude, radius,
        )
        m_idx = asarray(list(chain.from_iterable(
            range(n + 1) for n in range(degree + 1)
        )))
        n_idx = asarray(list(chain.from_iterable(
            [n]*(n+1)  for n in range(degree + 1)
        )))
        return (radius * rad_series[n_idx] * p_series * (
            coef_g * cos_series[m_idx] + coef_h * sin_series[m_idx]
        )).sum()


    @classmethod
    def eval_potential(cls, degree, coef_g, coef_h, latitude, longitude, radius):
        rad_series, sin_series, cos_series, p_series, _ = cls.get_series(
            degree, latitude, longitude, radius
        )
        return spharpot(
            radius, coef_g, coef_h, p_series, rad_series,
            sin_series, cos_series, degree=degree
        )

    @classmethod
    def get_series(cls, degree, latitude, longitude, radius):
        rad_series = relradpow(radius, degree, is_internal=cls.is_internal)
        sin_series, cos_series = lonsincos(longitude, degree)
        p_series, dp_series = legendre(latitude, degree)
        return rad_series, sin_series, cos_series, p_series, dp_series

    @classmethod
    def _ravel_inputs(cls, degree, coef_g, coef_h, latitude, longitude, radius):
        def _prod(values):
            return reduce(operator.mul, values, 1)

        coef_g = asarray(coef_g)
        coef_h = asarray(coef_h)
        latitude = asarray(latitude)
        longitude = asarray(longitude)
        radius = asarray(radius)

        shape = ()
        for array in [coef_g, coef_h]:
            if len(array.shape[:-1]) > len(shape):
                shape = array.shape[:-1]

        for array in [latitude, longitude, radius]:
            if len(array.shape) > len(shape):
                shape = array.shape

        size = _prod(shape)

        return (
            shape,
            size,
            _reshape_variable(shape, coef_g, 1).copy().reshape((size, coef_g.shape[-1])),
            _reshape_variable(shape, coef_h, 1).copy().reshape((size, coef_h.shape[-1])),
            _reshape_variable(shape, latitude).ravel(),
            _reshape_variable(shape, longitude).ravel(),
            _reshape_variable(shape, radius).ravel(),
        )


class SphericalHarmonicsPotentialTestMixIn(SphericalHarmonicsCommonMixIn):

    def test_coefficients(self):
        # This test checks response of the spherical harmonics to
        # individual coefficients set to 1 while the rest is kept 0.
        max_degree = 3
        coords = [
            (lat, lon, 6371.2) for lat, lon
            in product(range(-90, 91, 10), range(-180, 181, 20))
        ]

        for degree in range(max_degree + 1):
            size = ((degree + 1)*(degree + 2))//2
            offset = (degree*(degree + 1))//2

            for order in range(0, degree + 1):
                coef_g, coef_h = zeros(size), zeros(size)
                coef_g[order + offset] = 1.0

                for latitude, longitude, radius in coords:
                    assert_allclose(
                        self.eval_potential(
                            degree, coef_g, coef_h, latitude, longitude, radius
                        ),
                        self.reference_potential(
                            degree, coef_g, coef_h, latitude, longitude, radius
                        )
                    )

                coef_g, coef_h = zeros(size), zeros(size)
                coef_h[order + offset] = 1.0

                for latitude, longitude, radius in coords:
                    assert_allclose(
                        self.eval_potential(
                            degree, coef_g, coef_h, latitude, longitude, radius
                        ),
                        self.reference_potential(
                            degree, coef_g, coef_h, latitude, longitude, radius
                        )
                    )


    def test_potential(self):
        coords = [
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]
        for latitude, longitude, radius in coords:
            assert_allclose(
                self.eval_potential(
                    self.degree, self.coef_g, self.coef_h, latitude, longitude,
                    radius
                ),
                self.reference_potential(
                    self.degree, self.coef_g, self.coef_h, latitude, longitude,
                    radius
                ),
            )

    def test_potential_array_input(self):
        coords = asarray([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]).reshape((37, 37, 3))

        latitude = coords[..., 0]
        longitude = coords[..., 1]
        radius = coords[..., 2]

        assert_allclose(
            self.eval_potential(
                self.degree, self.coef_g, self.coef_h, latitude, longitude,
                radius
            ),
            self.reference_potential(
                self.degree, self.coef_g, self.coef_h, latitude, longitude,
                radius
            ),
        )


class TestSphericalHarmonicsPotentialInternal(TestCase, SphericalHarmonicsPotentialTestMixIn):
    is_internal = True
    degree = sifm.DEGREE
    coef_g = sifm.COEF_G
    coef_h = sifm.COEF_H


class TestSphericalHarmonicsPotentialExternal(TestCase, SphericalHarmonicsPotentialTestMixIn):
    is_internal = False
    degree = mma_external.DEGREE
    coef_g = mma_external.COEF_Q
    coef_h = mma_external.COEF_S


if __name__ == "__main__":
    main()
