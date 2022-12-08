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
from itertools import chain, product
from random import random
from math import pi, cos, sin, sqrt
from numpy import asarray, zeros, empty, stack
from numpy.testing import assert_allclose
from eoxmagmod._pymm import sphargrd
from eoxmagmod.tests.pymm_spharpot import SphericalHarmonicsCommonMixIn
from eoxmagmod.tests.data import sifm
from eoxmagmod.tests.data import mma_external


class SphericalHarmonicsGradientTestMixIn(SphericalHarmonicsCommonMixIn):

    @classmethod
    def eval_gradient(cls, degree, coeff, latitude, longitude, radius):
        rad_series, cos_sin_series, p_series, dp_series = cls.get_series(
            degree, latitude, longitude, radius
        )
        return sphargrd(
            latitude, coeff, p_series, dp_series, rad_series,
            cos_sin_series, is_internal=cls.is_internal,
            degree=degree
        )

    @classmethod
    def reference_gradient(cls, degree, coeff, latitude, longitude, radius):
        (
            shape, size, coeff, latitude, longitude, radius,
        ) = cls._ravel_inputs(degree, coeff, latitude, longitude, radius)

        result = empty((size, 3))

        for i in range(size):
            result[i] = cls.reference_gradient_scalar(
                degree, coeff[i], latitude[i], longitude[i], radius[i]
            )

        return result.reshape((*shape, 3))

    @classmethod
    def reference_gradient_scalar(cls, degree, coeff, latitude, longitude, radius):
        rad_series, cos_sin_series, p_series, dp_series = cls.get_series(
            degree, latitude, longitude, radius
        )
        m_idx = asarray(list(chain.from_iterable(
            range(n + 1) for n in range(degree + 1)
        )))
        n_idx = asarray(list(chain.from_iterable(
            [n]*(n+1)  for n in range(degree + 1)
        )))

        grad_lat = (rad_series[n_idx] * dp_series * (
            coeff[:, 0] * cos_sin_series[m_idx, 0] +
            coeff[:, 1] * cos_sin_series[m_idx, 1]
        )).sum()

        cos_lat = cos(latitude * pi / 180.0)
        if cos_lat > 1e-10:
            grad_lon = -(m_idx * rad_series[n_idx] * p_series * (
                coeff[:, 0] * cos_sin_series[m_idx, 1] -
                coeff[:, 1] * cos_sin_series[m_idx, 0]
            )).sum() / cos_lat
        else:
            sin_lat = sin(latitude * pi / 180.0)
            sin_lon = sin(longitude * pi / 180.0)
            cos_lon = cos(longitude * pi / 180.0)

            scale = [1.0]
            sqn1, ps1, ps0 = 1.0, 1.0, 1.0
            for i in range(2, degree + 1):
                # evaluate ratio between the Gauss-normalised and Schmidt
                # quasi-normalised associated Legendre functions.
                tmp = ((i-1)*(i-1)-1)/((2*i-1)*(2*i-3))
                ps1, ps0 = ps0, ps0*sin_lat - ps1*tmp
                sqn1 *= (2*i-1)/i
                scale.append(ps0 * sqn1 * sqrt((i*2)/(i+1)))
            scale = asarray(scale[:(degree + 1)])
            idx = asarray([
                1 + (i*(i + 1))//2 for i in range(1, degree + 1)
            ], dtype="int")
            grad_lon = -(scale * rad_series[1:] * (
                coeff[idx, 0]*sin_lon - coeff[idx, 1]*cos_lon
            )).sum()

        rad_scale = n_idx + 1 if cls.is_internal else -n_idx
        grad_rad = -(rad_scale * rad_series[n_idx] * p_series * (
            coeff[:, 0] * cos_sin_series[m_idx, 0] +
            coeff[:, 1] * cos_sin_series[m_idx, 1]
        )).sum()

        return asarray([grad_lat, grad_lon, grad_rad])

    def test_coefficients(self):
        max_degree = 3
        coords = [
            (lat, lon, 6371.2) for lat, lon
            in product(range(-90, 91, 10), range(-180, 180, 20))
        ]

        for degree in range(max_degree + 1):
            size = ((degree + 1)*(degree + 2))//2
            offset = (degree*(degree + 1))//2
            for order in range(0, degree + 1):
                coeff = zeros((size, 2))
                coeff[order + offset, 0] = 1.0

                for latitude, longitude, radius in coords:
                    assert_allclose(
                        self.eval_gradient(
                            degree, coeff, latitude, longitude, radius
                        ),
                        self.reference_gradient(
                            degree, coeff, latitude, longitude, radius
                        ),
                        atol=1e-14
                    )

                coeff = zeros((size, 2))
                coeff[order + offset, 1] = 1.0

                for latitude, longitude, radius in coords:
                    assert_allclose(
                        self.eval_gradient(
                            degree, coeff, latitude, longitude, radius
                        ),
                        self.reference_gradient(
                            degree, coeff, latitude, longitude, radius
                        ),
                        atol=1e-14
                    )

    def test_gradient(self):
        coords = [
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]

        for latitude, longitude, radius in coords[:1]:
            try:
                assert_allclose(
                    self.eval_gradient(
                        self.degree, self.coeff, latitude, longitude, radius
                    ),
                    self.reference_gradient(
                        self.degree, self.coeff, latitude, longitude, radius
                    ),
                )
            except AssertionError as exc:
                raise AssertionError(
                    f"point coordinates: ({latitude}, {longitude}, {radius})\n{exc}"
                ) from None
            break

    def test_gradient_array_input(self):
        coords = asarray([
            (lat, lon, 6371.2*(1.0 + random())) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ]).reshape((37, 37, 3))

        latitude = coords[..., 0]
        longitude = coords[..., 1]
        radius = coords[..., 2]

        assert_allclose(
            self.eval_gradient(
                self.degree, self.coeff, latitude, longitude, radius
            ),
            self.reference_gradient(
                self.degree, self.coeff, latitude, longitude, radius
            ),
        )

    def test_gradient_compared_with_finite_differences(self):
        # Compare gradient with the calculated finite differences.
        # The test fails if there in no match between the finite difference
        # approximation and the gradient evaluated by the spherical harmonics.
        eps = 0.1
        eps_deg = eps * (180.0 / (pi * 6371.2))
        radius = 6371.2

        def _compare_with_findiff(coord_centre, coord_lat, coord_lon, coord_rad):
            pot0 = self.eval_potential(
                self.degree, self.coeff, *coord_centre
            )
            pot_lat = self.eval_potential(
                self.degree, self.coeff, *coord_lat
            )
            pot_lon = self.eval_potential(
                self.degree, self.coeff, *coord_lon
            )
            pot_rad = self.eval_potential(
                self.degree, self.coeff, *coord_rad
            )

            grad_approx = asarray([
                (pot_lat - pot0)/eps, (pot_lon - pot0)/eps, (pot_rad - pot0)/eps
            ])
            grad_spharm = self.eval_gradient(
                self.degree, self.coeff, *coord_centre
            )

            try:
                assert_allclose(grad_spharm, grad_approx, rtol=1e-4, atol=1.0)
            except AssertionError as exc:
                latitude, longitude, radius = coord_centre
                raise AssertionError(
                    f"point coordinates: ({latitude}, {longitude}, {radius})\n{exc}"
                ) from None

        coords = list(product(range(-80, 81, 10), range(-180, 180, 20)))
        for latitude, longitude in coords:
            cos_lat = cos(latitude * pi / 180.0)
            _compare_with_findiff(
                (latitude, longitude, radius),
                (latitude + eps_deg, longitude, radius),
                (latitude, longitude + eps_deg / cos_lat, radius),
                (latitude, longitude, radius + eps)
            )

        for longitude in range(-180, 180, 20):
            _compare_with_findiff(
                (90, longitude, radius),
                (90 - eps_deg, longitude + 180, radius),
                (90 - eps_deg, longitude + 90, radius),
                (90, longitude, radius + eps),
            )

        for longitude in range(-180, 180, 20):
            _compare_with_findiff(
                (-90, longitude, radius),
                (-90 + eps_deg, longitude, radius),
                (-90 + eps_deg, longitude + 90, radius),
                (-90, longitude, radius + eps),
            )


class TestSphericalHarmonicsGradientInternal(TestCase, SphericalHarmonicsGradientTestMixIn):
    is_internal = True
    degree = sifm.DEGREE
    coeff = stack((sifm.COEF_G, sifm.COEF_H), axis=-1)



class TestSphericalHarmonicsGradientExternal(TestCase, SphericalHarmonicsGradientTestMixIn):
    is_internal = False
    degree = mma_external.DEGREE
    coeff = stack((mma_external.COEF_Q, mma_external.COEF_S), axis=-1)


if __name__ == "__main__":
    main()
