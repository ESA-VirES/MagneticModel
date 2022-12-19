#-------------------------------------------------------------------------------
#
#  Coefficients evaluated by 2D Fourier series - tests
#
# Author: Martin Paces <martin.paces@eox.at>
#
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
# pylint: disable=missing-docstring, invalid-name, too-few-public-methods

from unittest import TestCase, main
from numpy import pi, empty, arange, broadcast_to, sin, cos
from numpy.random import uniform
from numpy.testing import assert_allclose
from eoxmagmod._pymm import fourier2d
from eoxmagmod.tests.data import mio


class Fourier2DTestMixIn:
    coeff = mio.COEFF_I
    degrees = ((mio.PMIN, mio.PMAX), (mio.SMIN, mio.SMAX))
    scales = (1.0, 1.0)
    value_ranges = ((0.0, 2*pi), (0.0, 2*pi))

    def _get_xy(self, shape):
        x_range, y_range = self.value_ranges
        return uniform(*x_range, shape), uniform(*y_range, shape)

    def _test_fourier2d(self, x, y):
        result = self.eval_fourier2d(x, y)
        result_ref = self.eval_reference_fourier2d(x, y)
        assert_allclose(result, result_ref, atol=1e-15, rtol=1e-15)


    @classmethod
    def eval_fourier2d(cls, x, y):
        (min_degree1, _), (min_degree2, _) = cls.degrees
        scale1, scale2 = cls.scales
        return fourier2d(
            y, x, cls.coeff, min_degree2, min_degree1, scale2, scale1
        )

    @classmethod
    def eval_reference_fourier2d(cls, x, y):
        return cls._eval_coeff_fourier2d(x, y)

    @classmethod
    def _eval_coeff_fourier2d(cls, x, y):
        assert x.shape == y.shape
        shape = x.shape
        x, y = x.ravel(), y.ravel()
        result = empty((x.size, cls.coeff.shape[0]))
        for i in range(x.size):
            result[i, :] = cls._eval_coeff_fourier2d_single_point(x[i], y[i])
        return  result.reshape((*shape, result.shape[1]))

    @classmethod
    def _eval_coeff_fourier2d_single_point(cls, x, y):
        """ Function using the 2D Fourier series. """
        coeff = cls.coeff
        x_scale, y_scale = cls.scales
        sin_xy, cos_xy = cls._calculate_sincos_matrices(x_scale * x, y_scale * y)
        sin_xy = broadcast_to(sin_xy, (coeff.shape[0], *sin_xy.shape))
        cos_xy = broadcast_to(cos_xy, (coeff.shape[0], *cos_xy.shape))
        return (coeff[..., 0]*cos_xy + coeff[..., 1]*sin_xy).sum(axis=(1, 2))

    @classmethod
    def _calculate_sincos_matrices(cls, x, y):
        """ Get sin/cos matrices used by the 2D Fourier transform. """
        (n_min, n_max), (m_min, m_max) = cls.degrees
        xi = x * arange(n_min, n_max + 1)
        yj = y * arange(m_min, m_max + 1)
        n_col = n_max - n_min + 1
        n_row = m_max - m_min + 1
        xi_yj = (
            broadcast_to(xi, (n_row, n_col)) +
            broadcast_to(yj, (n_col, n_row)).transpose()
        )
        return sin(xi_yj), cos(xi_yj)

    def test_eval_0d(self):
        self._test_fourier2d(*self._get_xy(()))

    def test_eval_1d(self):
        self._test_fourier2d(*self._get_xy((10)))

    def test_eval_2d(self):
        self._test_fourier2d(*self._get_xy((4, 3)))


class TestFourier2D(TestCase, Fourier2DTestMixIn):
    pass


class TestFourier2DScaled(TestCase, Fourier2DTestMixIn):
    scales = (pi, pi)
    value_ranges = ((-1.0, 1.0), (-1.0, 1.0))


if __name__ == "__main__":
    main()
