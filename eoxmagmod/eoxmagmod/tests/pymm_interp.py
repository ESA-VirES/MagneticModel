#-------------------------------------------------------------------------------
#
#  Coefficients time interpolation - tests
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
from numpy import (
    nan, inf, linspace,
    digitize, asarray, prod, full, empty,
)
from numpy.random import random
from numpy.testing import assert_equal
from eoxmagmod._pymm import (
    interp, INTERP_C1, INTERP_C1D1,
)


class InterpTestMixIn():

    def call_interp(self, time, time0, coeff0):
        raise NotImplementedError

    def call_interp_ref(self, time, time0, coeff0):
        raise NotImplementedError

    def _test_interp(self, time, time0, coeff0):
        assert_equal(
            self.call_interp(time, time0, coeff0),
            self.call_interp_ref(time, time0, coeff0)
        )

    def test_special_values(self):
        self._test_interp([nan, -inf, inf], [0, 1], [1, -1])

    def test_dimension_mismatch1(self):
        with self.assertRaises(ValueError):
            self._test_interp([], [1, 3, 5], [2, 7])

    def test_dimension_mismatch2(self):
        with self.assertRaises(ValueError):
            self._test_interp([], [1, 3, 5], [[2, 7, 3, 5]])

    def test_wrong_time_dimenion1(self):
        with self.assertRaises(ValueError):
            self._test_interp([], 1, [2])

    def test_wrong_coeff_dimenion1(self):
        with self.assertRaises(ValueError):
            self._test_interp([], [1], 2)

    def test_wrong_time_dimenion2(self):
        with self.assertRaises(ValueError):
            self._test_interp([], [[1, 3, 5]], [[2, 7, 3]])

    def test_border_values(self):
        self._test_interp(
            [0, 1, 2, 3, 4, 5, 6],
            [1, 3, 5],
            [2, 7, 3]
        )

    def test_dicontinuity(self):
        self._test_interp(
            [0, 1, 2, 3, 4, 5, 6],
            [1, 3, 3, 5],
            [2, 7, 6, 3]
        )

    def test_shape_scalar_time_2d_coeff(self):
        self._test_interp(
            1,
            [1, 3, 5],
            [[[2, 7, 3], [7, 3, 2]], [[3, 2, 7], [2, 7, 3]]]
        )

    def test_shape_empty_1d_time_2d_coeff(self):
        self._test_interp(
            [],
            [1, 3, 5],
            [[[2, 7, 3], [7, 3, 2]], [[3, 2, 7], [2, 7, 3]]]
        )

    def test_shape_2d_time_2d_coeff(self):
        self._test_interp(
            [[0, 1, 2,], [4, 5, 6]],
            [1, 3, 5],
            [[[2, 7, 3], [7, 3, 2]], [[3, 2, 7], [2, 7, 3]]]
        )

    def test_shape_2d_time_1d_coeff(self):
        self._test_interp(
            [[0, 1, 2,], [4, 5, 6]],
            [1, 3, 5],
            [[2, 7, 3], [7, 3, 2]]
        )

    def test_shape_2d_time_empty_1d_coeff(self):
        self._test_interp(
            [[0, 1, 2,], [4, 5, 6]],
            [1, 3, 5],
            empty((0, 3)) # note the last dimension matching the times
        )

    def test_shape_2d_time_scalar_coeff(self):
        self._test_interp(
            [[0, 1, 2,], [4, 5, 6]],
            [1, 3, 5],
            [2, 7, 3],
        )

    def test_array_0d(self):
        times = linspace(0.1, 0.9, 5)
        coeff = random((3, times.size))
        self._test_interp(random(()), times, coeff)

    def test_array_1d(self):
        times = linspace(0.1, 0.9, 5)
        coeff = random((3, times.size))
        self._test_interp(random((5,)), times, coeff)

    def test_array_2d(self):
        times = linspace(0.1, 0.9, 5)
        coeff = random((3, times.size))
        self._test_interp(random((2, 3)), times, coeff)

#-------------------------------------------------------------------------------

class TestInterpLinear(TestCase, InterpTestMixIn):
    def call_interp(self, time, time0, coeff0):
        return interp(time, time0, coeff0, kind=INTERP_C1)

    def call_interp_ref(self, time, time0, coeff0):
        return interpolate_c1(time, time0, coeff0)


class TestInterpDefault(TestInterpLinear):
    def call_interp(self, time, time0, coeff0):
        return interp(time, time0, coeff0)


class TestInterpLinearD1(TestCase, InterpTestMixIn):
    def call_interp(self, time, time0, coeff0):
        return interp(time, time0, coeff0, kind=INTERP_C1D1)

    def call_interp_ref(self, time, time0, coeff0):
        return interpolate_c1d1(time, time0, coeff0)

#-------------------------------------------------------------------------------
# reference Python implementation

def interpolate_c1(time, time0, coeff0):
    def get_basis(time, time0, idx):
        alpha = (time - time0[idx]) / (time0[idx+1] - time0[idx])
        return (1.0 - alpha, alpha)
    return _interpolate_c1(time, time0, coeff0, get_basis)


def interpolate_c1d1(time, time0, coeff0):
    def get_basis(time, time0, idx):
        del time
        alpha = 1.0 / (time0[idx+1] - time0[idx])
        return (-alpha, alpha)
    return _interpolate_c1(time, time0, coeff0, get_basis)


def _interpolate_c1(time, time0, coeff0, get_basis):
    time, time0, coeff0 = to_array(time, time0, coeff0)
    time, coeff0, t_shape, c_shape = _simplify_shape(time, coeff0)

    if time0.size < 2:
        return full(t_shape + c_shape, nan)

    # lookup spline segment
    idx0 = find_interval(time, time0)
    idx = clip(idx0.copy(), 0, time0.size-2)

    # get basis and perform its products with the coefficients
    basis0, basis1 = get_basis(time, time0, idx)
    coeff = basis0 * coeff0[:, idx] + basis1 * coeff0[:, idx+1]

    # clear out of bounds values
    coeff[:, (idx0 != idx) & (time != time0[-1])] = nan

    return _restore_shape(coeff, t_shape, c_shape)


def _restore_shape(coeff, t_shape, c_shape):
    return asarray(coeff.T, order="C").reshape(t_shape + c_shape)


def _simplify_shape(time, coeff0):
    t_shape = time.shape
    time = time.ravel()

    c_shape = coeff0.shape[:-1]
    coeff0 = coeff0.reshape(
        (prod(c_shape) if c_shape else 1, coeff0.shape[-1])
    )

    return  time, coeff0, t_shape, c_shape


def to_array(*args):
    return tuple(asarray(arg, order='C') for arg in args)


def find_interval(time, time0):
    return digitize(time, time0) - 1


def clip(values, value_min, value_max):
    values[values < value_min] = value_min
    values[values > value_max] = value_max
    return values

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
