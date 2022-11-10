#-------------------------------------------------------------------------------
#
#  Bi-section interval search - tests
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
from numpy import searchsorted, nan, inf, linspace, isnan, asarray
from numpy.random import random
from numpy.testing import assert_equal
from eoxmagmod._pymm import (
    BISECT_SIDE_LEFT, BISECT_SIDE_RIGHT, bisect,
)


class BisectTestMixIn:

    def call_bisect(self, intervals, points):
        raise NotImplementedError

    def call_bisect_ref(self, intervals, points):
        raise NotImplementedError

    def _test_bisect(self, intervals, points):
        assert_equal(
            self.call_bisect(intervals, points),
            self.call_bisect_ref(intervals, points)
        )

    def test_special_values(self):
        self._test_bisect([0, 1], [nan, -inf, inf])

    def test_dimension_too_low(self):
        with self.assertRaises(ValueError):
            self.call_bisect(1, [])

    def test_dimension_too_high(self):
        with self.assertRaises(ValueError):
            self.call_bisect([[1, 2], [3, 4]], [])

    def test_border_values(self):
        self._test_bisect([1, 3, 5], [0, 1, 2, 3, 4, 5, 6])

    def test_zero_interval_values(self):
        self._test_bisect([1, 3, 3, 5], [0, 1, 2, 3, 4, 5, 6])

    def test_array_0d(self):
        self._test_bisect(linspace(0.1, 0.9, 5), random(()))

    def test_array_1d(self):
        self._test_bisect(linspace(0.1, 0.9, 5), random((5,)))

    def test_array_2d(self):
        self._test_bisect(linspace(0.1, 0.9, 5), random((2,3)))

#-------------------------------------------------------------------------------

class TestBisectLeft(TestCase, BisectTestMixIn):
    def call_bisect(self, intervals, points):
        return bisect(intervals, points, side=BISECT_SIDE_LEFT)

    def call_bisect_ref(self, intervals, points):
        points = asarray(points)
        idx = asarray(searchsorted(intervals, points, side="left") - 1)
        idx[isnan(points)] = -1
        return idx


class TestBisectRight(TestCase, BisectTestMixIn):
    def call_bisect(self, intervals, points):
        return bisect(intervals, points, side=BISECT_SIDE_RIGHT)

    def call_bisect_ref(self, intervals, points):
        points = asarray(points)
        idx = asarray(searchsorted(intervals, points, side="right") - 1)
        idx[isnan(points)] = -1
        return idx


class TestBisectDefault(TestBisectLeft):
    def call_bisect(self, intervals, points):
        return bisect(intervals, points)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
