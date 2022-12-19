#-------------------------------------------------------------------------------
#
#  Test of the MMA coefficients merge.
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
# The above copyright notice and this permission notice shall be included in all
# copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------
# pylint: disable=missing-docstring,attribute-defined-outside-init

from unittest import TestCase, main
from numpy import asarray, concatenate
from numpy.random import random
from numpy.testing import assert_equal
from eoxmagmod.magnetic_model.loader_mma import _merge_coefficients


class MMAMergeTest(TestCase):
    variable = "data"
    set1 = {
        "degree_min": 1,
        "degree_max": 2,
        "nm": asarray([(1, 0), (1, 1), (2, -1)]),
        "t": asarray([0.0, 0.5, 1.0, 1.5, 2.0]),
        "data": random((3, 5)),
    }
    set2 = {
        "degree_min": 1,
        "degree_max": 2,
        "nm": asarray([(1, 0), (1, 1), (2, -1)]),
        "t": asarray([2.0, 2.5, 3.0, 3.5, 4.0]),
        "data": random((3, 5)),
    }
    set12_merged = {
        **set1,
        **{
            "t": concatenate((set1["t"], set2["t"])),
            "data": concatenate((set1["data"], set2["data"]), axis=1),
        }
    }
    set2_min_degree_mismatch = {
        "degree_min": 2,
        "degree_max": 2,
        "nm": asarray([(1, 0), (1, 1), (2, -1)]),
        "t": asarray([2.0, 2.5, 3.0, 3.5, 4.0]),
        "data": random((3, 5)),
    }
    set2_max_degree_mismatch = {
        "degree_min": 1,
        "degree_max": 1,
        "nm": asarray([(1, 0), (1, 1), (2, -1)]),
        "t": asarray([2.0, 2.5, 3.0, 3.5, 4.0]),
        "data": random((3, 5)),
    }
    set2_nm_mismatch = {
        "degree_min": 1,
        "degree_max": 2,
        "nm": asarray([(1, 0), (1, 1), (1, -1)]),
        "t": asarray([2.0, 2.5, 3.0, 3.5, 4.0]),
        "data": random((3, 5)),
    }

    @staticmethod
    def _test_merge(tested, reference):
        assert_equal(tested["degree_min"], reference["degree_min"])
        assert_equal(tested["degree_max"], reference["degree_max"])
        assert_equal(tested["nm"], reference["nm"])
        assert_equal(tested["t"], reference["t"])
        assert_equal(tested["data"], reference["data"])

    def test_merge_single(self):
        self._test_merge(
            _merge_coefficients(self.set1, variable=self.variable),
            self.set1,
        )

    def test_merge_multiple(self):
        self._test_merge(
            _merge_coefficients(self.set1, self.set2, variable=self.variable),
            self.set12_merged
        )

    def test_min_degree_mismatch(self):
        with self.assertRaises(ValueError):
            _merge_coefficients(
                self.set1, self.set2_min_degree_mismatch,
                variable=self.variable
            )

    def test_max_degree_mismatch(self):
        with self.assertRaises(ValueError):
            _merge_coefficients(
                self.set1, self.set2_max_degree_mismatch,
                variable=self.variable
            )

    def test_mn_degree_mismatch(self):
        with self.assertRaises(ValueError):
            _merge_coefficients(
                self.set1, self.set2_nm_mismatch, variable=self.variable
            )

    def test_unordered_time(self):
        with self.assertRaises(ValueError):
            _merge_coefficients(self.set2, self.set1, variable=self.variable)


if __name__ == "__main__":
    main()
