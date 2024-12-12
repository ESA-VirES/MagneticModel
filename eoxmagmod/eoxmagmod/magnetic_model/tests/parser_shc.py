#-------------------------------------------------------------------------------
#
#  SHC format parser - test
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
# pylint: disable=missing-docstring,invalid-name

from unittest import TestCase, main
from numpy import abs as aabs
from numpy.testing import assert_equal
from eoxmagmod.magnetic_model.parser_shc import parse_shc_file, parse_shc_header
from eoxmagmod.data import (
    CHAOS_CORE_LATEST, CHAOS_CORE_PREDICTION_LATEST, CHAOS_STATIC_LATEST,
    IGRF12, IGRF13, IGRF14, IGRF_LATEST, SIFM, LCS1, MF7,
)


class TestSHCParser(TestCase):

    @staticmethod
    def parse(filename):
        with open(filename, encoding="utf8") as file_in:
            return parse_shc_file(file_in)

    @staticmethod
    def parse_header(filename):
        with open(filename, encoding="utf8") as file_in:
            return parse_shc_header(file_in)

    def _assert_valid(self, data, expected_data):
        tested_data = {
            key: data[key] for key in expected_data
        }
        self.assertEqual(tested_data, expected_data)
        self.assertEqual(data["t"].size, data["gh"].shape[1])
        self.assertEqual(data["nm"].shape[0], data["gh"].shape[0])
        self.assertEqual(data["nm"].shape[1], 2)
        self.assertEqual(data["nm"][..., 0].min(), data["degree_min"])
        self.assertEqual(data["nm"][..., 0].max(), data["degree_max"])
        self.assertTrue(aabs(data["nm"][..., 1]).max() <= data["degree_max"])

    def _test_header(self, filename):
        reference = self.parse(filename)
        tested = self.parse_header(filename)
        assert_equal(reference["t"], tested["t"])
        for remove_key in ["nm", "gh", "t"]:
            reference.pop(remove_key)
        tested.pop("t")
        self.assertEqual(reference, tested)

    def test_parse_shc_file_sifm(self):
        data = self.parse(SIFM)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 70,
            "spline_order": 2,
            "ntime": 2,
            "nstep": 1,
        })

    def test_parse_shc_header_sifm(self):
        self._test_header(SIFM)

    def test_parse_shc_file_igrf12(self):
        data = self.parse(IGRF12)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 13,
            "spline_order": 2,
            "ntime": 25,
            "nstep": 1,
        })

    def test_parse_shc_file_igrf13(self):
        data = self.parse(IGRF13)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 13,
            "spline_order": 2,
            "ntime": 26,
            "nstep": 1,
        })

    def test_parse_shc_file_igrf14(self):
        data = self.parse(IGRF14)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 13,
            "spline_order": 2,
            "ntime": 27,
            "nstep": 1,
        })

    def test_parse_shc_file_igrf_last(self):
        data = self.parse(IGRF_LATEST)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 13,
            "spline_order": 2,
            "ntime": 27,
            "nstep": 1,
        })

    def test_parse_shc_header_igrf12(self):
        self._test_header(IGRF12)

    def test_parse_shc_file_chaos_core_latest(self):
        data = self.parse(CHAOS_CORE_LATEST)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 20,
            "spline_order": 6,
            "ntime": 281,
            "nstep": 5,
        })

    def test_parse_shc_header_chaos_core_latest(self):
        self._test_header(CHAOS_CORE_LATEST)

    def test_parse_shc_file_chaos_core_prediction_latest(self):
        data = self.parse(CHAOS_CORE_PREDICTION_LATEST)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 20,
            "spline_order": 2,
            "ntime": 2,
            "nstep": 1,
        })

    def test_parse_shc_header_chaos_core_prediction_latest(self):
        self._test_header(CHAOS_CORE_PREDICTION_LATEST)

    def test_parse_shc_file_chaos_static(self):
        data = self.parse(CHAOS_STATIC_LATEST)
        self._assert_valid(data, {
            "degree_min": 21,
            "degree_max": 185,
            "spline_order": 1,
            "ntime": 1,
            "nstep": 0,
        })

    def test_parse_shc_header_chaos_static(self):
        self._test_header(CHAOS_STATIC_LATEST)

    def test_parse_shc_file_lcs1(self):
        data = self.parse(LCS1)
        self._assert_valid(data, {
            "degree_min": 1,
            "degree_max": 185,
            "spline_order": 1,
            "ntime": 1,
            "nstep": 1,
        })

    def test_parse_shc_file_mf7(self):
        data = self.parse(MF7)
        self._assert_valid(data, {
            "degree_min": 16,
            "degree_max": 133,
            "spline_order": 1,
            "ntime": 1,
            "nstep": 1,
        })


if __name__ == "__main__":
    main()
