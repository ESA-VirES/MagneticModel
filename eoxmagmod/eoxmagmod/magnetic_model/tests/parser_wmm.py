#-------------------------------------------------------------------------------
#
#  World Magnetic Model file format parser - test
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
# pylint: disable=missing-docstring

from unittest import TestCase, main

try:
    # Python 2
    from StringIO import StringIO
except ImportError:
    from io import StringIO

from numpy.testing import assert_allclose
from numpy import abs as aabs
from eoxmagmod.magnetic_model.parser_wmm import (
    parse_wmm_file, WMM_VALIDITY_PERIOD,
)
from eoxmagmod.data import WMM_2015


class TestWMMParser(TestCase):

    @staticmethod
    def parse(filename):
        with open(filename, encoding="utf8") as file_in:
            return parse_wmm_file(file_in)

    def _assert_valid(self, data, expected_data):
        tested_data = {
            key: data[key] for key in expected_data
        }
        assert_allclose(
            data["t"], [data["epoch"], data["epoch"] + WMM_VALIDITY_PERIOD]
        )
        self.assertEqual(tested_data, expected_data)
        self.assertEqual(data["t"].size, data["gh"].shape[1])
        self.assertEqual(data["nm"].shape[0], data["gh"].shape[0])
        self.assertEqual(data["nm"].shape[1], 2)
        self.assertEqual(data["nm"][..., 0].min(), data["degree_min"])
        self.assertEqual(data["nm"][..., 0].max(), data["degree_max"])
        self.assertTrue(aabs(data["nm"][..., 1]).max() <= data["degree_max"])

    def test_parse_wmm_test_sample(self):
        data = parse_wmm_file(StringIO(TEST_PRODUCT))
        assert_allclose(data["nm"], [
            (10, 1), (10, -2), (11, 3), (11, -4),
            (12, 5), (12, -6), (13, 7), (13, -7),
        ])
        assert_allclose(data["gh"], [
            (1, 1), (1, 1), (0, 5), (0, 5),
            (1, 6), (1, 6), (1, 6), (1, 6),
        ])
        self._assert_valid(data, {
            "name": "WMM-9999",
            "version": "DD/MM/YYYY",
            "epoch": 9999.0,
            "degree_min": 10,
            "degree_max": 13,
        })

    def test_parse_wmm_file_wmm2015(self):
        data = self.parse(WMM_2015)
        self._assert_valid(data, {
            "name": "WMM-2015v2",
            "version": "09/18/2018",
            "epoch": 2015.0,
            "degree_min": 1,
            "degree_max": 12,
        })


TEST_PRODUCT = ("""\
9999.0            WMM-9999        DD/MM/YYYY
9 0 0.0 1.0 0.0 1.0
9 1 0.0 0.0 0.0 0.0
10 1 1.0 0.0 0.0 0.0
10 2 0.0 1.0 0.0 0.0
11 3 0.0 0.0 1.0 0.0
11 4 0.0 0.0 0.0 1.0
12 5 1.0 0.0 1.0 0.0
12 6 0.0 1.0 0.0 1.0
13 7 1.0 1.0 1.0 1.0
999999999999999999999999999999999999999999999999
999999999999999999999999999999999999999999999999
""")


if __name__ == "__main__":
    main()
