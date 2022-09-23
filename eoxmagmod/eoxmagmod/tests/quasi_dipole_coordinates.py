#-------------------------------------------------------------------------------
#
#  Quasi-Dipole coordinates - test
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
#pylint: disable=missing-docstring

from unittest import TestCase, main
from numpy import array
from numpy.testing import assert_allclose
from eoxmagmod.quasi_dipole_coordinates import (
    eval_qdlatlon, eval_mlt, eval_subsol,
    eval_qdlatlon_with_base_vectors,
)
from eoxmagmod.tests.data import QUASI_DIPOLE_TEST_DATA


def load_test_data(filename):
    """ Load test data from a tab-separated values file. """
    def _load_test_data(file_in):
        header = next(file_in).strip().split("\t")
        records = array([
            [float(v) for v in line.strip().split("\t")] for line in file_in
        ])
        return {
            variable: records[..., idx] for idx, variable in enumerate(header)
        }

    with open(filename, encoding="ascii") as file_in:
        return _load_test_data(file_in)


class TestQuasiDipoleCoordinates(TestCase):
    test_data = load_test_data(QUASI_DIPOLE_TEST_DATA)

    def test_eval_qdlatlon(self):
        qdlat, qdlon = eval_qdlatlon(
            self.test_data["Latitude"],
            self.test_data["Longitude"],
            self.test_data["Radius"],
            self.test_data["DecimalYear"],
        )
        assert_allclose(qdlat, self.test_data["QDLatitude"], atol=1e-8)
        assert_allclose(qdlon, self.test_data["QDLongitude"], atol=1e-8)

    def test_eval_qdlatlon_with_base_vectors(self):
        qdlat, qdlon, f11, f12, f21, f22, f__ = eval_qdlatlon_with_base_vectors(
            self.test_data["Latitude"],
            self.test_data["Longitude"],
            self.test_data["Radius"],
            self.test_data["DecimalYear"],
        )
        assert_allclose(qdlat, self.test_data["QDLatitude"], atol=1e-8)
        assert_allclose(qdlon, self.test_data["QDLongitude"], atol=1e-8)
        assert_allclose(f11, self.test_data["F11"], atol=1e-8)
        assert_allclose(f12, self.test_data["F12"], atol=1e-8)
        assert_allclose(f21, self.test_data["F21"], atol=1e-8)
        assert_allclose(f22, self.test_data["F22"], atol=1e-8)
        assert_allclose(f__, self.test_data["F"], atol=1e-8)

    def test_eval_mlt(self):
        mlt = eval_mlt(self.test_data["QDLongitude"], self.test_data["MJD2000"])
        assert_allclose(mlt, self.test_data["MagneticLocalTime"], atol=1e-8)

    def test_eval_subsol(self):
        sollat, sollon = eval_subsol(self.test_data["MJD2000"])
        assert_allclose(sollat, self.test_data["SubsolarLatitude"], atol=1e-8)
        assert_allclose(sollon, self.test_data["SubsolarLongitude"], atol=1e-8)


if __name__ == "__main__":
    main()
