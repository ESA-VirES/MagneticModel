#-------------------------------------------------------------------------------
#
#  Solar position - test
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
from eoxmagmod.quasi_dipole_coordinates import eval_subsol
from eoxmagmod.solar_position import sunpos, sunpos_original
from eoxmagmod.tests.data import SUN_POSITION_TEST_DATA


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


class TestSunPosition(TestCase):
    test_data = load_test_data(SUN_POSITION_TEST_DATA)

    def test_sunpos(self):
        decl, rasc, lha, azimuth, zenith = sunpos(
            self.test_data["MJD2000"],
            self.test_data["Latitude"],
            self.test_data["Longitude"],
            self.test_data["Radius"],
        )
        assert_allclose(decl, self.test_data["Declination"], atol=1e-8)
        assert_allclose(rasc, self.test_data["RightAscension"], atol=1e-8)
        assert_allclose(lha, self.test_data["HourAngle"], atol=1e-8)
        assert_allclose(azimuth, self.test_data["Azimuth"], atol=1e-8)
        assert_allclose(zenith, self.test_data["Zenith"], atol=1e-8)

    def test_sunpos_original(self):
        # Note: the original code does handle correctly only days
        #       between 2000-03-01 and 2400-02-29.
        mask = (
            (self.test_data["MJD2000"] >= 60.0) &
            (self.test_data["MJD2000"] < 146157.0)
        )
        decl, rasc, lha, azimuth, zenith = sunpos_original(
            self.test_data["MJD2000"][mask],
            self.test_data["Latitude"][mask],
            self.test_data["Longitude"][mask],
            self.test_data["Radius"][mask],
            pres=0.0,
        )
        assert_allclose(decl, self.test_data["Declination"][mask], atol=1e-8)
        assert_allclose(rasc, self.test_data["RightAscension"][mask], atol=1e-8)
        assert_allclose(lha, self.test_data["HourAngle"][mask], atol=1e-8)
        assert_allclose(azimuth, self.test_data["Azimuth"][mask], atol=1e-8)
        assert_allclose(zenith, self.test_data["Zenith"][mask], atol=1e-8)


    def test_sunpos_subsol_comparison(self):
        # The two different Solar models are expected to be aligned.
        decl, _, gha, _, _, = sunpos(self.test_data["MJD2000"], 0, 0, 0)
        sollat, sollon = eval_subsol(self.test_data["MJD2000"])
        assert_allclose(decl, sollat, atol=4e-1)
        assert_allclose(-gha, sollon, atol=2e-1)


if __name__ == "__main__":
    main()
