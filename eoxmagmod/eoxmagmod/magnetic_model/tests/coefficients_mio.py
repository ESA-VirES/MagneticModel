#-------------------------------------------------------------------------------
#
#  Swarm MIO_SHA_2* product spherical harmonic coefficients.
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
from numpy import inf
from numpy.testing import assert_allclose
from eoxmagmod._pytimeconv import decimal_year_to_mjd2000
from eoxmagmod.magnetic_model.coefficients_mio import SparseSHCoefficientsMIO
from eoxmagmod.magnetic_model.tests.data import SWARM_MIO_SHA_2_TEST_DATA
from eoxmagmod.magnetic_model.parser_mio import parse_swarm_mio_file


class MIOSHCoeffMixIn(object):
    is_internal = None
    degree = 2
    validity = (-inf, inf)

    def test_degree(self):
        self.assertEqual(self.coefficients.degree, self.degree)

    def test_validity(self):
        assert_allclose(self.coefficients.validity, self.validity)

    def test_is_internal(self):
        self.assertEqual(self.coefficients.is_internal, self.is_internal)

    def eval_coeff(self, time):
        return self.coefficients(
            time, lat_ngp=self.lat_ngp, lon_ngp=self.lon_ngp
        )


class TestSparseSHCoefficientsMIOInternal(TestCase, MIOSHCoeffMixIn):
    is_internal = True

    def setUp(self):
        with open(SWARM_MIO_SHA_2_TEST_DATA, "rb") as file_in:
            data = parse_swarm_mio_file(file_in)

        self.lat_ngp = data["lat_NGP"]
        self.lon_ngp = data["lon_NGP"]
        self.coefficients = SparseSHCoefficientsMIO(
            data["nm"], data["gh"],
            ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
            is_internal=True,
        )

    def test_callable(self):
        coeff, degree = self.eval_coeff(decimal_year_to_mjd2000(2018.5))
        assert_allclose(coeff, [
            (0, 0),
            (0.11415226, 0),
            (-0.22361224, -0.02624283),
            (-0.23244221, 0),
            (-0.34379927, -0.81524805),
            (0.01412196, -0.18049787),
        ], atol=1e-8)
        self.assertEqual(degree, self.degree)


class TestSparseSHCoefficientsMIOExternal(TestCase, MIOSHCoeffMixIn):
    is_internal = False

    def setUp(self):
        with open(SWARM_MIO_SHA_2_TEST_DATA, "rb") as file_in:
            data = parse_swarm_mio_file(file_in)

        self.lat_ngp = data["lat_NGP"]
        self.lon_ngp = data["lon_NGP"]
        self.coefficients = SparseSHCoefficientsMIO(
            data["nm"], data["qs"],
            ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
            is_internal=False,
        )

    def test_callable(self):
        coeff, degree = self.eval_coeff(decimal_year_to_mjd2000(2018.5))
        assert_allclose(coeff, [
            (0, 0),
            (-1.44080531, 0),
            (-0.27711333, -0.25918073),
            (-0.63645813, 0),
            (-0.22583255, -2.02305245),
            (0.05941826, -0.24358544),
        ], atol=1e-8)
        self.assertEqual(degree, self.degree)


if __name__ == "__main__":
    main()
