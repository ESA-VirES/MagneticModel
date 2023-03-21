#-------------------------------------------------------------------------------
#
#  Time conversion utilities - test
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
# pylint: disable=missing-docstring, invalid-name, too-few-public-methods

from unittest import TestCase, main
from datetime import date, datetime
from numpy import vectorize, inf, nan
from numpy.random import uniform
from numpy.testing import assert_allclose
from eoxmagmod.tests.util import FunctionTestMixIn
from eoxmagmod.time_util import (
    decimal_year_to_mjd2000_simple,
    mjd2000_to_decimal_year_simple,
    mjd2000_to_year_fraction_simple,
    datetime_to_decimal_year,
)


class TestDatetimeToDecimalYear(FunctionTestMixIn, TestCase):
    NAME = "datetime_to_decimal_year"

    @staticmethod
    def eval(input_):
        return datetime_to_decimal_year(input_)

    ACCEPTED = [
        (datetime(2001, 1, 12), 2001.0301369863014),
        (datetime(2012, 8, 31), 2012.6639344262296),
        (datetime(2014, 8, 31), 2014.66301369863),
        (datetime(2024, 12, 31, 23, 59, 59, 999), 2024.9999999684085),
    ]

    REJECTED = [
        (None, TypeError),
        (date(2001, 1, 12), TypeError),
    ]


class TestMjd2000ToYearFractionSimple(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_mjd2000_to_year_fraction_simple)(value)

    @staticmethod
    def eval(value):
        return mjd2000_to_year_fraction_simple(value)

    @staticmethod
    def _assert(tested, expected):
        assert_allclose(tested, expected, rtol=1e-14, atol=1e-11)

    def test_mjd2000_to_year_fraction_far_range(self):
        values = uniform(-730487., 730485., (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_mjd2000_to_year_fraction_near_range(self):
        values = uniform(-36524., 36525., (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_mjd2000_to_year_fraction_sanity_check(self):
        self._assert(self.eval(0.), 0.0)
        self._assert(self.eval([-365.25, 365.25]), [0.0, 0.0])
        self._assert(self.eval([6757.125, 7487.625]), [0.5, 0.5])


class TestMjd2000ToDecimalYearSimple(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_mjd2000_to_decimal_year_simple)(value)

    @staticmethod
    def eval(value):
        return mjd2000_to_decimal_year_simple(value)

    @staticmethod
    def _assert(tested, expected):
        assert_allclose(tested, expected, rtol=1e-14, atol=1e-11)

    def test_mjd2000_to_decimal_year_far_range(self):
        values = uniform(-730487., 730485., (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_mjd2000_to_decimal_year_near_range(self):
        values = uniform(-36524., 36525., (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_mjd2000_to_decimal_year_sanity_check(self):
        self._assert(self.eval(0.), 2000.0)
        self._assert(self.eval([-365.25, 365.25]), [1999.0, 2001])
        self._assert(self.eval([6757.125, 7487.625]), [2018.5, 2020.5])

    def test_mjd2000_to_decimal_year_special_values(self):
        self._assert(self.eval([-inf, inf, nan]), [-inf, inf, nan])


class TestDecimalYearToMjd2000Simple(TestCase):

    @staticmethod
    def reference(value):
        return vectorize(_decimal_year_to_mjd2000_simple)(value)

    @staticmethod
    def eval(value):
        return decimal_year_to_mjd2000_simple(value)

    @staticmethod
    def _assert(tested, expected):
        assert_allclose(tested, expected, rtol=1e-14, atol=6e-8)

    def test_decimal_year_to_mjd2000_far_range(self):
        values = uniform(0.0, 4000.0, (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_decimal_year_to_mjd2000_near_range(self):
        values = uniform(1900, 2100.0, (100, 100))
        self._assert(self.eval(values), self.reference(values))

    def test_decimal_year_to_mjd2000_sanity_check(self):
        self._assert(self.eval(2000.), 0.0)
        self._assert(self.eval([1999., 2001.0]), [-365.25, 365.25])
        self._assert(self.eval([2018.5, 2020.5]), [6757.125, 7487.625])

    def test_decimal_year_to_mjd2000_special_values(self):
        self._assert(self.eval([-inf, inf, nan]), [-inf, inf, nan])

#-------------------------------------------------------------------------------
# reference implementation

def _mjd2000_to_year_fraction_simple(mjd2k):
    """ Convert Modified Julian Date 2000 to (Julian) year fraction.
    """
    return (mjd2k / 365.25) % 1


def _mjd2000_to_decimal_year_simple(mjd2k):
    """ Convert Modified Julian Date 2000 to decimal year.
    """
    return 2000.0 + mjd2k / 365.25


def _decimal_year_to_mjd2000_simple(decimal_year):
    """ Covert decimal year to Modified Julian Date 2000.
    """
    return (decimal_year - 2000.0) * 365.25

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
