#-------------------------------------------------------------------------------
#
#  Time utilities
#
# Project: VirES
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
# all
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

from datetime import datetime
from numpy import asarray, floor
from ._pytimeconv import (
    decimal_year_to_mjd2000,
    mjd2000_to_decimal_year,
    mjd2000_to_year_fraction,
)

__all__ = [
    "decimal_year_to_mjd2000",
    "mjd2000_to_decimal_year",
    "mjd2000_to_year_fraction",
    "mjd2000_to_decimal_year_simple",
    "mjd2000_to_year_fraction_simple",
    "decimal_year_to_mjd2000_simple",
    "datetime_to_decimal_year",
]


def mjd2000_to_decimal_year_simple(mjd2000):
    """ Convert Modified Julian Date since 2000 to (Julian) decimal year
    using the simplified conversion method:
        decimal_year = 2000.0 + mjd2000/365.25
    """
    return 2000.0 + asarray(mjd2000) / 365.25


def mjd2000_to_year_fraction_simple(mjd2000):
    """ Convert Modified Julian Date since 2000 to (Julian) year fraction
    using the simplified conversion method:
        year_fraction = mjd2000/365.25 - floor(mjd2000/365.25)
    """
    mjy2000 = asarray(mjd2000) / 365.25
    return mjy2000 - floor(mjy2000)


def decimal_year_to_mjd2000_simple(decimal_year):
    """ Convert (Julian) decimal year to Modified Julian Date since 2000
    using the simplified conversion method:
        mjd2000 = (decimal_year - 2000.0) * 365.25
    """
    return (asarray(decimal_year) - 2000.0) * 365.25


def datetime_to_decimal_year(time):
    """ Convert time given by a `datetime.datetime` object to a decimal year
    value.
    """
    if not isinstance(time, datetime):
        raise TypeError("The input must be a datetime object.")

    year_start = datetime(year=time.year, month=1, day=1)
    next_year_start = datetime(year=time.year+1, month=1, day=1)

    year_elapsed = (time - year_start).total_seconds()
    year_total = (next_year_start - year_start).total_seconds()

    return time.year + year_elapsed / year_total
