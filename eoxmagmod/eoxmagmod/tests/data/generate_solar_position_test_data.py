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

import sys
from itertools import product
from numpy import array
from numpy.random import uniform
from eoxmagmod.solar_position import sunpos

EARTH_RADIUS = 6371.2   # km
START_TIME = -10957.0   # MJD2000 / 1970-01-01T00:00:00Z
END_TIME = 10958.0      # MJD2000 / 2030-01-01T00:00:00Z


def generate_test_data(file_out):
    """ Generate test dataset. """
    tmp = array(list(product(range(-90, 91, 5), range(-180, 181, 10))))
    lat, lon = tmp[..., 0], tmp[..., 1]
    rad = uniform(0*EARTH_RADIUS, 2*EARTH_RADIUS, lat.shape)
    time_mjd2000 = uniform(START_TIME, END_TIME, lat.shape)
    decl, rasc, lha, azim, znth = sunpos(time_mjd2000, lat, lon, rad)

    header = [
        "MJD2000", "Latitude", "Longitude", "Radius",
        "Declination", "RightAscension", "HourAngle", "Azimuth", "Zenith",
    ]
    records = zip(time_mjd2000, lat, lon, rad, decl, rasc, lha, azim, znth)

    file_out.write("\t".join(header) + "\r\n")
    for record in records:
        file_out.write("\t".join(f"{value:.14e}" for value in record) + "\r\n")


if __name__ == "__main__":
    generate_test_data(sys.stdout)
