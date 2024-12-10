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

import sys
from numpy import array
from eoxmagmod import mjd2000_to_decimal_year
from eoxmagmod.quasi_dipole_coordinates import (
    eval_mlt, eval_subsol, eval_qdlatlon_with_base_vectors,
)


def load_coordinates(file_in):

    _ = next(file_in)

    times, lats, lons, rads = [], [], [], []

    for line in file_in:
        time, _, lat, lon, rad = [float(v) for v in line.split()[:5]]
        times.append(time)
        lats.append(lat)
        lons.append(lon)
        rads.append(rad)

    return array(times), array(lats), array(lons), array(rads)


def update_test_data(file_out, time_mjd2000, lat, lon, rad):
    """ Generate test dataset. """
    time_decimal_year = mjd2000_to_decimal_year(time_mjd2000)

    qdlat, qdlon, f11, f12, f21, f22, f__ = eval_qdlatlon_with_base_vectors(
        lat, lon, rad, time_decimal_year
    )
    mlt = eval_mlt(qdlon, time_mjd2000)
    sol_lat, sol_lon = eval_subsol(time_mjd2000)

    header = [
        "MJD2000", "DecimalYear", "Latitude", "Longitude", "Radius",
        "QDLatitude", "QDLongitude", "F11", "F12", "F21", "F22", "F",
        "MagneticLocalTime", "SubsolarLatitude", "SubsolarLongitude",
    ]
    records = zip(
        time_mjd2000, time_decimal_year, lat, lon, rad,
        qdlat, qdlon, f11, f12, f21, f22, f__,
        mlt, sol_lat, sol_lon
    )

    file_out.write("\t".join(header) + "\n")
    for record in records:
        file_out.write("\t".join(f"{value:.14e}" for value in record) + "\n")


def main(filename):
    with open(filename, encoding="ascii") as file_in:
        time, lat, lon, rad = load_coordinates(file_in)
    update_test_data(sys.stdout, time, lat, lon, rad)


if __name__ == "__main__":
    main(sys.argv[1])
