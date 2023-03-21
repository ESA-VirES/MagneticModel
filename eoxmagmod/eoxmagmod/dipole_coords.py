#-------------------------------------------------------------------------------
#
#  Magnetic Dipole Coordinates
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

from math import pi, sin, cos
from numpy import array, dot
from eoxmagmod._pymm import (
    convert, vrot_sph2cart, vrot_cart2sph,
    GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
)

__all__ = [
    "get_dipole_rotation_matrix",
    "convert_to_dipole",
    "vrot_from_dipole",
]

DEG2RAD = pi / 180.0


def get_dipole_rotation_matrix(latitude, longitude):
    """ Get rotation matrix for given north pole coordinates.
    """
    sin_lat, cos_lat = sin(DEG2RAD*latitude), cos(DEG2RAD*latitude)
    sin_lon, cos_lon = sin(DEG2RAD*longitude), cos(DEG2RAD*longitude)
    return array([
        [sin_lat*cos_lon, -sin_lon, cos_lat*cos_lon],
        [sin_lat*sin_lon, cos_lon, cos_lat*sin_lon],
        [-cos_lat, 0, sin_lat],
    ])


def convert_to_dipole(coords, lat_ngp, lon_ngp,
                      coord_type_in=GEOCENTRIC_SPHERICAL):
    """ Convert coordinates (by default geocentric spherical)
    to dipole coordinates defined by the given latitude and longitude
    of the geomagnetic pole.
    The dipole coordinates are a simple rotated coordinate frame in which
    the North pole (dipole latitude == 0) is aligned with the geomagnetic pole
    and the prime meridian (dipole longitude == 0) is the meridian
    passing trough the geomagnetic pole.
    """
    rotation_matrix = get_dipole_rotation_matrix(lat_ngp, lon_ngp)
    coords = convert(coords, coord_type_in, GEOCENTRIC_CARTESIAN)
    coords = dot(coords, rotation_matrix)
    coords = convert(coords, GEOCENTRIC_CARTESIAN, GEOCENTRIC_SPHERICAL)
    return coords


def vrot_from_dipole(vectors, lat_ngp, lon_ngp, lat_dipole, lon_dipole,
                     lat_out=None, lon_out=None,
                     coord_type_out=GEOCENTRIC_CARTESIAN):
    """ Rotate vectors from dipole (NEC) coordinate frame to the
    Cartesian (XYZ) (default), geocentric spherical or geodetic (NEC)
    coordinate frame.
    """
    rotation_matrix = get_dipole_rotation_matrix(lat_ngp, lon_ngp).transpose()
    vectors = vrot_sph2cart(vectors, lat_dipole, lon_dipole)
    vectors = dot(vectors, rotation_matrix)
    if coord_type_out == GEOCENTRIC_CARTESIAN:
        # coordinates are already in the Cartesian coordinates
        return vectors
    return vrot_cart2sph(vectors, lat_out, lon_out)
