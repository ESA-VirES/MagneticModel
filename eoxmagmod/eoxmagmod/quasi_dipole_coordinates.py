#-------------------------------------------------------------------------------
#
#  Quasi-Dipole Apex Coordinates / Geomagnetism Library
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2015-2024 EOX IT Services GmbH
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

from os.path import isfile
from . import _pyqd
from .data import APEX_LATEST

__all__ = [
    "eval_qdlatlon",
    "eval_qdlatlon_with_base_vectors",
    "eval_mlt",
    "eval_subsol",
    "QDIPOLE_VERSION",
]

QDIPOLE_VERSION = _pyqd.QDIPOLE_VERSION

def eval_qdlatlon(gclat, gclon, gcrad, time, fname=APEX_LATEST):
    """
          Evaluate magnetic quasi-dipole coordinates a single or multiple input
          points.

          Inputs:
            gclat - geocentric latitude(s).
            gclon - geocentric longitude(s).
            gcrad - geocentric radial coordinate(s) in km.
            time  - decimal year time(s)
            fname - file-name of the model text file.

          Outputs:
            qdlat - quasi-dipole latitude(s).
            qdlon - quasi-dipole longitude(s).
    """
    if not isfile(fname):
        raise IOError(f"File not found! fname={fname!r}")
    return _pyqd.eval_qdlatlon(gclat, gclon, gcrad, time, fname)


def eval_qdlatlon_with_base_vectors(gclat, gclon, gcrad, time, fname=APEX_LATEST):
    """
          Evaluate magnetic quasi-dipole coordinates a single or multiple input
          coordinates.

          Inputs:
            gclat - geocentric latitude(s).
            gclon - geocentric longitude(s).
            gcrad - geocentric radial coordinate(s) in km.
            time  - decimal year time(s)
            fname - file-name of the model text file.

          Outputs:
            qdlat - quasi-dipole latitude(s).
            qdlon - quasi-dipole longitude(s).
            f11 - base vector F1 component 1
            f12 - base vector F1 component 2
            f21 - base vector F2 component 1
            f22 - base vector F2 component 2
            f - | F1 x F2 | value
    """
    if not isfile(fname):
        raise IOError(f"File not found! fname={fname!r}")
    return _pyqd.eval_qdlatlon(gclat, gclon, gcrad, time, fname, True)


def eval_mlt(qdlon, time):
    """
          Evaluate magnetic local time for given quasi dipole longitudes.

          Inputs:
            qdlon - quasi-dipole longitudes(s).
            time  - MJD2000 time(s)

          Outputs:
            mlt - magnetic local time(s).
    """
    return _pyqd.eval_mlt(qdlon, time)


def eval_subsol(time):
    """
          Evaluate sub-solar point coordinates.

          Inputs:
            time  - MJD2000 time(s)

          Outputs:
            gdlat - sub-solar point latitude(s).
            gdlon - sub-solar point longitude(s).
    """
    return _pyqd.eval_subsol(time)
