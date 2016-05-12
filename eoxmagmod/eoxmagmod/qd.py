#-------------------------------------------------------------------------------
#
#  Quasi-Dipole Apex Coordinates / Geomagnetism Library
#
# Project: Earth magnetic field in Python.
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2015 EOX IT Services GmbH
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

import os.path
import _pyqd

# location of the data files
dirname = os.path.dirname(__file__)
dirname = os.path.join(dirname, 'data')
DATA_APEX_2015 = os.path.join(dirname, 'apexsh_1995-2015.txt')
DATA_APEX_2020 = os.path.join(dirname, 'apexsh_1980-2020.txt')

def eval_apex(gclat, gclon, gcrad, time, fname=DATA_APEX_2020):
    """
        qdlat, qdlon, mlt = eval_apex(gclat, gclon, gcrad, time, fname)

          Evaluate magnetic quasi-dipole coordinates and the magnetic
          local time for a single or multiple input coordinates.

          Inputs:
            gclat - geocentric latitude(s).
            gclon - geocentric longitude(s).
            gcrad - geocentric radial coordinate(s) in km.
            time  - decimal year time(s)
            fname - file-name of the model text file.

          Outputs:
            qdlat - quasi-dipole latitude(s).
            qdlon - quasi-dipole longitudes(s).
            mlt - magnetic local time(s).
    """

    if not os.path.isfile(fname):
        raise IOError("File not found! fname=%r" % fname)

    return _pyqd.eval_apex(gclat, gclon, gcrad, time, fname)[:3]
