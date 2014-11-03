#-------------------------------------------------------------------------------
#
#  World Magnetic Model 2010 / Geomagnetism Library
#
# Project: Earth magnetic field in Python.
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2014 EOX IT Services GmbH
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

import re
import numpy as np
from base import MagneticModel
import _pywmm
from _pywmm import (
    convert, legendre, lonsincos, relradpow, spharpot, sphargrd, geomag,
    vrot_sph2geod, vrot_sph2cart,
)

class GeomagWMM2010(MagneticModel):
    """ Magnetic Model based on the WMM2010 Geomagnetism library."""

    def eval(self, arr_in, date, coord_type_in=MagneticModel.GEODETIC_ABOVE_WGS84,
                coord_type_out=None, secvar=False, mode=MagneticModel.GRADIENT):
        """Evaluate Magnetic Model for a given set of spatio-teporal
        coordinates.

        Input:
            arr_in - numpy array of (x, y, z) coordinates.
                     The type of the x, y, z, values depends on the selected
                     input coordinate system.

            date - The time as a decimal year value (e.g., 2015.23).
            coord_type_in - shall be set to one of the valid coordinate systems:
                        MagneticModel.GEODETIC_ABOVE_WGS84 (default)
                        MagneticModel.GEODETIC_ABOVE_EGM96
                        MagneticModel.GEOCENTRIC_SPHERICAL
                        MagneticModel.GEOCENTRIC_CARTESIAN

            coord_type_out - coordinate system of the output vector
                        MagneticModel.GEODETIC_ABOVE_WGS84
                        MagneticModel.GEODETIC_ABOVE_EGM96 (the same as WGS84)
                        MagneticModel.GEOCENTRIC_SPHERICAL
                        MagneticModel.GEOCENTRIC_CARTESIAN
                    Ouput coordinate system defaults to the input coordinate
                    system.

            secvar - if True secular variation of the magentic field is
                     calculated rather than the magnetic fied.

            mode - type of output to be produced. The possible values are:
                        MagneticModel.POTENTIAL
                        MagneticModel.GRADIENT (default)
                        MagneticModel.POTENTIAL_AND_GRADIENT

        Output:

            arr_out - output numpy array with the same shape as the input
                      array contaning the calcuated magentic field parameters.
        """
        if mode not in dict(self.EVAL_MODES):
            raise ValueError("Invalid mode value!")

        if coord_type_in not in dict(self.COORD_TYPES):
            raise ValueError("Invalid input coordinate type!")

        if coord_type_out is None:
            coord_type_out = coord_type_in

        if coord_type_out not in dict(self.COORD_TYPES):
            raise ValueError("Invalid output coordinate type!")

        if secvar:
            degree = self.degree_secvar
            coef_g, coef_h = self.coef_secvar
        else:
            degree = self.degree_static
            coef_g, coef_h = self.get_coef_static(date)

        return geomag(arr_in, degree, coef_g, coef_h, coord_type_in, coord_type_out, mode)


def read_model_wmm2010(fname):
    """ Read model parameters from a coeficient file in the WMM2010 format."""
    prm = {'sources': [fname]}

    with file(fname, 'r') as fid:

        line = next(fid).rstrip()
        epoch, prm['name'], prm['version'] = line.split()
        prm['epoch'] = float(epoch)
        prm['headers'] = [line]

        # PASS1 - parse the input file
        degree = 0
        lcoef = []
        for lidx, line in enumerate(fid):
            try:
                n, m, g, h, dg, dh = line.split()
                n, m = int(n), int(m)
                g, h, dg, dh = float(g), float(h), float(dg), float(dh)
                if m > n:
                    raise ValueError
            except ValueError:
                if not re.match(r'^\s*9+\s*$', line):
                    raise ValueError("Invalid line #%d: %s"%(lidx+1, line.rstrip()))
            else:
                lcoef.append((n, m, g, h, dg, dh))
                degree = max(degree, n)

        # PASS2 - fill the coeficient arrays
        nterm = ((degree+2)*(degree+1))/2
        coef_g = np.zeros(nterm)
        coef_h = np.zeros(nterm)
        coef_dg = np.zeros(nterm)
        coef_dh = np.zeros(nterm)

        for n, m, g, h, dg, dh in lcoef:
            idx = m + ((n+1)*n)/2
            coef_g[idx] = g
            coef_h[idx] = h
            coef_dg[idx] = dg
            coef_dh[idx] = dh

        prm['degree_static'] = degree
        prm['degree_secvar'] = degree
        prm['coef_static_g'] = coef_g
        prm['coef_static_h'] = coef_h
        prm['coef_secvar_g'] = coef_dg
        prm['coef_secvar_h'] = coef_dh

    return MagneticModel(prm)
