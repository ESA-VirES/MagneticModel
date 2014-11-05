#-------------------------------------------------------------------------------
#
#  World Magnetic Model 2010 / Geomagnetism Library
#  - IGRF format reader
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

import numpy as np
from shc import MagneticModelSHCPP

def read_model_igrf11(fname):
    """ Read model parameters from a coeficient file in the IGRF11 format."""

    prm = {'sources': [fname], 'headers': []}

    with file(fname, 'r') as fid:

        # parse text headers (comments)
        lidx = 0
        for line in fid:
            lidx += 1
            if line.strip().startswith("#"):
                prm['headers'].append(line)
                if lidx == 1:
                    line = line.strip()
                    prm['name'] = line.strip()[1:].lstrip()
            else:
                break

        # parse model headers
        prm['labels'] = line.split()

        # parse time labels
        line = next(fid)
        lidx += 1
        item = line.split()
        time = [float(v) for v in item[3:-1]]
        time.append(float(item[-1][:2]+item[-1][-2:]))
        time = np.array(time)

        # parse coeficients
        degree = 0
        lcoef = []
        for line in fid:
            lidx += 1
            item = line.split()
            item = (item[0], int(item[1]), int(item[2]),
                        np.array([float(v) for v in item[3:]]))
            lcoef.append(item)
            degree = max(item[1], degree)
            if item[0] not in ('g', 'h'):
               raise ValueError("Unexpected row label '%s' at line %d!"%(item[0], lidx))

        # fill the coeficient arrays
        ntime = time.size
        nterm = ((degree+1)*(degree+2))/2
        coef_g = np.zeros((nterm, ntime))
        coef_h = np.zeros((nterm, ntime))

        for label, i, j, coef in lcoef:
            idx = j + (i*(i+1))/2
            if label == 'g':
                coef_g[idx,:] = coef
            elif label == 'h':
                coef_h[idx,:] = coef

        # fix the last column
        dt = time[-1] - time[-2]
        coef_g[:,-1] = coef_g[:,-2] + dt*coef_g[:,-1]
        coef_h[:,-1] = coef_h[:,-2] + dt*coef_h[:,-1]

        prm.update({
            'degree_min': 0,
            'degree': degree,
            'spline_order': 1,
            'nstep': 0,
            'time': time,
            'coef_h': coef_h,
            'coef_g': coef_g,
        })

        return MagneticModelSHCPP(prm)
