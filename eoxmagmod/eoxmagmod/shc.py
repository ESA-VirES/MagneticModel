#-------------------------------------------------------------------------------
#
#  World Magnetic Model 2010 / Geomagnetism Library
#  - SHC format reader
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
from base import MagneticModel, DATA_CHAOS5_CORE, DATA_CHAOS5_STATIC

class MagneticModelSHCPP(MagneticModel):
    """ SHC piecewise polynomial magnetic model class.
    """

    def __init__(self, model_prm):
        """ Model constructor """
        super(MagneticModelSHCPP, self).__init__(model_prm)

    @property
    def validity(self):
        return (np.min(self.prm['time']), np.max(self.prm['time']))

    @property
    def degree_static(self):
        return self.prm['degree']

    @property
    def degree_secvar(self):
        return self.prm['degree']

    def get_coef_secvar(self, date):
        """Get secular variation coeficients."""
        dates = self.prm['time']
        degree = self.prm['degree']
        degree_min = self.prm['degree_min']
        src_coef_g = self.prm['coef_h']
        src_coef_h = self.prm['coef_g']
        nterms = ((degree+2)*(degree+1))/2
        idxoff = (degree_min*(degree_min+1))/2
        # TODO: proper spline interpolation

        coef_g = np.zeros(nterms)
        coef_h = np.zeros(nterms)

        if dates.size > 1:
            # lookup the interval
            try:
                idx = (dates <= date).nonzero()[0].max()
            except ValueError:
                idx = 0

            a = 1.0/(dates[idx+1]-dates[idx])
            coef_g[idxoff:] = a*(src_coef_g[:,idx+1] - src_coef_g[:,idx])
            coef_h[idxoff:] = a*(src_coef_h[:,idx+1] - src_coef_h[:,idx])

        return coef_g, coef_h

    def get_coef_static(self, date):
        """ Calculate model static coeficients for a date specified by a decimal year value.
        """
        dates = self.prm['time']
        degree = self.prm['degree']
        degree_min = self.prm['degree_min']
        src_coef_g = self.prm['coef_g']
        src_coef_h = self.prm['coef_h']
        nterms = ((degree+2)*(degree+1))/2
        idxoff = (degree_min*(degree_min+1))/2
        # TODO: proper spline interpolation

        coef_g = np.zeros(nterms)
        coef_h = np.zeros(nterms)

        if dates.size == 1:
            coef_g[idxoff:] = src_coef_g[:,0]
            coef_h[idxoff:] = src_coef_h[:,0]

        elif dates.size > 1:
            # lookup the interval
            try:
                idx = (dates <= date).nonzero()[0].max()
            except ValueError:
                idx = 0

            a1 = (date-dates[idx])/(dates[idx+1]-dates[idx])
            a0 = 1.0 - a1

            coef_g[idxoff:] = a0*src_coef_g[:,idx] + a1*src_coef_g[:,idx+1]
            coef_h[idxoff:] = a0*src_coef_h[:,idx] + a1*src_coef_h[:,idx+1]

        return coef_g, coef_h

    def print_info(self):
        prm = self.prm
        print "Magnetic Model:"
        print "\tclass:    ", self.__class__.__name__
        print "\tname:     ", prm.get('name', 'n/a')
        print "\tversion:  ", prm.get('version', 'n/a')
        print "\tvalidity: ", self.validity
        print "\tdegree ", prm['degree']
        print "\tsource(s):"
        for src in prm.get('sources', []):
            print "\t\t", src


def read_model_shc(fname=DATA_CHAOS5_CORE):
    """ Read model parameters from a coeficient file in the SHC format."""

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
        header = [int(v) for v in line.split()]
        degree_min, degree, ntime, spline_order, nstep = header

        # parse time header
        line = next(fid)
        lidx += 1
        time = np.array([float(v) for v in line.split()])
        if time.size != ntime:
            ValueError("The list of times does not match the file header!")

        nterm = ((degree+1)*(degree+2)-(degree_min)*(degree_min+1))/2

        coef_g = np.zeros((nterm, ntime))
        coef_h = np.zeros((nterm, ntime))

        for line in fid:
            lidx += 1
            line = line.split()
            i, j = [int(v) for v in line[:2]]
            idx = abs(j)+((i)*(i+1)-(degree_min)*(degree_min+1))/2
            coef = np.array([float(v) for v in line[2:]])
            if j < 0:
                coef_h[idx,:] = coef
            else:
                coef_g[idx,:] = coef

        prm.update({
            'degree_min': degree_min,
            'degree': degree,
            'spline_order': spline_order,
            'nstep': nstep,
            'time': time,
            'coef_h': coef_h,
            'coef_g': coef_g,
        })

        return MagneticModelSHCPP(prm)
