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
from base import MagneticModel, DATA_WMM_2010, DATA_WMM_2015

class MagneticModelSimple(MagneticModel):
    """ Simple magnetic model class.
        To be used for models with a single epoch and secular variation.
    """

    def __init__(self, model_prm):
        """ Model constructor """
        super(MagneticModelSimple, self).__init__(model_prm)

    @property
    def epoch(self):
        return self.prm['epoch']

    @property
    def validity(self):
        return (self.prm['epoch'], self.prm['valid_until'])

    @property
    def degree_static(self):
        return self.prm['degree_static']

    @property
    def degree_secvar(self):
        return self.prm['degree_secvar']

    def get_coef_static(self, date):
        """ Calculate model static coeficients for a date specified by a decimal year value.
        """
        prm = self.prm
        ddate = date - self.epoch
        nterm = min(prm['coef_static_g'].size, prm['coef_secvar_g'].size)

        coef_static_g = np.copy(prm['coef_static_g'])
        coef_static_h = np.copy(prm['coef_static_h'])
        coef_static_g[:nterm] += ddate * prm['coef_secvar_g'][:nterm]
        coef_static_h[:nterm] += ddate * prm['coef_secvar_h'][:nterm]

        return coef_static_g, coef_static_h

    def get_coef_secvar(self, date):
        """Get secular variation coeficients."""
        return self.prm['coef_secvar_g'], self.prm['coef_secvar_h']

    def print_info(self):
        def _analyse_coeficiens(degree, coef_g, coef_h, prefix):
            nz_g = coef_g.nonzero()
            nz_h = coef_h.nonzero()
            nz_max = max(np.max(nz_g[0]), np.max(nz_h[0]))
            degree_real = 0
            for dg in xrange(degree, 0, -1):
                if nz_max >= ((dg*(dg+1))/2):
                    degree_real = dg
                    break
            n_all = coef_g.size + coef_h.size
            n_zero = n_all - (nz_g[0].size + nz_h[0].size)
            sparsity = float(n_zero) / float(n_all)

            print "%sdegree:      %d"%(prefix, degree)
            print "%strue degree: %d"%(prefix, degree_real)
            print "%ssparsity:    %g%% (%d of %d)"%(prefix, 100*sparsity, n_zero, n_all)

        prm = self.prm
        print "Magnetic Model:"
        print "\tclass:    ", self.__class__.__name__
        print "\tname:     ", prm.get('name', 'n/a')
        print "\tversion:  ", prm.get('version', 'n/a')
        print "\tepoch:    ", self.epoch
        print "\tvalidity: ", self.validity
        print "\tsource(s):"
        for src in prm.get('sources', []):
            print "\t\t", src
        print "\tstatic model:"
        _analyse_coeficiens(prm['degree_static'], prm['coef_static_g'],
                            prm['coef_static_h'], "\t\t")
        print "\tsecvar model:"
        _analyse_coeficiens(prm['degree_secvar'], prm['coef_secvar_g'],
                            prm['coef_secvar_h'], "\t\t")



def read_model_wmm2010():
    """ Read WMM2015 model coefficients."""
    return read_model_wmm(DATA_WMM_2010)


def read_model_wmm2015():
    """ Read WMM2015 model coefficients."""
    return read_model_wmm(DATA_WMM_2010)


def read_model_wmm(fname):
    """ Read model coefficients from a WMM coefficient file. """
    prm = {'sources': [fname]}

    with file(fname, 'r') as fid:

        line = next(fid).rstrip()
        epoch, prm['name'], prm['version'] = line.split()
        prm['epoch'] = float(epoch)
        prm['valid_until'] = prm['epoch'] + 5
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

    return MagneticModelSimple(prm)
