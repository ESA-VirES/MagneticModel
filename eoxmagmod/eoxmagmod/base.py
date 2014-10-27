#-------------------------------------------------------------------------------
#
#  Magnetic Model
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

import os.path
import numpy as np

# location of the data files
dirname = os.path.dirname(__file__)
DATA_WMM_2010 = os.path.join(dirname, 'data/WMM.COF')
DATA_EMM_2010_STATIC = os.path.join(dirname, 'data/EMM-720_V3p0_static.cof')
DATA_EMM_2010_SECVAR = os.path.join(dirname, 'data/EMM-720_V3p0_secvar.cof')

# coordinate systems and their trasnformation
from _pywmm import (
    GEODETIC_ABOVE_WGS84, GEODETIC_ABOVE_EGM96,
    GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
)

COORD_TYPES = (
    (GEODETIC_ABOVE_WGS84, "GEODETIC_ABOVE_WGS84"),
    (GEODETIC_ABOVE_EGM96, "GEODETIC_ABOVE_EGM96"),
    (GEOCENTRIC_SPHERICAL, "GEOCENTRIC_SPHERICAL"),
    (GEOCENTRIC_CARTESIAN, "GEOCENTRIC_CARTESIAN"),
)


class MagneticModel(object):
    """ Base Magnetic model class """

    # Supported coordinate systems:
    # Geodetic lat/lon coordinates with height above WGS84 ellipsoid. (lat, lon, height/elevation)
    GEODETIC_ABOVE_WGS84 = GEODETIC_ABOVE_WGS84
    # Geodetic lat/lon coordinates with height above EGM96 geoid. (lat, lon, height/elevation)
    GEODETIC_ABOVE_EGM96 = GEODETIC_ABOVE_EGM96
    # Geocentric spherical coordinates. (lat/azimut(-180,180), lon/elevation(-90,90), radius)
    GEOCENTRIC_SPHERICAL = GEOCENTRIC_SPHERICAL
    # Geocentric cartesian coordinates (x, y, z), x-axis points to prime-meridian/equator intersection
    GEOCENTRIC_CARTESIAN = GEOCENTRIC_CARTESIAN
    COORD_TYPES = COORD_TYPES

    def __init__(self, model_prm):
        """ Model constructor """
        if isinstance(model_prm, dict):
            self.prm = model_prm
        elif isinstance(model_prm, MagneticModel):
            self.prm = model_prm.prm

    @property
    def epoch(self):
        return self.prm['epoch']

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


    @property
    def coef_secvar(self):
        """Get secular variation coeficients."""
        return self.prm['coef_secvar_g'], self.prm['coef_secvar_h']

    @staticmethod
    def get_intensity(arr):
        """Calculate intensities for an array of vectors."""
        return np.sqrt((arr*arr).sum(axis=arr.ndim-1))

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
        print "\tepoch:    ", prm.get('epoch', 'n/a')
        print "\tsource(s):"
        for src in prm.get('sources', []):
            print "\t\t", src
        print "\tstatic model:"
        _analyse_coeficiens(prm['degree_static'], prm['coef_static_g'],
                            prm['coef_static_h'], "\t\t")
        print "\tsecvar model:"
        _analyse_coeficiens(prm['degree_secvar'], prm['coef_secvar_g'],
                            prm['coef_secvar_h'], "\t\t")

