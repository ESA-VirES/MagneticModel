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

from .base import (
    MagneticModel,
    DATA_WMM_2010,
    DATA_EMM_2010_STATIC,
    DATA_EMM_2010_SECVAR,
    DATA_CHAOS5_CORE,
    DATA_CHAOS5_STATIC,
    DATA_IGRF11,
    GEODETIC_ABOVE_WGS84,
    GEODETIC_ABOVE_EGM96,
    GEOCENTRIC_SPHERICAL,
    GEOCENTRIC_CARTESIAN,
    POTENTIAL,
    GRADIENT,
    POTENTIAL_AND_GRADIENT,
    vnorm,
    vincdecnorm,
    convert,
    legendre,
    lonsincos,
    relradpow,
    spharpot,
    sphargrd,
    vrot_sph2geod,
    vrot_sph2cart,
)
from .emm import read_model_emm2010
from .wmm import read_model_wmm2010
from .shc import read_model_shc
from .igrf import read_model_igrf11
from .qd import eval_apex, DATA_APEX_2015

__all__ = [
    'MagneticModel',
    'read_model_wmm2010',
    'read_model_emm2010',
    'read_model_shc',
    'read_model_igrf11',
    'vnorm',
    'vincdecnorm',
    'vrot_sph2geod',
    'vrot_sph2cart',
    'convert',
    'legendre',
    'lonsincos',
    'relradpow',
    'spharpot',
    'sphargrd',
    'DATA_WMM_2010',
    'DATA_EMM_2010_STATIC',
    'DATA_EMM_2010_SECVAR',
    'DATA_CHAOS5_CORE',
    'DATA_CHAOS5_STATIC',
    'DATA_IGRF11',
    'GEODETIC_ABOVE_WGS84',
    'GEODETIC_ABOVE_EGM96',
    'GEOCENTRIC_SPHERICAL',
    'GEOCENTRIC_CARTESIAN',
    'POTENTIAL',
    'GRADIENT',
    'POTENTIAL_AND_GRADIENT',
    'eval_apex',
    'DATA_APEX_2015',
]

__version__ = '0.2.0dev'
__author__ = 'Martin Paces (martin.paces@eox.at)'
__copyright__ = 'Copyright (C) 2014 EOX IT Services GmbH'
__licence__ = 'EOX licence (MIT style)'
