#-------------------------------------------------------------------------------
#
#  EOX Magnetic Model Library
#
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

from .data import (
    # Data items are renamed to provide backward compatibility.
    WMM_2010 as DATA_WMM_2010,
    WMM_2015 as DATA_WMM_2015,
    EMM_2010_STATIC as DATA_EMM_2010_STATIC,
    EMM_2010_SECVAR as DATA_EMM_2010_SECVAR,
    CHAOS5_CORE as DATA_CHAOS5_CORE,
    CHAOS5_CORE_V4 as DATA_CHAOS5_CORE_V4,
    CHAOS5_STATIC as DATA_CHAOS5_STATIC,
    CHAOS6_CORE as DATA_CHAOS6_CORE,
    CHAOS6_CORE_X3 as DATA_CHAOS6_CORE_X3,
    CHAOS6_STATIC as DATA_CHAOS6_STATIC,
    IGRF11 as DATA_IGRF11,
    IGRF12 as DATA_IGRF12,
    SIFM as DATA_SIFM,
    APEX_2015 as DATA_APEX_2015,
    APEX_2020 as DATA_APEX_2020,
)
from .base import MagneticModel
from .util import (
    vnorm,
    vrotate,
    vincdecnorm,
)
from ._pywmm import (
    GEODETIC_ABOVE_WGS84,
    GEODETIC_ABOVE_EGM96,
    GEOCENTRIC_SPHERICAL,
    GEOCENTRIC_CARTESIAN,
    POTENTIAL,
    GRADIENT,
    POTENTIAL_AND_GRADIENT,
    convert,
    vrot_sph2geod,
    vrot_sph2cart,
    vrot_cart2sph,
    legendre,
    lonsincos,
    relradpow,
    spharpot,
    sphargrd,
    sheval,
)
from .emm import read_model_emm2010
from .wmm import read_model_wmm, read_model_wmm2010, read_model_wmm2015
from .shc import read_model_shc
from .igrf import read_model_igrf11
from .quasi_dipole_coordinates import (
    eval_qdlatlon, eval_mlt, eval_subsol,
    eval_qdlatlon_with_base_vectors,
)
from .solar_position import sunpos, sunpos_original
from .dipole_coords import (
    get_dipole_rotation_matrix, convert_to_dipole, vrot_from_dipole,
)
from .sheval_dipole import sheval_dipole

__all__ = [
    'MagneticModel',
    'read_model_wmm',
    'read_model_wmm2010',
    'read_model_wmm2015',
    'read_model_emm2010',
    'read_model_shc',
    'read_model_igrf11',
    'vnorm',
    'vincdecnorm',
    'vrotate',
    'vrot_sph2geod',
    'vrot_sph2cart',
    'vrot_cart2sph',
    'convert',
    'legendre',
    'lonsincos',
    'relradpow',
    'spharpot',
    'sphargrd',
    'sheval',
    'get_dipole_rotation_matrix',
    'vrot_from_dipole',
    'convert_to_dipole',
    'sheval_dipole',
    'DATA_WMM_2010',
    'DATA_WMM_2015',
    'DATA_EMM_2010_STATIC',
    'DATA_EMM_2010_SECVAR',
    'DATA_CHAOS5_CORE',
    'DATA_CHAOS5_CORE_V4',
    'DATA_CHAOS5_STATIC',
    'DATA_IGRF11',
    'DATA_IGRF12',
    'DATA_SIFM',
    'GEODETIC_ABOVE_WGS84',
    'GEODETIC_ABOVE_EGM96',
    'GEOCENTRIC_SPHERICAL',
    'GEOCENTRIC_CARTESIAN',
    'POTENTIAL',
    'GRADIENT',
    'POTENTIAL_AND_GRADIENT',
    'eval_qdlatlon',
    'eval_mlt',
    'eval_subsol',
    'DATA_APEX_2015',
    'DATA_APEX_2020',
    'sunpos',
    'sunpos_original',
]

__version__ = '0.5.0'
__author__ = 'Martin Paces (martin.paces@eox.at)'
__copyright__ = 'Copyright (C) 2014 EOX IT Services GmbH'
__licence__ = 'EOX licence (MIT style)'
