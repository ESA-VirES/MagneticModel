#-------------------------------------------------------------------------------
#
#  EOX Magnetic Model Library
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2014-2024 EOX IT Services GmbH
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

try:
    from .util import (
        vnorm,
        vrotate,
        vincdecnorm,
    )
    from .time_util import (
        decimal_year_to_mjd2000,
        mjd2000_to_decimal_year,
        mjd2000_to_year_fraction,
    )
    from ._pymm import (
        EARTH_RADIUS,
        GEODETIC_ABOVE_WGS84,
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
        loncossin,
        relradpow,
        spharpot,
        sphargrd,
        sheval,
    )
    from .quasi_dipole_coordinates import (
        eval_qdlatlon, eval_mlt, eval_subsol,
        eval_qdlatlon_with_base_vectors,
        QDIPOLE_VERSION,
    )
    from .solar_position import sunpos, sunpos_original
    from .dipole_coords import (
        get_dipole_rotation_matrix, convert_to_dipole, vrot_from_dipole,
    )
    from .dipole import sheval_dipole
    from .magnetic_time import mjd2000_to_magnetic_universal_time
    from .magnetic_model.loader_shc import load_model_shc, load_model_shc_combined
    from .magnetic_model.loader_igrf import load_model_igrf
    from .magnetic_model.loader_wmm import load_model_wmm
    from .magnetic_model.loader_emm import load_model_emm
    from .magnetic_model.loader_mma import (
        load_model_swarm_mma_2c_internal,
        load_model_swarm_mma_2c_external,
        load_model_swarm_mma_2f_geo_internal,
        load_model_swarm_mma_2f_geo_external,
        load_model_swarm_mma_2f_sm_internal,
        load_model_swarm_mma_2f_sm_external,
    )
    from .magnetic_model.loader_mio import (
        load_model_swarm_mio_internal,
        load_model_swarm_mio_external,
    )
    from .magnetic_model.field_lines import trace_field_line
    from .magnetic_model.model_composed import ComposedGeomagneticModel

    __all__ = [
        'vnorm',
        'vincdecnorm',
        'vrotate',
        'vrot_sph2geod',
        'vrot_sph2cart',
        'vrot_cart2sph',
        'convert',
        'legendre',
        'loncossin',
        'relradpow',
        'spharpot',
        'sphargrd',
        'sheval',
        'get_dipole_rotation_matrix',
        'vrot_from_dipole',
        'convert_to_dipole',
        'sheval_dipole',
        'EARTH_RADIUS',
        'GEODETIC_ABOVE_WGS84',
        'GEOCENTRIC_SPHERICAL',
        'GEOCENTRIC_CARTESIAN',
        'POTENTIAL',
        'GRADIENT',
        'POTENTIAL_AND_GRADIENT',
        'eval_qdlatlon',
        'eval_qdlatlon_with_base_vectors',
        'eval_mlt',
        'eval_subsol',
        'QDIPOLE_VERSION',
        'sunpos',
        'sunpos_original',
        'decimal_year_to_mjd2000',
        'mjd2000_to_decimal_year',
        'mjd2000_to_year_fraction',
        'mjd2000_to_magnetic_universal_time',
        'load_model_shc',
        'load_model_shc_combined',
        'load_model_igrf',
        'load_model_wmm',
        'load_model_emm',
        'load_model_swarm_mma_2c_internal',
        'load_model_swarm_mma_2c_external',
        'load_model_swarm_mma_2f_geo_internal',
        'load_model_swarm_mma_2f_geo_external',
        'load_model_swarm_mma_2f_sm_internal',
        'load_model_swarm_mma_2f_sm_external',
        'load_model_swarm_mio_internal',
        'load_model_swarm_mio_external',
        'trace_field_line',
        'ComposedGeomagneticModel',
    ]
except ImportError:
    pass

__version__ = '0.14.0'
__author__ = 'Martin Paces (martin.paces@eox.at)'
__copyright__ = 'Copyright (C) 2014-2024 EOX IT Services GmbH'
__licence__ = 'EOX licence (MIT style)'
