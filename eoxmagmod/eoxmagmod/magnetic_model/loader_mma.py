#-------------------------------------------------------------------------------
#
#  Swarm MMA_SHA_2C model loaders.
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2018 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in all
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
# pylint: disable=invalid-name

from spacepy import pycdf
from .model import DipoleSphericalHarmomicGeomagneticModel
from .coefficients import (
    SparseSHCoefficientsTimeDependent, CombinedSHCoefficients,
)
from .parser_mma import (
    read_swarm_mma_2c_internal, read_swarm_mma_2c_external,
)

# North geomagnetic coordinates used by the MMA products (IGRF-11, 2010.0)
MMA2C_NGP_LATITUDE = 90 - 9.92   # deg.
MMA2C_NGP_LONGITUDE = 287.78 - 360.0  # deg.


def load_model_swarm_mma_2c_internal(path, lat_ngp=MMA2C_NGP_LATITUDE,
                                     lon_ngp=MMA2C_NGP_LONGITUDE):
    """ Load internal (secondary field) model from a Swarm MMA_SHA_2C product.
    """
    return DipoleSphericalHarmomicGeomagneticModel(
        load_coeff_swarm_mma_2c_internal(path),
        north_pole=(lat_ngp, lon_ngp),
    )


def load_model_swarm_mma_2c_external(path, lat_ngp=MMA2C_NGP_LATITUDE,
                                     lon_ngp=MMA2C_NGP_LONGITUDE):
    """ Load external (primary field) model from a Swarm MMA_SHA_2C product.
    """
    return DipoleSphericalHarmomicGeomagneticModel(
        load_coeff_swarm_mma_2c_external(path),
        north_pole=(lat_ngp, lon_ngp),
    )


def load_coeff_swarm_mma_2c_internal(path):
    """ Load internal model coefficients from a Swarm MMA_SHA_2C product file.
    """
    with pycdf.CDF(path) as cdf:
        data = read_swarm_mma_2c_internal(cdf)

    return CombinedSHCoefficients(*[
        SparseSHCoefficientsTimeDependent(
            item["nm"], item["gh"], item["t"], is_internal=True
        ) for item in data
    ])


def load_coeff_swarm_mma_2c_external(path):
    """ Load external model coefficients from a Swarm MMA_SHA_2C product file.
    """
    with pycdf.CDF(path) as cdf:
        data = read_swarm_mma_2c_external(cdf)

    return CombinedSHCoefficients(*[
        SparseSHCoefficientsTimeDependent(
            item["nm"], item["qs"], item["t"], is_internal=False
        ) for item in data
    ])
