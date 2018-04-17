#-------------------------------------------------------------------------------
#
#  Swarm MMA_SHA_2C coefficients loader.
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

from spacepy import pycdf
from eoxmagmod.magnetic_model.parser_mma import (
    read_swarm_mma_2c_internal, read_swarm_mma_2c_external,
)
from .coefficients import (
    SparseSHCoefficientsTimeDependent, CombinedSHCoefficients,
)


def load_swarm_mma_2c_internal(path):
    """ Load internal model coefficients from a Swarm MMA_SHA_2C product file.
    """
    with pycdf.CDF(path) as cdf:
        data = read_swarm_mma_2c_internal(cdf)

    return CombinedSHCoefficients(*[
        SparseSHCoefficientsTimeDependent(
            item["nm"], item["gh"], item["t"], is_internal=True
        ) for item in data
    ])


def load_swarm_mma_2c_external(path):
    """ Load external model coefficients from a Swarm MMA_SHA_2C product file.
    """
    with pycdf.CDF(path) as cdf:
        data = read_swarm_mma_2c_external(cdf)

    return CombinedSHCoefficients(*[
        SparseSHCoefficientsTimeDependent(
            item["nm"], item["qs"], item["t"], is_internal=False
        ) for item in data
    ])
