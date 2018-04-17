#-------------------------------------------------------------------------------
#
#  Swarm MIO_SHA_2* coefficients loaders
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

from numpy import arange
from eoxmagmod.magnetic_model.parser_mio import parse_swarm_mio_file
from .coefficients_mio import SparseSHCoefficientsMIO
from .parser_mio import parse_swarm_mio_file

MIO_EARTH_RADIUS = 6371.2 # km


def load_swarm_mio_internal(path):
    """ Load internal model coefficients and other parameters
    from a Swarm MIO_SHA_2* product file.
    """
    with open(path, "rb") as file_in:
        data = parse_swarm_mio_file(file_in)

    return SparseSHCoefficientsMIO(
        data["nm"], data["gh"],
        ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
        lat_ngp=data["lat_NGP"],
        lon_ngp=data["lon_NGP"],
        mio_radius=(data["height"] + MIO_EARTH_RADIUS),
        wolf_ratio=data["wolf_ratio"],
        is_internal=True,
    )


def load_swarm_mio_external(path, above_ionosphere=True):
    """ Load external model coefficients from a Swarm MIO_SHA_2* product file.
    Use the `above_ionosphere` to pick the right model variant.
    """
    with open(path, "rb") as file_in:
        data = parse_swarm_mio_file(file_in)

    indices = data["nm"]
    coefficients = data["qs"]
    mio_radius = data["height"] + MIO_EARTH_RADIUS
    if above_ionosphere:
        is_internal = False
    else:
        is_internal = True
        coefficients = convert_external_mio_coeff(
            data["degree_max"], indices, coefficients, mio_radius
        )

    return SparseSHCoefficientsMIO(
        indices, coefficients,
        ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
        lat_ngp=data["lat_NGP"],
        lon_ngp=data["lon_NGP"],
        mio_radius=mio_radius,
        wolf_ratio=data["wolf_ratio"],
        is_internal=is_internal,
    )


def convert_external_mio_coeff(degree, indices, coefficients, mio_radius):
    """ Convert external coefficients to internal ones. """
    nrrad = -mio_radius/MIO_EARTH_RADIUS
    order = arange(degree + 1, dtype='float')
    scale = (order/(order + 1)) * nrrad**(2*order + 1)
    return (scale[indices[:, 0]] * coefficients.transpose()).transpose()
