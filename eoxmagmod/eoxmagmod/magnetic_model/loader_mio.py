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
from .model_mio import DipoleMIOGeomagneticModel, MIO_EARTH_RADIUS
from .coefficients_mio import SparseSHCoefficientsMIO
from .parser_mio import parse_swarm_mio_file


def load_model_swarm_mio_internal(path, f107=0.0):
    """ Load internal (secondary field) model from a Swarm MIO_SHA_2* product.

    The loader requires F10.7 index source (see `DipoleMIOGeomagneticModel`
    for details).
    """
    coefficients, params = load_coeff_swarm_mio_internal(path)
    return _create_mio_model(coefficients, params, f107)


def load_model_swarm_mio_external(path, f107=0.0, above_ionosphere=True):
    """ Load external (primary field) model from a Swarm MIO_SHA_2* product.

    The loader requires F10.7 index source (see `DipoleMIOGeomagneticModel`
    for details).
    """
    coefficients, params = load_coeff_swarm_mio_external(path, above_ionosphere)
    return _create_mio_model(coefficients, params, f107)


def _create_mio_model(coefficients, params, f107):
    return DipoleMIOGeomagneticModel(
        coefficients, north_pole=(params["lat_NGP"], params["lon_NGP"]),
        f107=f107, wolf_ratio=params["wolf_ratio"], height=params["height"],
    )


def load_coeff_swarm_mio_internal(path):
    """ Load internal model coefficients and other parameters
    from a Swarm MIO_SHA_2* product file.
    """
    with open(path, "rb") as file_in:
        data = parse_swarm_mio_file(file_in)

    return SparseSHCoefficientsMIO(
        data["nm"], data["gh"],
        ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
        is_internal=True,
    ), data


def load_coeff_swarm_mio_external(path, above_ionosphere=True):
    """ Load external model coefficients from a Swarm MIO_SHA_2* product file.
    Use the `above_ionosphere` to pick the right model variant.
    """
    with open(path, "rb") as file_in:
        data = parse_swarm_mio_file(file_in)

    indices = data["nm"]
    coefficients = data["qs"]
    if above_ionosphere:
        is_internal = False
    else:
        is_internal = True
        coefficients = convert_external_mio_coeff(
            data["degree_max"], indices, coefficients, data["height"]
        )

    return SparseSHCoefficientsMIO(
        indices, coefficients,
        ps_extent=(data["pmin"], data["pmax"], data["smin"], data["smax"]),
        is_internal=is_internal,
    ), data


def convert_external_mio_coeff(degree, indices, coefficients, height):
    """ Convert external coefficients to internal ones. """
    nrrad = -(1.0 + height/MIO_EARTH_RADIUS)
    order = arange(degree + 1, dtype='float')
    scale = (order/(order + 1)) * nrrad**(2*order + 1)
    return (scale[indices[:, 0]] * coefficients.transpose()).transpose()
