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
from numpy import concatenate, array_equal
from .model import (
    SphericalHarmomicGeomagneticModel,
    DipoleSphericalHarmomicGeomagneticModel,
)
from .coefficients import (
    SparseSHCoefficientsTimeDependent, CombinedSHCoefficients,
)
from .parser_mma import (
    read_swarm_mma_2c_internal, read_swarm_mma_2c_external,
    read_swarm_mma_2f_geo_internal, read_swarm_mma_2f_geo_external,
    read_swarm_mma_2f_sm_internal, read_swarm_mma_2f_sm_external,
)

__all__ = [
    "load_model_swarm_mma_2c_internal",
    "load_model_swarm_mma_2c_external",
    "load_model_swarm_mma_2f_geo_internal",
    "load_model_swarm_mma_2f_geo_external",
    "load_model_swarm_mma_2f_sm_internal",
    "load_model_swarm_mma_2f_sm_external",
    "load_coeff_swarm_mma_2c_internal",
    "load_coeff_swarm_mma_2c_external",
    "load_coeff_swarm_mma_2f_geo_internal",
    "load_coeff_swarm_mma_2f_geo_external",
    "load_coeff_swarm_mma_2f_sm_internal",
    "load_coeff_swarm_mma_2f_sm_external",
]

# North geomagnetic coordinates used by the MMA products (IGRF-11, 2010.0)
MMA2C_NGP_LATITUDE = 90 - 9.92   # deg.
MMA2C_NGP_LONGITUDE = 287.78 - 360.0  # deg.


def load_model_swarm_mma_2c_internal(*paths, lat_ngp=MMA2C_NGP_LATITUDE,
                                     lon_ngp=MMA2C_NGP_LONGITUDE):
    """ Load internal (secondary field) model from a Swarm MMA_SHA_2C product.
    """
    return DipoleSphericalHarmomicGeomagneticModel(
        load_coeff_swarm_mma_2c_internal(*paths),
        north_pole=(lat_ngp, lon_ngp),
    )


def load_model_swarm_mma_2c_external(*paths, lat_ngp=MMA2C_NGP_LATITUDE,
                                     lon_ngp=MMA2C_NGP_LONGITUDE):
    """ Load external (primary field) model from a Swarm MMA_SHA_2C product.
    """
    return DipoleSphericalHarmomicGeomagneticModel(
        load_coeff_swarm_mma_2c_external(*paths),
        north_pole=(lat_ngp, lon_ngp),
    )


def load_model_swarm_mma_2f_geo_internal(*paths):
    """ Load geographic frame internal (secondary field) model from a Swarm
    MMA_SHA_2F product.
    """
    return SphericalHarmomicGeomagneticModel(
        load_coeff_swarm_mma_2f_geo_internal(*paths)
    )


def load_model_swarm_mma_2f_geo_external(*paths):
    """ Load geographic frame internal (primary field) model from a Swarm
    MMA_SHA_2F product.
    """
    return SphericalHarmomicGeomagneticModel(
        load_coeff_swarm_mma_2f_geo_external(*paths)
    )


def load_model_swarm_mma_2f_sm_internal(*paths):
    """ Load solar magnetic frame internal (secondary field) model from a Swarm
    MMA_SHA_2F product.
    """
    # FIXME: solar-magnetic frame model
    return SphericalHarmomicGeomagneticModel(
        load_coeff_swarm_mma_2f_sm_internal(*paths),
    )


def load_model_swarm_mma_2f_sm_external(*paths):
    """ Load solar magnetic frame internal (primary field) model from a Swarm
    MMA_SHA_2F product.
    """
    # FIXME: solar-magnetic frame model
    return SphericalHarmomicGeomagneticModel(
        load_coeff_swarm_mma_2f_sm_external(*paths)
    )


def load_coeff_swarm_mma_2c_internal(*paths):
    """ Load internal model coefficients from one or more Swarm MMA_SHA_2C product files.
    Note that the models are loaded and merged in the order of the arguments.
    """
    return _load_coeff_mma_multi_set(
        paths, read_swarm_mma_2c_internal, "gh", is_internal=True
    )


def load_coeff_swarm_mma_2c_external(*paths):
    """ Load external model coefficients from one or more Swarm MMA_SHA_2C product files.
    Note that the models are loaded and merged in the order of the arguments.
    """
    return _load_coeff_mma_multi_set(
        paths, read_swarm_mma_2c_external, "qs", is_internal=False
    )


def load_coeff_swarm_mma_2f_geo_internal(*paths):
    """ Load internal geographic frame model coefficients from a Swarm
    MMA_SHA_2F product file.
    """
    return _load_coeff_mma_single_set(
        paths, read_swarm_mma_2f_geo_internal, "gh", is_internal=True
    )


def load_coeff_swarm_mma_2f_geo_external(*paths):
    """ Load external geographic frame model coefficients from a Swarm
    MMA_SHA_2F product file.
    """
    return _load_coeff_mma_single_set(
        paths, read_swarm_mma_2f_geo_external, "qs", is_internal=False
    )


def load_coeff_swarm_mma_2f_sm_internal(*paths):
    """ Load internal solar magnetic frame model coefficients from a Swarm
    MMA_SHA_2F product file.
    """
    return _load_coeff_mma_single_set(
        paths, read_swarm_mma_2f_sm_internal, "gh", is_internal=True
    )


def load_coeff_swarm_mma_2f_sm_external(*paths):
    """ Load external solar magnetic frame model coefficients from a Swarm
    MMA_SHA_2F product file.
    """
    return _load_coeff_mma_single_set(
        paths, read_swarm_mma_2f_sm_external, "qs", is_internal=False
    )


def _load_coeff_mma_multi_set(paths, cdf_reader, variable, is_internal):
    data = [_read_data(path, cdf_reader) for path in paths]
    data = [
        _merge_coefficients(*items, variable=variable)
        for items in zip(*data)
    ]
    return CombinedSHCoefficients(*[
        SparseSHCoefficientsTimeDependent(
            item["nm"], item[variable], item["t"], is_internal=is_internal
        ) for item in data
    ])


def _load_coeff_mma_single_set(paths, cdf_reader, variable, is_internal):
    data = [_read_data(path, cdf_reader) for path in paths]
    data = _merge_coefficients(*data, variable=variable)
    return SparseSHCoefficientsTimeDependent(
        data["nm"], data[variable], data["t"], is_internal=is_internal
    )


def _read_data(path, cdf_reader):
    with pycdf.CDF(path) as cdf:
        return cdf_reader(cdf)


def _merge_coefficients(*sets, variable):
    head, *tail = sets

    # sanity check
    previous_item = head
    for item in tail:
        if (
            previous_item["degree_min"] != item["degree_min"] or
            previous_item["degree_max"] != item["degree_max"]
        ):
            raise ValueError("Degree range mismatch!")
        if not array_equal(previous_item["nm"], item["nm"]):
            raise ValueError("Incompatible sets of coefficients!")
        if previous_item["t"][-1] > item["t"][0]:
            raise ValueError("Unordered sets of coefficients!")
        previous_item = item

    return {
        "degree_min": head["degree_min"],
        "degree_max": head["degree_max"],
        "nm": head["nm"],
        "t": concatenate([item["t"] for item in sets], axis=0),
        variable: concatenate([item[variable] for item in sets], axis=1),
    }
