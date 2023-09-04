#-------------------------------------------------------------------------------
#
#  SHC file format model loader
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2018-2022 EOX IT Services GmbH
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

from ..time_util import decimal_year_to_mjd2000
from .util import parse_file
from .model import SphericalHarmomicGeomagneticModel
from .coefficients import (
    SHCoefficients,
    SparseSHCoefficientsTimeDependent,
    SparseSHCoefficientsTimeDependentDecimalYear,
    SparseSHCoefficientsConstant,
    CombinedSHCoefficients,
    ComposedSHCoefficients,
)
from .parser_shc import parse_shc_file

__all__ = [
    "load_model_shc_combined",
    "load_model_shc",
    "load_coeff_shc_combined",
    "load_coeff_shc_composed",
    "load_coeff_shc",
]


def load_model_shc_combined(*paths, **kwargs):
    """ Load model with coefficients combined from multiple SHC files.
    E.g., combining core and lithospheric models.
    """
    return SphericalHarmomicGeomagneticModel(
        load_coeff_shc_combined(*paths, **kwargs)
    )


def load_model_shc(*paths, **kwargs):
    """ Load composed model from one or more SHC files.
    E.g., composing multiple core models, each with a different time extent.
    """
    return SphericalHarmomicGeomagneticModel(
        load_coeff_shc_composed(*paths, **kwargs)
    )


def load_coeff_shc_combined(*paths, **kwargs):
    """ Load coefficients combined from multiple SHC files.
    To be used for a combination of complementing models (core + lithosphere).
    """
    return CombinedSHCoefficients(*[
        _load_coeff_shc_or_coeff(path, **kwargs) for path in paths
    ])


def load_coeff_shc_composed(*paths, **kwargs):
    """ Load coefficients composed of multiple SHC files.
    To be used for a temporal sequence of models.
    """
    return ComposedSHCoefficients(*[
        _load_coeff_shc_or_coeff(path, **kwargs) for path in paths
    ])


def _load_coeff_shc_or_coeff(input_, **kwargs):
    if isinstance(input_, SHCoefficients):
        return input_
    return load_coeff_shc(input_, **kwargs)


def load_coeff_shc(path, interpolate_in_decimal_years=False, **kwargs):
    """ Load coefficients from an SHC file.

    The `interpolate_in_decimal_years` flag forces interpolation
    of time dependent models to be performed in the decimal years
    rather then in the default.
    """
    data = parse_file(parse_shc_file, path)


    keys = ["validity_start", "validity_end"]

    options = {
        key: data[key]
        for key in keys if key in data
    }
    options.update(kwargs) # extend or override the default model options

    if not "to_mjd2000" in options:
        options["to_mjd2000"] = decimal_year_to_mjd2000

    times = data["t"]

    kwargs["spline_oder"] = kwargs.get("spline_oder", min(data["spline_order"], 2))

    if len(times) == 1:
        return SparseSHCoefficientsConstant(
            data["nm"], data["gh"][:, 0], **options
        )

    if interpolate_in_decimal_years:
        coeff_class = SparseSHCoefficientsTimeDependentDecimalYear
    else:
        coeff_class = SparseSHCoefficientsTimeDependent

    return coeff_class(
        data["nm"], data["gh"], times, **options
    )
