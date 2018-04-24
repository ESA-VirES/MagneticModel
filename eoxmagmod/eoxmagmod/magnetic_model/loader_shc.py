#-------------------------------------------------------------------------------
#
#  SHC file format model loader
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

from .util import parse_file
from .model import SphericalHarmomicGeomagneticModel
from .coefficients import (
    SparseSHCoefficientsTimeDependentDecimalYear,
    SparseSHCoefficientsConstant,
    CombinedSHCoefficients,
)
from .parser_shc import parse_shc_file


def load_model_shc_combined(*paths):
    """ Load model with coefficients combined from multiple SHC files. """
    return SphericalHarmomicGeomagneticModel(load_coeff_shc_combined(*paths))


def load_model_shc(path):
    """ Load model from an SHC file. """
    return SphericalHarmomicGeomagneticModel(load_coeff_shc(path))


def load_coeff_shc_combined(*paths):
    """ Load coefficients combined from multiple SHC files. """
    return CombinedSHCoefficients(*[load_coeff_shc(path) for path in paths])


def load_coeff_shc(path):
    """ Load coefficients from an SHC file. """
    data = parse_file(parse_shc_file, path)

    options = {
        key: data[key]
        for key in ("validity_start", "validity_end") if key in data
    }

    times = data["t"]
    if len(times) == 1:
        return SparseSHCoefficientsConstant(
            data["nm"], data["gh"][:, 0], **options
        )
    else:
        return SparseSHCoefficientsTimeDependentDecimalYear(
            data["nm"], data["gh"], times, **options
        )
