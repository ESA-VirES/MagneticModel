#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Coefficients specific to Swarm MIO model.
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

from math import pi
from collections import namedtuple
from numpy import asarray, zeros
from .coefficients import SparseSHCoefficients, coeff_size
from ..time_util import (
    mjd2000_to_year_fraction as mjd2000_to_year_fraction_default,
)
from .._pymm import fourier2d

__all__ = [
    "SparseSHCoefficientsMIO",
]

SCALE_SEASONAL = 2*pi
SCALE_DIURNAL = 2*pi/24.

DegreeRanges = namedtuple("DegreeRanges", ["pmin", "pmax", "smin", "smax"])


class SparseSHCoefficientsMIO(SparseSHCoefficients):
    """ Time dependent 2D Fourier series Swarm MIO coefficients.
        Parameters:
            indices - array if the nm indices
            coefficients - gh or qs coefficients
            ps_extent - (pmin, pmax, smin, smax) tuple
            is_internal - set False for an external model
    """

    def get_f2fs_coeff_set(self, **parameters):
        """ Return coefficient set which can be passed to sheval2dfs. """
        _, coeff, nm_, _ = self.subset_degree(
            parameters.get("min_degree", -1),
            parameters.get("max_degree", -1),
        )
        return (
            coeff, nm_,
            self.ps_extent.smin, self.ps_extent.pmin,
            SCALE_SEASONAL, SCALE_DIURNAL
        )

    def __init__(self, indices, coefficients, ps_extent,
                 mjd2000_to_year_fraction=mjd2000_to_year_fraction_default,
                 **kwargs):
        SparseSHCoefficients.__init__(self, indices, coefficients, **kwargs)
        pmin, pmax, smin, smax = ps_extent
        if pmin > pmax or smin > smax:
            raise Exception(f"Invalid ps_extent {ps_extent}!")
        self.ps_extent = DegreeRanges(pmin, pmax, smin, smax)
        self.mjd2000_to_year_fraction = mjd2000_to_year_fraction

    def __call__(self, time, mut, **parameters):
        time = asarray(time)
        mut = asarray(mut)
        degree, coeff, _, index = self.subset_degree(
            parameters.get("min_degree", -1), parameters.get("max_degree", -1)
        )
        coeff_full = zeros((*time.shape, coeff_size(degree), 2))
        coeff_full[..., index[:, 0], index[:, 1]] = fourier2d(
            self.mjd2000_to_year_fraction(time), mut,
            coeff, self.ps_extent.smin, self.ps_extent.pmin,
            SCALE_SEASONAL, SCALE_DIURNAL
        )
        return coeff_full, degree
