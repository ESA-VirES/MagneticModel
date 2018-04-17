#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Coefficients.
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
# pylint: disable=too-few-public-methods,abstract-method

from numpy import (
    inf, array, zeros, dot, digitize, argsort, abs as aabs,
)


class SHCoefficients(object):
    """ Abstract base class for the spherical harmonic coefficients. """

    def __init__(self, is_internal=True, **kwargs):
        self.is_internal = is_internal
        self.validity = (
            kwargs.get("validity_start", -inf),
            kwargs.get("validity_end", +inf)
        )
        self.scale = kwargs.get("scale", 1.0)

    @property
    def degree(self):
        """ Get model degree. """
        raise NotImplementedError

    def __call__(self, time, **parameters):
        """ Return the matrix of the full model coefficients. """
        raise NotImplementedError


class CombinedSHCoefficients(SHCoefficients):
    """ Model composed of multiple coefficient sets. """

    def __init__(self, *items):
        if len(items) < 1:
            raise ValueError(
                "The composed model must be composed from at least "
                "coefficient set."
            )

        item = items[0]
        is_internal = item.is_internal
        validity_start, validity_end = item.validity
        degree = item.degree

        for item in items[1:]:
            if is_internal != item.is_internal:
                raise ValueError(
                    "Mixing of external and iternal coefficient sets!"
                )
            new_start, new_end = item.validity
            validity_start = max(validity_start, new_start)
            validity_end = min(validity_end, new_end)
            degree = max(degree, item.degree)

        SHCoefficients. __init__(
            self, is_internal=is_internal, validity_start=validity_start,
            validity_end=validity_end,
        )

        self._degree = degree
        self._items = items

    @property
    def degree(self):
        return self._degree

    def __call__(self, time, **parameters):
        degree = self.degree
        coeff_full = zeros((coeff_size(degree), 2))
        for item in self._items:
            item_coeff, item_degree, _ = item(time, **parameters)
            coeff_full[:coeff_size(item_degree), :] += item_coeff
        return coeff_full, degree, self.is_internal



class SparseSHCoefficients(SHCoefficients):
    """ Base class for sparse spherical harmonic coefficients. """
    def __init__(self, indices, **kwargs):
        SHCoefficients.__init__(self, **kwargs)
        n_idx, m_idx = indices[..., 0], indices[..., 1]
        self._degree_index = aabs(m_idx) + (n_idx*(n_idx + 1))//2
        self._coeff_index = (m_idx < 0).astype('int')
        self._degree = n_idx.max()

    @property
    def degree(self):
        return self._degree


class SparseSHCoefficientsConstant(SparseSHCoefficients):
    """ Time invariant sparse spherical harmonic coefficients. """
    def __init__(self, indices, coefficients, **kwargs):
        SparseSHCoefficients.__init__(self, indices, **kwargs)
        self._coeff = coefficients

    def __call__(self, time, **parameters):
        degree = self.degree
        coeff_full = zeros((coeff_size(degree), 2))
        coeff_full[self._degree_index, self._coeff_index] = self._coeff
        return coeff_full, degree, self.is_internal


class SparseSHCoefficientsTimeDependent(SparseSHCoefficients):
    """ Time dependent sparse spherical harmonic coefficients. """
    def __init__(self, indices, coefficients, times, **kwargs):
        order = argsort(times)
        self._coeff = coefficients[:, order]
        self._times = times[order]

        if kwargs.get("validity_start") is None:
            kwargs["validity_start"] = self._times[0]
        if kwargs.get("validity_end") is None:
            kwargs["validity_end"] = self._times[-1]

        SparseSHCoefficients.__init__(self, indices, **kwargs)

    def __call__(self, time, **parameters):
        degree = self.degree
        coeff_full = zeros((coeff_size(degree), 2))
        coeff_full[self._degree_index, self._coeff_index] = (
            self._interpolate_coefficients(time, self._coeff)
        )
        return coeff_full, degree, self.is_internal

    def _interpolate_coefficients(self, time, coeff):
        """ Return interpolated coefficients. """
        idx, basis = self._interpolation_basis(time)
        return dot(coeff[:, idx], basis).reshape(coeff.shape[0])

    def _interpolation_basis(self, time):
        """ Return indices and corresponding non-zero basis functions. """
        times = self._times
        idx1 = digitize([time], times)[0]
        idx1 = min(times.size - 1, max(1, idx1))
        idx0 = idx1 - 1
        alpha = (time -times[idx0])/(times[idx1] - times[idx0])
        return [idx0, idx1], array([1.0 - alpha, alpha]).reshape((2, 1))


def coeff_size(degree):
    """ Size of the full coefficient array. """
    return ((degree + 2)*(degree + 1))//2
