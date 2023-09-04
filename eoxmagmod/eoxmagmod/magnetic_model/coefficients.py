#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Coefficients.
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
# pylint: disable=too-few-public-methods,abstract-method

from numpy import (
    inf, asarray, zeros, ones, argsort, abs as aabs, stack,
)
from ..time_util import decimal_year_to_mjd2000, mjd2000_to_decimal_year
from .._pymm import interp, INTERP_C0, INTERP_C1
from .util import (
    get_nonoverlapping_intervals, aggregate_intersected_intervals,
    coeff_size, convert_value,
)

__all__ = [
    "SHCoefficients",
    "ComposedSHCoefficients",
    "CombinedSHCoefficients",
    "SparseSHCoefficients",
    "SparseSHCoefficientsTimeDependent",
    "SparseSHCoefficientsConstant",
    "SparseSHCoefficientsTimeDependentDecimalYear",
]


class SHCoefficients:
    """ Abstract base class for all spherical harmonic coefficients. """
    time_scales = ()

    def decompose(self):
        """ Return a sequence of decomposed coefficient sets. """
        raise NotImplementedError

    def __init__(self, is_internal=True, **kwargs):
        self.is_internal = is_internal
        self.validity = self._get_converted_validity(**kwargs)

    @staticmethod
    def _get_converted_validity(validity_start=None, validity_end=None,
                                to_mjd2000=None, **_):
        """ Get validity range converted to MJD2000 using the optionally
        provided conversion function.
        """

        def _get_converted_value_or_default(value, default, conversion_function):
            if value is None:
                value = default
            if conversion_function is not None:
                value = conversion_function(value)
            return value

        return (
            _get_converted_value_or_default(validity_start, -inf, to_mjd2000),
            _get_converted_value_or_default(validity_end, +inf, to_mjd2000),
        )


    def is_valid(self, time):
        """ Check if the time is within the coefficients validity range. """
        validity_start, validity_end = self.validity
        return (validity_start <= time) & (time <= validity_end)

    @property
    def degree(self):
        """ Get [maximum] model degree. """
        raise NotImplementedError

    @property
    def min_degree(self):
        """ Get minimum model degree.
        Below this degree all model coefficients are zero.
        """
        raise NotImplementedError

    def __call__(self, time, **parameters):
        """ Return the matrix of the full model coefficients and its degree. """
        raise NotImplementedError


class ComposedSHCoefficients(SHCoefficients):
    """ Model composed of a sequence of multiple time-non-overlapping
    coefficient sets, i.e., there is not more than one model set applicable
    for the requested time instant.

    The order in which the coefficients are composed matters. The first set
    of coefficients matching the requested time instance is evaluated.
    """

    def decompose(self):
        result = []
        for item in self:
            result.extend(get_nonoverlapping_intervals(item.decompose()))
        return result

    @property
    def time_scales(self):
        """ Time scales of the composed models. """
        time_scales = set()
        for item in self:
            time_scales.update(item.time_scales)
        return tuple(time_scales)

    def __new__(cls, *items):
        if len(items) == 1:
            return items[0]
        return super(ComposedSHCoefficients, cls).__new__(cls)

    def __iter__(self):
        return iter(self._items)

    def __init__(self, *items):

        if len(items) < 1:
            raise ValueError(
                "The composed model requires at least one set of coefficients!"
            )

        validity_start, validity_end = inf, -inf
        last_item_validity_end = inf
        degree = -inf
        min_degree = inf
        is_internal = items[0].is_internal

        for item in sorted(items, key=lambda item: item.validity):
            item_validity_start, item_validity_end = item.validity

            if item_validity_start > last_item_validity_end:
                # time-gaps between the models are not allowed
                raise ValueError(
                    "The composed model does not allow time-gaps between "
                    "the sets of coefficients!"
                )

            last_item_validity_end = item_validity_end

            validity_start = min(validity_start, item_validity_start)
            validity_end = max(validity_end, item_validity_end)

            degree = max(degree, item.degree)
            min_degree = min(min_degree, item.min_degree)

            if is_internal != item.is_internal:
                raise ValueError(
                    "Mixing of external and internal coefficient sets"
                    "is not allowed!"
                )

        SHCoefficients.__init__(
            self,
            is_internal=is_internal,
            validity_start=validity_start,
            validity_end=validity_end,
        )

        self._degree = degree
        self._min_degree = min_degree
        self._items = items

    @property
    def degree(self):
        return self._degree

    @property
    def min_degree(self):
        return self._min_degree

    def __call__(self, time, **parameters):
        time = asarray(time)
        if time.ndim == 0:
            return self._get_coeff_single_time(time, **parameters)
        return self._get_coeff_multi_time(time, **parameters)

    def _get_coeff_multi_time(self, time, **parameters):
        mask = ones(time.shape, 'bool')

        # process parts matched by different coefficient sets ...
        degree = 0
        results = []
        for item in self._items:
            item_mask = mask & item.is_valid(time)
            item_coeff, item_degree = item(time[item_mask], **parameters)
            results.append((item_mask, item_coeff))
            degree = max(degree, item_degree)
            mask &= ~item_mask
        item_coeff, item_degree = self._items[0](time[mask], **parameters)
        results.append((mask, item_coeff))
        degree = max(degree, item_degree)

        # merge parts
        coeff = zeros((*time.shape, coeff_size(degree), 2))
        for item_mask, item_coeff in results:
            coeff[item_mask, :item_coeff.shape[-2], :] = item_coeff

        return coeff, degree

    def _get_coeff_single_time(self, time, **parameters):
        for item in self._items:
            if item.is_valid(time):
                return item(time, **parameters)
        return self._items[0](time, **parameters)


class CombinedSHCoefficients(SHCoefficients):
    """ Coefficients combined of multiple time-overlapping coefficient sets.
    These sets are evaluated together to form a single set of coefficients
    and for a time instant all combined coefficient sets are applicable.

    The combined coefficients can be used e.g., to merge time-variable core and
    constant lithospheric model coefficients.
    """

    def decompose(self):
        aggergates = None
        for item in self:
            aggergates = aggregate_intersected_intervals(
                item.decompose(), aggergates
            )
        result = []
        for validity, items in aggergates:
            result.append((validity, CombinedSHCoefficients(*items)))
        return result

    @property
    def time_scales(self):
        """ Time scales of the combined models. """
        time_scales = set()
        for item in self:
            time_scales.update(item.time_scales)
        return tuple(time_scales)

    def __new__(cls, *items):
        if len(items) == 1:
            return items[0]
        return super(CombinedSHCoefficients, cls).__new__(cls)

    def __iter__(self):
        return iter(self._items)

    def __init__(self, *items):

        if len(items) < 1:
            raise ValueError(
                "The combined model requires at least one set of coefficients!"
            )

        item = items[0]
        is_internal = item.is_internal
        validity_start, validity_end = item.validity
        degree = item.degree
        min_degree = item.min_degree

        for item in items[1:]:
            if is_internal != item.is_internal:
                raise ValueError(
                    "Mixing of external and internal coefficient sets"
                    "is not allowed!"
                )
            new_start, new_end = item.validity
            validity_start = max(validity_start, new_start)
            validity_end = min(validity_end, new_end)
            degree = max(degree, item.degree)
            min_degree = min(min_degree, item.min_degree)

        SHCoefficients. __init__(
            self,
            is_internal=is_internal,
            validity_start=validity_start,
            validity_end=validity_end,
        )

        self._degree = degree
        self._min_degree = min_degree
        self._items = items

    @property
    def degree(self):
        return self._degree

    @property
    def min_degree(self):
        return self._min_degree

    def get_coefficient_sets(self, **parameters):
        """ Extract coefficients sets which can be passed to C evaluation. """
        coeff_sets = []
        for item in self._items:
            coeff_sets.extend(item.get_coefficient_sets(**parameters))
        return coeff_sets

    @property
    def convert_time(self):
        """ Extract time-conversion function, if available. """
        for item in self._items:
            convert_time = getattr(item, "convert_time", None)
            if convert_time:
                return convert_time
        return None

    def __call__(self, time, **parameters):
        time = asarray(time)
        # evaluate coefficient sets ...
        degree = 0
        results = []
        for item in self._items:
            item_coeff, item_degree = item(time, **parameters)
            results.append(item_coeff)
            degree = max(degree, item_degree)
        # merge coefficient sets ...
        coeff = zeros((*time.shape, coeff_size(degree), 2))
        for item_coeff in results:
            coeff[..., :item_coeff.shape[-2], :] += item_coeff
        return coeff, degree



class SparseSHCoefficients(SHCoefficients):
    """ Base class for sparse spherical harmonic coefficients.
    Each sparse coefficient has a degree and order defining its position
    in the full coefficient array.
    """
    def decompose(self):
        return [(self.validity, self)]

    def __init__(self, indices, coefficients, **kwargs):
        SHCoefficients.__init__(self, **kwargs)
        n_idx, m_idx = indices[..., 0], indices[..., 1]
        self._degree = n_idx.max()
        self._min_degree = n_idx.min()
        self._nm = indices.astype('int32')
        self._index = stack((
            aabs(m_idx) + (n_idx*(n_idx + 1))//2,
            (m_idx < 0).astype('int'),
        ), -1)
        self._coeff = coefficients

    @property
    def degree(self):
        return self._degree

    @property
    def min_degree(self):
        return self._min_degree

    def subset_degree(self, min_degree, max_degree):
        """ Get subset of the coefficients for the give minimum and maximum
        degrees.
        """
        default_max_degree = self._degree
        default_min_degree = self._min_degree

        nm_ = self._nm
        index = self._index
        coeff = self._coeff

        if max_degree < 0:
            max_degree = default_max_degree

        if min_degree < 0:
            min_degree = default_min_degree

        if min_degree > default_min_degree or max_degree < default_max_degree:
            idx, = (
                (nm_[:, 0] <= max_degree) & (nm_[:, 0] >= min_degree)
            ).nonzero()
            coeff = coeff[idx]
            nm_ = nm_[idx]
            index = index[idx]
            degree = default_max_degree

            if index.shape[0] > 0:
                degree = nm_[:, 0].max()
            else:
                degree = 0
        else:
            degree = default_max_degree

        return degree, coeff, nm_, index


class SparseSHCoefficientsTimeDependent(SparseSHCoefficients):
    """ Time dependent sparse spherical harmonic coefficients evaluated
    by spine interpolation of a time series of coefficients snapshots
    interpolated in the MJD2000 time domain.
    """
    time_scales = ("MJD2000",)

    def __init__(self, indices, coefficients, times, spline_order=2, **kwargs):
        indices = asarray(indices)
        coefficients = asarray(coefficients)
        times = asarray(times)

        # check array dimensions
        if times.ndim != 1:
            raise ValueError(f"Invalid times array dimension {times.shape}!")

        if coefficients.ndim != 2:
            raise ValueError(f"Invalid coefficients array dimension {coefficients.shape}!")

        if indices.ndim != 2 and indices.shape[1] != 2:
            raise ValueError(f"Invalid indices array dimension {indices.shape}!")

        if indices.shape[0] != coefficients.shape[0]:
            raise ValueError(
                "Shape mismatch of the indices and coefficients arrays "
                f"({indices.shape} vs. {coefficients.shape})!"
            )

        if times.shape[0] != coefficients.shape[1]:
            raise ValueError(
                "Shape mismatch of the times and coefficients arrays "
                f"({times.shape} vs. {coefficients.shape})!"
            )

        # assure sort times
        order = argsort(times)
        times = times[order]
        coefficients = coefficients[:, order]

        # extract temporal validity
        kwargs['validity_start'] = kwargs.get('validity_start', times[0])
        kwargs['validity_end'] = kwargs.get('validity_end', times[-1])

        SparseSHCoefficients.__init__(self, indices, coefficients, **kwargs)

        self._times = convert_value(times, kwargs.get("to_mjd2000"))
        self._spline_order = spline_order

        try:
            self._interp_options = {
                1: {"kind": INTERP_C0},
                2: {"kind": INTERP_C1},
            }[spline_order]
        except KeyError:
            raise ValueError("Unsupported spline order {spline_order}!") from None

    def get_coefficient_sets(self, **parameters):
        """ Extract coefficients sets which can be passed to C evaluation. """
        _, coeff, nm_, _ = self.subset_degree(
            parameters.get("min_degree", -1),
            parameters.get("max_degree", -1),
        )
        return [(self._times, coeff, nm_, self._spline_order)]

    def __call__(self, time, **parameters):
        time = convert_value(asarray(time), getattr(self, "convert_time", None))
        degree, coeff, _, index = self.subset_degree(
            parameters.get("min_degree", -1),
            parameters.get("max_degree", -1),
        )
        coeff_full = zeros((*time.shape, coeff_size(degree), 2))
        coeff_full[..., index[:, 0], index[:, 1]] = interp(
            time, self._times, coeff, extrapolate=True, **self._interp_options
        )
        return coeff_full, degree


class SparseSHCoefficientsConstant(SparseSHCoefficientsTimeDependent):
    """ Time invariant sparse spherical harmonic coefficients. """
    time_scale = ()

    def __init__(self, indices, coefficients, **kwargs):
        indices = asarray(indices)
        coefficients = asarray(coefficients)

        # check array dimensions
        if coefficients.ndim != 1:
            raise ValueError(f"Invalid coefficients array dimension {coefficients.shape}!")

        if indices.ndim != 2 and indices.shape[1] != 2:
            raise ValueError(f"Invalid indices array dimension {indices.shape}!")

        if indices.shape[0] != coefficients.shape[0]:
            raise ValueError(
                "Shape mismatch of the indices and coefficients arrays "
                f"({indices.shape} vs. {coefficients.shape})!"
            )

        # extract temporal validity
        kwargs['validity_start'] = kwargs.get('validity_start', -inf)
        kwargs['validity_end'] = kwargs.get('validity_end', +inf)

        SparseSHCoefficientsTimeDependent.__init__(
            self, indices, coefficients.reshape((*coefficients.shape, 1)),
            times=[-inf], spline_order=1, **kwargs
        )


class SparseSHCoefficientsTimeDependentDecimalYear(SparseSHCoefficientsTimeDependent):
    """ Time dependent sparse spherical harmonic coefficients evaluated
    by spine interpolation of a time series of coefficients snapshots.
    with time in decimal year interpolated in the decimal years time domain.
    """
    time_scales = ("MJD2000",)

    def __init__(self, indices, coefficients, times,
                 to_mjd2000=decimal_year_to_mjd2000,
                 to_decimal_year=mjd2000_to_decimal_year,
                 **kwargs):
        SparseSHCoefficientsTimeDependent.__init__(
            self, indices, coefficients, times, **kwargs
        )
        self.convert_time = to_decimal_year
        # Fix the validity range to be in the expected MJD2000.
        self.validity = self._get_converted_validity(
            *self.validity, to_mjd2000=to_mjd2000
        )
