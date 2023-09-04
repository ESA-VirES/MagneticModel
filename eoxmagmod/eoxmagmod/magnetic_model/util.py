#-------------------------------------------------------------------------------
#
#  shared utilities
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

from bisect import bisect_right
from numpy import inf
from numpy.lib.stride_tricks import as_strided

__all__ = [
    "coeff_size",
    "convert_value",
    "parse_file",
    "reshape_times_and_coordinates",
    "mask_array",
    "reshape_array",
    "get_nonoverlapping_intervals",
    "aggregate_intersected_intervals",
]


def coeff_size(degree):
    """ Calculate size of the full coefficient array from the given degree. """
    return ((degree + 2)*(degree + 1))//2


def convert_value(value, conversion_function):
    """ Convert value using the optional conversion function. """
    if conversion_function is not None:
        value = conversion_function(value)
    return value


def parse_file(parser, file_, *args, **kwargs):
    """ Parse file with the provided file parser. """
    if isinstance(file_, str):
        with open(file_, encoding="ascii") as file_in:
            return parser(file_in, *args, **kwargs)
    else:
        return parser(file_, *args, **kwargs)


def reshape_times_and_coordinates(time, coords):
    """ Reshape coordinates to match the times. """
    ndim_common = min(time.ndim, coords.ndim-1)
    ndim_result = max(time.ndim, coords.ndim-1)

    if time.shape[:ndim_common] != coords.shape[:ndim_common]:
        raise ValueError(
            "Incompatible dimensions of the time and coordinates arrays."
        )

    if coords.ndim - 1 < ndim_result:
        coords = as_strided(
            coords,
            shape=(*time.shape, coords.shape[-1]),
            strides=(
                *coords.strides[:-1],
                *((0,) * (ndim_result - ndim_common)),
                coords.ndim[-1]
            )
        )

    if 0 < time.ndim < ndim_result:
        time = as_strided(
            time,
            shape=coords.shape[:-1],
            strides=(*time.strides, *((0,) * (ndim_result - ndim_common)))
        )

    return time, coords


def mask_array(data, mask, min_ndim=0):
    """ Apply mask if data is a muti-value array or pass through
    a single value.
    """
    return data[mask] if data.ndim > min_ndim else data


def reshape_array(shape, array):
    """ Reshape array to match the given shape. """
    ndim_common = min(len(shape), array.ndim)

    if (
        shape[:ndim_common] != array.shape[:ndim_common]
        or len(shape) < array.ndim
    ):
        raise ValueError(
            "Incompatible dimensions of the source and reshaped array arrays."
        )

    if array.ndim < len(shape):
        array = as_strided(
            array,
            shape=shape,
            strides=(*array.strides, *((0,) * (len(shape) - array.ndim)))
        )

    return array


def get_nonoverlapping_intervals(intervals):
    """ Get a list of non-overlapping intervals from the given list of
    intervals and their associated objects.

    Input:
        intervals - list of ((start, end), obj) tuples

    Output:
        list of ((start, end), obj) tuples
    """
    starts, ends, items = [], [], []

    def _insert(idx, start, end, item):
        starts.insert(idx, start)
        ends.insert(idx, end)
        items.insert(idx, item)

    for (start, end), item in intervals:
        if start > end: # skip invalid intervals
            continue

        while start < end:
            idx = bisect_right(ends, start)
            if idx == len(starts):
                if idx == 0 or ends[idx-1] < end:
                    _insert(idx, start, end, item)
                break
            if starts[idx] > start:
                if starts[idx] >= end:
                    _insert(idx, start, end, item)
                    break
                _insert(idx, start, starts[idx], item)
                idx += 1
            start = ends[idx]

    return list(zip(zip(starts, ends), items))


def aggregate_intersected_intervals(intervals0, aggregates=None):
    """ Aggregate associated objects of the intersected non-overlapping
    intervals.

    Inputs:
        intervals - list of ((start, end), obj) tuples
        aggregates - list of ((start, end), [obj, ...]) tuples or None

    Output:
        list of ((start, end), [obj, ...]) tuples
    """
    if aggregates is None:
        aggregates = [((-inf, inf), [])]

    it_aggregates = iter(aggregates)
    it_intervals0 = iter(intervals0)

    result = []

    def _insert(start, end, items, item):
        result.append(((start, end), [*items, item]))

    def _intersect(start1, end1, start0, end0):
        if start1 <= end0 and start0 <= end1:
            return max(start1, start0), min(end1, end0)
        return None

    try:
        (start1, end1), items = next(it_aggregates)
        (start0, end0), item = next(it_intervals0)

        while True:
            intersection = _intersect(start1, end1, start0, end0)

            if intersection:
                _insert(*intersection, items, item)

            if end0 < end1:
                (start0, end0), item = next(it_intervals0)
            else:
                (start1, end1), items = next(it_aggregates)

    except StopIteration:
        pass

    return result
