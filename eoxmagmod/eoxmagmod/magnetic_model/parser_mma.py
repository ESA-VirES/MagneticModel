#-------------------------------------------------------------------------------
#
#  Swarm MMA_SHA_2C product file format parser
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

from numpy import array

CDF_EPOCH_TYPE = 31 # spacepy.pycdf.const.CDF_EPOCH.value
CDF_EPOCH_2000 = 63113904000000.0
CDF_EPOCH_TO_DAYS = 1.0/86400000.0


def read_swarm_mma_2c_file(cdf):
    """ Read Swarm MIO_SHA_2C product CDF file and returns a dictionary
    containing the parsed model data.

    The function expect a spacepy.pycdf.CDF object.
    """
    return {
        "gh": (
            read_swarm_mma_2c_coefficients(cdf, "gh_1", "gh"),
            read_swarm_mma_2c_coefficients(cdf, "gh_2", "gh"),
        ),
        "qs": (
            read_swarm_mma_2c_coefficients(cdf, "qs_1", "qs"),
            read_swarm_mma_2c_coefficients(cdf, "qs_2", "qs"),
        ),
    }


def read_swarm_mma_2c_coefficients(cdf, source_variable, target_variable):
    """ Read a single set of Swarm MIO_SHA_2C coefficients.

    The function expect a spacepy.pycdf.CDF object.
    """
    time = cdf.raw_var("t_" + source_variable)
    n_idx = cdf["nm_" + source_variable][0, :, 0]
    return {
        "degree_min": n_idx.min(),
        "degree_max": n_idx.max(),
        "t": _cdf_rawtime_to_mjd2000(time[0], time.type()),
        "nm": array(cdf["nm_" + source_variable][0]),
        target_variable: array(cdf[source_variable][0].transpose()),
    }


def _cdf_rawtime_to_mjd2000(raw_time, cdf_type):
    """ Convert an array of CDF raw time values to array of MJD2000 values.
    """
    if cdf_type == CDF_EPOCH_TYPE:
        return (raw_time - CDF_EPOCH_2000) * CDF_EPOCH_TO_DAYS
    else:
        raise TypeError("Unsupported CDF time type %r !" % cdf_type)
