#-------------------------------------------------------------------------------
#
#  Magnetic Model
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

from numpy import asarray, empty, nditer
from .._pywmm import GRADIENT, GEOCENTRIC_SPHERICAL, sheval
from ..sheval_dipole import sheval_dipole


class GeomagneticModel(object):
    """ Abstract base class of the Earth magnetic field model. """
    # list of the required model parameters
    parameters = ("time", "location")

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):
        """ Evaluate magnetic field for the given MJD2000 times and coordinates.
        """
        raise NotImplementedError

    @property
    def validity(self):
        """ Get model's validity range as a tuple of two MJD2000 times.
        In case of an unconstrained validity rage (-inf, +inf) tuple is
        returned.
        """
        raise NotImplementedError


class SphericalHarmomicGeomagneticModel(GeomagneticModel):
    """ Earth magnetic field model calculated by the Spherical Harmonic
    Expansion.
    """

    def __init__(self, coefficients):
        self.coefficients = coefficients

    @property
    def validity(self):
        return self.coefficients.validity

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):
        time = asarray(time)
        if len(time.shape) == 0:
            _eval = self._eval_single_time
        else:
            _eval = self._eval_multi_time

        return _eval(
            time, location, input_coordinate_system, output_coordinate_system,
            **options
        )

    def _eval_multi_time(self, time, coords, input_coordinate_system,
                         output_coordinate_system, **options):
        """ Evaluate spherical harmonic for multiple times. """
        result = empty(coords.shape)
        if result.size > 0:
            iterator = nditer(
                [
                    time, coords[..., 0], coords[..., 1], coords[..., 2],
                    result[..., 0], result[..., 1], result[..., 2],
                ],
                op_flags=[
                    ['readonly'], ['readonly'], ['readonly'], ['readonly'],
                    ['writeonly'], ['writeonly'], ['writeonly'],
                ],
            )
            for time_, coord0, coord1, coord2, vect0, vect1, vect2 in iterator:
                vect0[...], vect1[...], vect2[...] = self._eval_single_time(
                    time_, [coord0, coord1, coord2], input_coordinate_system,
                    output_coordinate_system, **options
                )
        return result

    def _eval_single_time(self, time, coords, input_coordinate_system,
                          output_coordinate_system, **options):
        """ Evaluate spherical harmonic for a single time."""
        is_internal = self.coefficients.is_internal
        coeff, degree = self.coefficients(time, **options)
        return sheval(
            coords, degree, coeff[..., 0], coeff[..., 1],
            is_internal=is_internal, mode=GRADIENT,
            coord_type_in=input_coordinate_system,
            coord_type_out=output_coordinate_system,
            scale_gradient=-asarray(options.get("scale", 1.0))
        )


class DipoleSphericalHarmomicGeomagneticModel(SphericalHarmomicGeomagneticModel):
    """ Earth magnetic field model calculated by the Spherical Harmonic
    Expansion in the dipole coordinates.
    The dipole coordinates are defined by the time-dependent north pole
    latitude and longitude.
    """

    def __init__(self, coefficients, north_pole):
        SphericalHarmomicGeomagneticModel.__init__(self, coefficients)
        if callable(north_pole):
            # time dependent north pole
            self.north_pole = north_pole
        else:
            # time invariant north pole
            self.north_pole = lambda _: north_pole


    def _eval_single_time(self, time, coords, input_coordinate_system,
                          output_coordinate_system, **options):
        lat_ngp, lon_ngp = self.north_pole(time)
        is_internal = self.coefficients.is_internal
        coeff, degree = self.coefficients(
            time, lat_ngp=lat_ngp, lon_ngp=lon_ngp, **options
        )
        return sheval_dipole(
            coords, degree, coeff[..., 0], coeff[..., 1], lat_ngp, lon_ngp,
            is_internal=is_internal, mode=GRADIENT,
            coord_type_in=input_coordinate_system,
            coord_type_out=output_coordinate_system,
            scale_gradient=-asarray(options.get("scale", 1.0))
        )
