#-------------------------------------------------------------------------------
#
#  Magnetic Model
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
# pylint: disable=too-many-locals, too-many-arguments, no-self-use


from numpy import asarray, empty, full, nan, nditer
from .._pymm import GRADIENT, GEOCENTRIC_SPHERICAL, sheval, shevaltemp
from ..dipole import rotate_vectors_from_dipole
from ..dipole_coords import convert_to_dipole
from .util import reshape_times_and_coordinates
from .coefficients import (
    SparseSHCoefficientsTimeDependent,
    CombinedSHCoefficients,
)

__all__ = [
    "GeomagneticModel",
    "SphericalHarmomicGeomagneticModel",
    "DipoleSphericalHarmomicGeomagneticModel",
]


class GeomagneticModel:
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
        self._components = self._decompose_coefficients(coefficients)

    def _decompose_coefficients(self, coefficients):
        components = []
        for validity, item in coefficients.decompose():
            method = self._select_evaluation_method(item)
            components.append((validity, method, item))
        return components

    def _select_evaluation_method(self, coefficients):
        """ Choose the best evaluation method. """
        if isinstance(coefficients, SparseSHCoefficientsTimeDependent):
            return self._eval_interpolated

        if isinstance(coefficients, CombinedSHCoefficients):
            if len(coefficients.time_scales) > 1:
                # mixed time-scales - falling back to the default
                return self._eval_default

            try:
                coefficients.get_coefficient_sets()
            except AttributeError:
                # coefficient cannot be extracted - falling back to the default
                return self._eval_default

            return self._eval_interpolated

        # default evaluation method
        return self._eval_default

    @property
    def validity(self):
        return self.coefficients.validity

    @property
    def degree(self):
        """ Get maximum degree of the model. """
        return self.coefficients.degree

    @property
    def min_degree(self):
        """ Get minimum degree of the model. """
        return self.coefficients.min_degree

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):

        # reshape time and location to a compatible shape
        time, location = reshape_times_and_coordinates(
            asarray(time), asarray(location)
        )

        def _is_valid(start, end):
            return (time >= start) & (time <= end)

        def _call_eval(method, coefficients, time, location):
            return self._eval(
                method, coefficients, time, location,
                input_coordinate_system, output_coordinate_system, **options
            )

        def _get_nans():
            return full(location.shape, nan)

        if time.ndim == 0:
            # time is a scalar value
            for (start, end), method, coefficients in self._components:
                if _is_valid(start, end):
                    return _call_eval(method, coefficients, time, location)
            return _get_nans()

        # time is an array
        result = _get_nans()
        for (start, end), method, coefficients in self._components:
            mask = _is_valid(start, end)
            result[mask, :] = _call_eval(
                method, coefficients, time[mask], location[mask, :]
            )
        return result

    def _eval(self, method, coefficients, time, location,
              input_coordinate_system, output_coordinate_system, **options):
        return method(
            coefficients, time, location,
            input_coordinate_system, output_coordinate_system, **options
        )

    def _eval_interpolated(self, coefficients, time, location,
                           input_coordinate_system, output_coordinate_system,
                           **options):
        """ Optimized SH expansion with time-series interpolation. """
        convert_time = getattr(coefficients, "convert_time", None)
        if convert_time:
            time = convert_time(time)
        return shevaltemp(
            time,
            location,
            coef_set_list=coefficients.get_coefficient_sets(**options),
            coord_type_in=input_coordinate_system,
            coord_type_out=output_coordinate_system,
            mode=GRADIENT,
            is_internal=coefficients.is_internal,
            scale_gradient=-asarray(options.get("scale", 1.0)),
        )

    def _eval_default(self, coefficients, time, location,
                      input_coordinate_system, output_coordinate_system,
                      **options):
        """ Default iterated single-time interpolation. """

        eval_model = (
            self._eval_multi_time if time.ndim > 0 else self._eval_single_time
        )

        return eval_model(
            coefficients, time, location,
            input_coordinate_system, output_coordinate_system, **options
        )

    def _eval_multi_time(self, coefficients, time, coords,
                         input_coordinate_system, output_coordinate_system,
                         **options):
        """ Evaluate spherical harmonic model for multiple times. """
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
                    coefficients, time_, [coord0, coord1, coord2],
                    input_coordinate_system, output_coordinate_system,
                    **options
                )
        return result

    def _eval_single_time(self, coefficients, time, coords,
                          input_coordinate_system, output_coordinate_system,
                          **options):
        """ Evaluate spherical harmonic for a single time."""
        coeff, degree = coefficients(time, **options)
        return sheval(
            coords, coeff,
            mode=GRADIENT, degree=degree, is_internal=coefficients.is_internal,
            coord_type_in=input_coordinate_system,
            coord_type_out=output_coordinate_system,
            scale_gradient=-asarray(options.get("scale", 1.0))
        )


class DipoleSphericalHarmomicGeomagneticModel(SphericalHarmomicGeomagneticModel):
    """ Earth magnetic field model calculated by the Spherical Harmonic
    Expansion in the dipole coordinates.
    The dipole coordinates are defined north latitude and longitude of
    the dipole axis.
    """

    def __init__(self, coefficients, north_pole):
        super().__init__(coefficients)
        self.north_pole = north_pole

    def _eval(self, method, coefficients, time, location,
              input_coordinate_system, output_coordinate_system, **options):

        lat_ngp, lon_ngp = self.north_pole
        scale = options.pop('scale', None)

        location_dipole = convert_to_dipole(
            location, lat_ngp, lon_ngp, input_coordinate_system
        )

        result = super()._eval(
            method, coefficients, time, location_dipole,
            GEOCENTRIC_SPHERICAL, GEOCENTRIC_SPHERICAL,
            **options
        )

        result = rotate_vectors_from_dipole(
            result, lat_ngp, lon_ngp, location_dipole, location,
            input_coordinate_system, output_coordinate_system
        )

        if scale is not None:
            result *= scale

        return result
