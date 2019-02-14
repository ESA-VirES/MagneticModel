#-------------------------------------------------------------------------------
#
#  MIO Magnetic Model
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
# pylint: disable=too-many-arguments,too-many-locals

from numpy import asarray, empty, nditer, logical_not, isscalar
from .._pymm import GEOCENTRIC_SPHERICAL, convert
from .model import GeomagneticModel, DipoleSphericalHarmomicGeomagneticModel

MIO_HEIGHT = 110.0 # km
MIO_EARTH_RADIUS = 6371.2 # km
MIO_WOLF_RATIO = 0.014850


class DipoleMIOPrimaryGeomagneticModel(GeomagneticModel):
    """ Composed model switching between MIO primary fields evaluation
    above and below the Ionosphere.
    """
    parameters = ("time", "location", "f107", "subsolar_point")

    @property
    def degree(self):
        return max(
            self.model_below_ionosphere.degree,
            self.model_above_ionosphere.degree
        )

    def __init__(self, model_below_ionosphere, model_above_ionosphere,
                 height=MIO_HEIGHT, earth_radius=MIO_EARTH_RADIUS):
        self.model_below_ionosphere = model_below_ionosphere
        self.model_above_ionosphere = model_above_ionosphere
        self.height = height
        self.earth_radius = earth_radius

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):
        location = convert(
            location, input_coordinate_system, GEOCENTRIC_SPHERICAL
        )

        radius_ionosphere = self.height + self.earth_radius
        radius = location[..., 2]
        mask_below_ionosphere = radius <= radius_ionosphere

        def _eval_masked_location(model, time, location, mask):
            result[mask] = model.eval(
                time, location[mask], GEOCENTRIC_SPHERICAL,
                output_coordinate_system, **options
            )

        def _eval_masked_time_and_location(model, time, location, mask):
            return _eval_masked_location(model, time[mask], location, mask)

        if isscalar(time):
            _eval_masked = _eval_masked_location
        else:
            _eval_masked = _eval_masked_time_and_location

        result = empty(location.shape)

        _eval_masked(
            self.model_below_ionosphere, time, location,
            mask_below_ionosphere
        )

        _eval_masked(
            self.model_above_ionosphere, time, location,
            logical_not(mask_below_ionosphere)
        )

        return result

    @property
    def validity(self):
        start_above, stop_above = self.model_above_ionosphere.validity
        start_below, stop_below = self.model_below_ionosphere.validity
        return (max(start_above, start_below), min(stop_above, stop_below))


class DipoleMIOGeomagneticModel(DipoleSphericalHarmomicGeomagneticModel):
    """ Swarm MIO model calculated by the Spherical Harmonic
    Expansion in the dipole coordinates.
    The dipole coordinates are defined by the north pole latitude and longitude.
    The model requires time dependent F10.7 index values.

    Options:
        coefficients - spherical harmonic coefficients
        north_pole - a tuple of dipole north pole lat/lon coordinates
                     or callable returning a tuple of dipole north pole
                     lat/lon coordinates for a given MJD2000 time.
        wolf_ratio - Wolf ratio (F10.7 scale)
        height - radius of the ionosphere (a + h)
        earth_radius - mean Earth radius used by the MIO model.
    """
    # list of the required model parameters
    parameters = ("time", "location", "f107", "subsolar_point")

    def __init__(self, coefficients, north_pole, wolf_ratio=MIO_WOLF_RATIO,
                 height=MIO_HEIGHT, earth_radius=MIO_EARTH_RADIUS):
        DipoleSphericalHarmomicGeomagneticModel.__init__(
            self, coefficients, north_pole
        )
        self.wolf_ratio = wolf_ratio
        self.height = height
        self.earth_radius = earth_radius

    def _eval_multi_time(self, time, coords, input_coordinate_system,
                         output_coordinate_system, f107=0.0,
                         lat_sol=None, lon_sol=None, **options):
        # TOOD: figure out a less clumsy way to iterate the optional arrays
        if lat_sol is None or lon_sol is None:
            return self._eval_multi_time_short(
                time, coords, input_coordinate_system, output_coordinate_system,
                f107, **options
            )
        else:
            return self._eval_multi_time_long(
                time, coords, input_coordinate_system, output_coordinate_system,
                f107, lat_sol, lon_sol, **options
            )

    def _eval_multi_time_short(self, time, coords, input_coordinate_system,
                               output_coordinate_system, f107, **options):
        result = empty(coords.shape)
        if result.size > 0:
            iterator = nditer(
                [
                    result[..., 0], result[..., 1], result[..., 2],
                    time, coords[..., 0], coords[..., 1], coords[..., 2],
                    asarray(f107),
                ],
                op_flags=[
                    ['writeonly'], ['writeonly'], ['writeonly'],
                    ['readonly'], ['readonly'], ['readonly'], ['readonly'],
                    ['readonly'],
                ],
            )
            for item in iterator:
                (
                    vect0, vect1, vect2,
                    time_, coord0, coord1, coord2, f107_,
                ) = item
                vect0[...], vect1[...], vect2[...] = self._eval_single_time(
                    time_, [coord0, coord1, coord2], input_coordinate_system,
                    output_coordinate_system, f107=f107_, **options
                )
        return result

    def _eval_multi_time_long(self, time, coords, input_coordinate_system,
                              output_coordinate_system, f107,
                              lat_sol, lon_sol, **options):
        result = empty(coords.shape)
        if result.size > 0:
            iterator = nditer(
                [
                    result[..., 0], result[..., 1], result[..., 2],
                    time, coords[..., 0], coords[..., 1], coords[..., 2],
                    asarray(f107), asarray(lat_sol), asarray(lon_sol),
                ],
                op_flags=[
                    ['writeonly'], ['writeonly'], ['writeonly'],
                    ['readonly'], ['readonly'], ['readonly'], ['readonly'],
                    ['readonly'], ['readonly'], ['readonly'],
                ],
            )
            for item in iterator:
                (
                    vect0, vect1, vect2,
                    time_, coord0, coord1, coord2, f107_, sslat, sslon,
                ) = item
                vect0[...], vect1[...], vect2[...] = self._eval_single_time(
                    time_, [coord0, coord1, coord2], input_coordinate_system,
                    output_coordinate_system, f107=f107_,
                    lat_sol=sslat, lon_sol=sslon, **options
                )
        return result

    def _eval_single_time(self, time, coords, input_coordinate_system,
                          output_coordinate_system, f107=0.0, **options):
        options["scale"] = (
            1.0 + self.wolf_ratio * f107
        ) * asarray(options.get("scale", 1.0))
        return DipoleSphericalHarmomicGeomagneticModel._eval_single_time(
            self, time, coords, input_coordinate_system,
            output_coordinate_system, **options
        )
