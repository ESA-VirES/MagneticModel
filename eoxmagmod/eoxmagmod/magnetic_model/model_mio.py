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
# pylint: disable=too-many-arguments

from numpy import nan, asarray, empty, isnan, full
from .._pymm import GRADIENT, GEOCENTRIC_SPHERICAL, convert, sheval2dfs
from ..magnetic_time import mjd2000_to_magnetic_universal_time
from .coefficients_mio import SparseSHCoefficientsMIO
from .model import GeomagneticModel, DipoleSphericalHarmomicGeomagneticModel
from .util import reshape_times_and_coordinates, reshape_array, mask_array

__all__ = ["DipoleMIOGeomagneticModel", "MIOPrimaryGeomagneticModel"]

MIO_HEIGHT = 110.0 # km
MIO_EARTH_RADIUS = 6371.2 # km
MIO_WOLF_RATIO = 0.014850


class MIOPrimaryGeomagneticModel(GeomagneticModel):
    """ Composed model switching between MIO primary fields evaluation
    above and below the Ionosphere.
    """
    parameters = ("time", "location", "f107", "subsolar_point")

    @property
    def degree(self):
        """ Get maximum degree of the model. """
        return max(
            self.model_below_ionosphere.degree,
            self.model_above_ionosphere.degree
        )

    @property
    def min_degree(self):
        """ Get minimum degree of the model. """
        return max(
            self.model_below_ionosphere.min_degree,
            self.model_above_ionosphere.min_degree
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
        time = asarray(time)
        location = convert(
            location, input_coordinate_system, GEOCENTRIC_SPHERICAL
        )

        radius_ionosphere = self.height + self.earth_radius
        radius = location[..., 2]
        mask_below_ionosphere = radius <= radius_ionosphere

        result = empty(location.shape)

        def _eval_model(model, mask, **options):

            def _mask_option(key):
                data = options.get(key)
                if data is not None:
                    options[key] = mask_array(asarray(data), mask)

            _mask_option('f107')
            _mask_option('lat_sol')
            _mask_option('lon_sol')

            result[mask] = model.eval(
                mask_array(time, mask), location[mask], GEOCENTRIC_SPHERICAL,
                output_coordinate_system, **options
            )

        _eval_model(
            self.model_below_ionosphere, mask_below_ionosphere, **options
        )
        _eval_model(
            self.model_above_ionosphere, ~mask_below_ionosphere, **options
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
        del height, earth_radius

        if not isinstance(coefficients, SparseSHCoefficientsMIO):
            raise TypeError(
                f"Invalid coefficients type {coefficients.__class__.__name__}! "
                f"{SparseSHCoefficientsMIO.__name__} is expected."
            )

        DipoleSphericalHarmomicGeomagneticModel.__init__(
            self, coefficients, north_pole
        )
        self.wolf_ratio = wolf_ratio

    def eval(self, time, location,
             input_coordinate_system=GEOCENTRIC_SPHERICAL,
             output_coordinate_system=GEOCENTRIC_SPHERICAL,
             **options):

        def _subset(data, mask):
            if data is None:
                return None
            data = asarray(data)
            if data.ndim > 0:
                return data[mask]
            return data

        # reshape time and location to a compatible shape
        time, location = reshape_times_and_coordinates(
            asarray(time), asarray(location)
        )

        # MIO scaling factor
        mask = True
        mio_scale = options.pop('f107', None)
        if mio_scale is not None:
            mio_scale = 1.0 + self.wolf_ratio * asarray(mio_scale)
            mask = ~isnan(mio_scale)

        start, end = self.validity
        mask = (time >= start) & (time <= end) & mask
        result = full(location.shape, nan)
        result[mask, :] = self._eval(
            self._eval_fourier2d, self.coefficients,
            time[mask], location[mask, :],
            input_coordinate_system, output_coordinate_system,
            lat_sol=_subset(options.pop('lat_sol', None), mask),
            lon_sol=_subset(options.pop('lon_sol', None), mask),
            **options,
        )

        if mio_scale is not None:
            if mio_scale.ndim > 0:
                mio_scale = reshape_array(location.shape, mio_scale)
            result *= mio_scale

        return result

    def _eval_fourier2d(self, coefficients, time, location,
                        input_coordinate_system, output_coordinate_system,
                        **options):
        """ Default SH expansion with MIO 2D Fourier series coefficients.  """

        year_fraction = coefficients.mjd2000_to_year_fraction(time)

        magnetic_universal_time = mjd2000_to_magnetic_universal_time(
            time, *(self.north_pole),
            lat_sol=options.pop('lat_sol', None),
            lon_sol=options.pop('lon_sol', None),
        )

        return sheval2dfs(
            year_fraction, magnetic_universal_time, location,
            coef_set=coefficients.get_f2fs_coeff_set(**options),
            coord_type_in=input_coordinate_system,
            coord_type_out=output_coordinate_system,
            mode=GRADIENT,
            is_internal=coefficients.is_internal,
            scale_gradient=-asarray(options.get("scale", 1.0)),
        )
