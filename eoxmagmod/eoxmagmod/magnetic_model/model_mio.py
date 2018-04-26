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

from numpy import asarray
from .model import DipoleSphericalHarmomicGeomagneticModel, GEOCENTRIC_SPHERICAL

MIO_HEIGHT = 110.0 # km
MIO_EARTH_RADIUS = 6371.2 # km
MIO_WOLF_RATIO = 0.014850


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
        f107 - F10.7 index value or a callable returning F10.7 index value
               for a given given MJD2000 time.
        wolf_ratio - Wolf ratio (F10.7 scale)
        height - radius of the ionosphere (a + h)
        earth_radius - mean Earth radius used by the MIO model.
    """

    def __init__(self, coefficients, north_pole, f107, wolf_ratio=MIO_WOLF_RATIO,
                 height=MIO_HEIGHT, earth_radius=MIO_EARTH_RADIUS):
        DipoleSphericalHarmomicGeomagneticModel.__init__(
            self, coefficients, north_pole
        )
        if callable(f107):
            # time dependent F10.7 index
            self.f107 = f107
        else:
            # time invariant F10.7 index
            self.f107 = lambda _: f107
        self.wolf_ratio = wolf_ratio
        self.height = height
        self.earth_radius = earth_radius

    def _eval_single_time(self, time, coords,
                          input_coordinate_system=GEOCENTRIC_SPHERICAL,
                          output_coordinate_system=GEOCENTRIC_SPHERICAL,
                          **options):
        options["scale"] = (
            1.0 + self.wolf_ratio * self.f107(time)
        ) * asarray(options.get("scale", 1.0))
        return DipoleSphericalHarmomicGeomagneticModel._eval_single_time(
            self, time, coords, input_coordinate_system,
            output_coordinate_system, **options
        )
