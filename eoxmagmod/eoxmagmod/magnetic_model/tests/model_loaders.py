#-------------------------------------------------------------------------------
#
#  Model loaders - tests
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
# pylint: disable=missing-docstring,no-self-use

from unittest import TestCase, main
from itertools import product
from numpy import inf, array, empty, nditer, asarray
from numpy.random import uniform
from numpy.testing import assert_allclose
from eoxmagmod import decimal_year_to_mjd2000
from eoxmagmod.magnetic_model.loader_shc import (
    load_model_shc, load_model_shc_combined,
)
from eoxmagmod.magnetic_model.loader_igrf import load_model_igrf
from eoxmagmod.magnetic_model.loader_wmm import load_model_wmm
from eoxmagmod.magnetic_model.loader_emm import load_model_emm
from eoxmagmod.magnetic_model.loader_mma import (
    load_model_swarm_mma_2c_internal,
    load_model_swarm_mma_2c_external,
)
from eoxmagmod.magnetic_model.loader_mio import (
    load_model_swarm_mio_internal,
    load_model_swarm_mio_external,
)
from eoxmagmod.data import (
    EMM_2010_STATIC, EMM_2010_SECVAR, WMM_2010, WMM_2015,
    CHAOS5_CORE, CHAOS5_CORE_V4, CHAOS5_STATIC,
    CHAOS6_CORE, CHAOS6_CORE_X3, CHAOS6_STATIC,
    IGRF11, IGRF12, SIFM,
)
from eoxmagmod.magnetic_model.tests.data import (
    SWARM_MMA_SHA_2C_TEST_DATA,
    SWARM_MIO_SHA_2_TEST_DATA,
)
from eoxmagmod.magnetic_model.model import (
    SphericalHarmomicGeomagneticModel,
    DipoleSphericalHarmomicGeomagneticModel,
)
from eoxmagmod.magnetic_model.model_mio import DipoleMIOGeomagneticModel
from eoxmagmod._pywmm import GEOCENTRIC_SPHERICAL, sheval, GRADIENT
from eoxmagmod.sheval_dipole import sheval_dipole


class SHModelTestMixIn(object):
    scale = 1.0
    range_lat = range(-90, 91, 5)
    range_lon = range(-180, 181, 10)
    validity = None
    model_class = SphericalHarmomicGeomagneticModel
    options = {}

    @property
    def model(self):
        if not hasattr(self, "_model"):
            self._model = self.load()
        return self._model

    @property
    def coordinates(self):
        return array([
            (lat, lon, 6371.2*uniform(1.0, 2.0)) for lat, lon
            in product(self.range_lat, self.range_lon)
        ])

    @staticmethod
    def _constrain_validity(time_min, time_max):
        return (
            max(time_min, decimal_year_to_mjd2000(1920)),
            min(time_max, decimal_year_to_mjd2000(2120)),
        )

    @property
    def time(self):
        return uniform(*self._constrain_validity(*self.validity))

    @property
    def times(self):
        time_min, time_max = self._constrain_validity(*self.validity)
        return uniform(
            time_min, time_max, len(self.range_lat) * len(self.range_lon)
        )

    def eval_model(self, times, coords):
        return self.model.eval(times, coords, **self.options)

    def eval_reference(self, times, coords):
        result = empty(coords.shape)
        iterator = nditer(
            [
                times, coords[..., 0], coords[..., 1], coords[..., 2],
                result[..., 0], result[..., 1], result[..., 2],
            ],
            op_flags=[
                ['readonly'], ['readonly'], ['readonly'], ['readonly'],
                ['writeonly'], ['writeonly'], ['writeonly'],
            ],
        )
        for time, coord0, coord1, coord2, vect0, vect1, vect2 in iterator:
            vect0[...], vect1[...], vect2[...] = self._eval_reference(
                time, [coord0, coord1, coord2]
            )
        return result

    def _eval_reference(self, time, coords):
        is_internal = self.model.coefficients.is_internal
        coeff, degree = self.model.coefficients(time, **self.options)
        return sheval(
            coords, degree, coeff[..., 0], coeff[..., 1],
            is_internal=is_internal, mode=GRADIENT,
            coord_type_in=GEOCENTRIC_SPHERICAL,
            coord_type_out=GEOCENTRIC_SPHERICAL,
            scale_gradient=-asarray(self.scale),
        )

    def test_class(self):
        self.assertIsInstance(self.model, self.model_class)

    def test_validity(self):
        assert_allclose(self.model.validity, self.validity)

    def test_eval_single_time(self):
        time = self.time
        coords = self.coordinates
        #self.eval_model(time, coords)
        assert_allclose(
            self.eval_model(time, coords),
            self.eval_reference(time, coords),
        )

    def test_eval_multi_time(self):
        times = self.times
        coords = self.coordinates
        #self.eval_model(times, coords)
        assert_allclose(
            self.eval_model(times, coords),
            self.eval_reference(times, coords),
        )


class DipoleSHModelTestMixIn(SHModelTestMixIn):
    model_class = DipoleSphericalHarmomicGeomagneticModel

    def _eval_reference(self, time, coords):
        is_internal = self.model.coefficients.is_internal
        lat_ngp, lon_ngp = self.model.north_pole(time)
        coeff, degree = self.model.coefficients(time)
        return sheval_dipole(
            coords, degree, coeff[..., 0], coeff[..., 1], lat_ngp, lon_ngp,
            is_internal=is_internal, mode=GRADIENT,
            coord_type_in=GEOCENTRIC_SPHERICAL,
            coord_type_out=GEOCENTRIC_SPHERICAL,
            scale_gradient=-asarray(self.scale),
        )

class DipoleMIOSHModelTestMixIn(SHModelTestMixIn):
    f107 = 70.0
    model_class = DipoleMIOGeomagneticModel

    def _eval_reference(self, time, coords):
        is_internal = self.model.coefficients.is_internal
        scale = -(
            1.0 + self.model.f107(time) * self.model.wolf_ratio
        ) * asarray(self.scale)
        lat_ngp, lon_ngp = self.model.north_pole(time)
        coeff, degree = self.model.coefficients(time, lat_ngp, lon_ngp)
        return sheval_dipole(
            coords, degree, coeff[..., 0], coeff[..., 1], lat_ngp, lon_ngp,
            is_internal=is_internal, mode=GRADIENT,
            coord_type_in=GEOCENTRIC_SPHERICAL,
            coord_type_out=GEOCENTRIC_SPHERICAL,
            scale_gradient=scale
        )

#-------------------------------------------------------------------------------

class TestWMM2010(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((2010.0, 2015.0))
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]
    def load(self):
        return load_model_wmm(WMM_2010)


class TestWMM2015(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((2015.0, 2020.0))
    def load(self):
        return load_model_wmm(WMM_2015)


class TestEMM2010(TestCase, SHModelTestMixIn):
    # The EMM models is huge and the test ranges have to be reduced.
    range_lat = range(-90, 91, 30)
    range_lon = range(-180, 181, 60)
    validity = decimal_year_to_mjd2000((2010.0, 2015.0))
    options = {"max_degree": 300}
    def load(self):
        return load_model_emm(EMM_2010_STATIC, EMM_2010_SECVAR)


class TestIGRF11(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((1900.0, 2015.0))
    def load(self):
        return load_model_igrf(IGRF11)


class TestIGRF12(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((1900.0, 2020.0))
    def load(self):
        return load_model_shc(IGRF12)


class TestSIFM(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((2013.4976, 2015.4962))
    def load(self):
        return load_model_shc(SIFM)


class TestCHAOS5Static(TestCase, SHModelTestMixIn):
    validity = (-inf, inf)
    def load(self):
        return load_model_shc(CHAOS5_STATIC)


class TestCHAOS5Core(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((1997.0021, 2015.0007))
    def load(self):
        return load_model_shc(CHAOS5_CORE)


class TestCHAOS5CoreV4(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((1997.1020, 2016.1027))
    def load(self):
        return load_model_shc(CHAOS5_CORE_V4)


class TestCHAOS5Combined(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((1997.1020, 2016.1027))
    def load(self):
        return load_model_shc_combined(CHAOS5_STATIC, CHAOS5_CORE_V4)


class TestCHAOS6Static(TestCase, SHModelTestMixIn):
    validity = (-inf, inf)
    def load(self):
        return load_model_shc(CHAOS6_STATIC)


class TestCHAOS6Core(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((1997.102, 2016.6023))
    def load(self):
        return load_model_shc(CHAOS6_CORE)


class TestCHAOS6CoreX3(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((1997.102, 2017.6016))
    def load(self):
        return load_model_shc(CHAOS6_CORE_X3)


class TestCHAOS6Combined(TestCase, SHModelTestMixIn):
    validity = decimal_year_to_mjd2000((1997.102, 2017.6016))
    def load(self):
        return load_model_shc_combined(CHAOS6_CORE_X3, CHAOS6_STATIC)

#-------------------------------------------------------------------------------

class TestMMA2CInternal(TestCase, DipoleSHModelTestMixIn):
    validity = (6179.125, 6209.875)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]
    def load(self):
        return load_model_swarm_mma_2c_internal(SWARM_MMA_SHA_2C_TEST_DATA)


class TestMMA2CExternal(TestCase, DipoleSHModelTestMixIn):
    validity = (6179.125, 6209.875)
    def load(self):
        return load_model_swarm_mma_2c_external(SWARM_MMA_SHA_2C_TEST_DATA)


class TestMIOSecondary(TestCase, DipoleMIOSHModelTestMixIn):
    validity = (-inf, inf)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]
    def load(self):
        return load_model_swarm_mio_internal(
            SWARM_MIO_SHA_2_TEST_DATA, self.f107
        )


class TestMIOPrimaryAboveIoSph(TestCase, DipoleMIOSHModelTestMixIn):
    validity = (-inf, inf)
    def load(self):
        return load_model_swarm_mio_external(
            SWARM_MIO_SHA_2_TEST_DATA, self.f107, above_ionosphere=True
        )


class TestMIOPrimaryBelowIoSph(TestCase, DipoleMIOSHModelTestMixIn):
    validity = (-inf, inf)
    def load(self):
        return load_model_swarm_mio_external(
            SWARM_MIO_SHA_2_TEST_DATA, self.f107, above_ionosphere=False
        )

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
