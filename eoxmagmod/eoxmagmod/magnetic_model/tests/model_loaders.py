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
from numpy import nan, inf, isinf, array, empty, full, nditer, asarray
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
    load_model_swarm_mma_2f_geo_internal,
    load_model_swarm_mma_2f_geo_external,
    load_model_swarm_mma_2f_sm_internal,
    load_model_swarm_mma_2f_sm_external,
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
    SWARM_MMA_SHA_2F_TEST_DATA,
    SWARM_MIO_SHA_2_TEST_DATA,
)
from eoxmagmod.magnetic_model.model import (
    SphericalHarmomicGeomagneticModel,
    DipoleSphericalHarmomicGeomagneticModel,
)
from eoxmagmod.magnetic_model.model_mio import (
    DipoleMIOPrimaryGeomagneticModel,
    DipoleMIOGeomagneticModel,
)
from eoxmagmod._pywmm import GEOCENTRIC_SPHERICAL, sheval, GRADIENT
from eoxmagmod.sheval_dipole import sheval_dipole


class SHModelTestMixIn(object):
    parameters = ("time", "location")
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
    def time_before(self):
        return self.validity[0] - 1.0

    @property
    def time_after(self):
        return self.validity[1] + 1.0

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
            if self._is_valid_time(time):
                vect0[...], vect1[...], vect2[...] = self._eval_reference(
                    time, [coord0, coord1, coord2]
                )
            else:
                vect0[...], vect1[...], vect2[...] = nan, nan, nan
        return result

    def _is_valid_time(self, time):
        validity_start, validity_end = self.model.validity
        return validity_start <= time <= validity_end

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

    def test_parameters(self):
        self.assertEqual(self.model.parameters, self.parameters)

    def test_eval_single_time(self):
        time = self.time
        coords = self.coordinates
        assert_allclose(
            self.eval_model(time, coords),
            self.eval_reference(time, coords),
        )

    def _test_eval_single_time_invalid(self, time):
        coords = self.coordinates
        assert_allclose(
            self.eval_model(time, coords),
            full(coords.shape, nan)
        )

    def test_eval_single_time_nan(self):
        self._test_eval_single_time_invalid(nan)

    def test_eval_single_time_before(self):
        time = self.time_before
        if not isinf(time):
            self._test_eval_single_time_invalid(time)

    def test_eval_single_time_after(self):
        time = self.time_after
        if not isinf(time):
            self._test_eval_single_time_invalid(time)

    def test_eval_multi_time(self):
        times = self.times
        coords = self.coordinates
        assert_allclose(
            self.eval_model(times, coords),
            self.eval_reference(times, coords),
        )

    def test_eval_multi_time_invalid(self):
        times = array([
            time for time in [nan, self.time_before, self.time, self.time_after]
            if not isinf(time)
        ])
        coords = array([(0, 0, 6371.2) for _ in times])
        assert_allclose(
            self.eval_model(times, coords),
            self.eval_reference(times, coords),
        )

    def test_eval_reference_values(self):
        times, coords, results = self.reference_values
        assert_allclose(self.eval_model(times, coords), results)

    def test_eval_empty_coords(self):
        assert_allclose(
            self.eval_model(self.time, empty((0, 3))),
            empty((0, 3))
        )

    def test_eval_empty_time_and_coords(self):
        assert_allclose(
            self.eval_model(empty(0), empty((0, 3))),
            empty((0, 3))
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
    parameters = ("time", "location", "f107", "subsolar_point")
    f107 = 70.0
    model_class = DipoleMIOGeomagneticModel

    def eval_model(self, times, coords):
        return self.model.eval(times, coords, f107=self.f107, **self.options)

    def _eval_reference(self, time, coords):
        return self._eval_reference_mio(self.model, time, coords)

    def _eval_reference_mio(self, model, time, coords):
        is_internal = model.coefficients.is_internal
        scale = -(
            1.0 + self.f107 * model.wolf_ratio
        ) * asarray(self.scale)
        lat_ngp, lon_ngp = model.north_pole(time)
        coeff, degree = model.coefficients(time, lat_ngp, lon_ngp)
        return sheval_dipole(
            coords, degree, coeff[..., 0], coeff[..., 1], lat_ngp, lon_ngp,
            is_internal=is_internal, mode=GRADIENT,
            coord_type_in=GEOCENTRIC_SPHERICAL,
            coord_type_out=GEOCENTRIC_SPHERICAL,
            scale_gradient=scale
        )

#-------------------------------------------------------------------------------

class TestWMM2010(TestCase, SHModelTestMixIn):
    reference_values = (
        4566.0, (30.0, 40.0, 8000.0),
        (15123.605974201277, 431.1067254253052, 14617.02644010297)
    )
    validity = decimal_year_to_mjd2000((2010.0, 2015.0))
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_wmm(WMM_2010)


class TestWMM2015(TestCase, SHModelTestMixIn):
    reference_values = (
        6392.0, (30.0, 40.0, 8000.0),
        (15124.716592471135, 533.1027023540182, -14728.4938691708)
    )
    validity = decimal_year_to_mjd2000((2015.0, 2020.0))

    def load(self):
        return load_model_wmm(WMM_2015)


class TestEMM2010(TestCase, SHModelTestMixIn):
    reference_values = (
        4566.0, (30.0, 40.0, 8000.0),
        (15124.606019372684, 442.2376840179962, -14612.282120230499)
    )
    # The EMM models is huge and the test ranges have to be reduced.
    range_lat = range(-90, 91, 30)
    range_lon = range(-180, 181, 60)
    validity = decimal_year_to_mjd2000((2010.0, 2015.0))
    options = {"max_degree": 300}

    def load(self):
        return load_model_emm(EMM_2010_STATIC, EMM_2010_SECVAR)


class TestIGRF11(TestCase, SHModelTestMixIn):
    reference_values = (
        -14609.5, (30.0, 40.0, 8000.0),
        (15265.918081037888, -142.6442876878355, -14044.282413158882)
    )
    validity = decimal_year_to_mjd2000((1900.0, 2015.0))

    def load(self):
        return load_model_igrf(IGRF11)


class TestIGRF12(TestCase, SHModelTestMixIn):
    reference_values = (
        -15522.5, (30.0, 40.0, 8000.0),
        (15259.57386772841, -159.00767967612023, -14015.952721753336)
    )
    validity = decimal_year_to_mjd2000((1900.0, 2020.0))

    def load(self):
        return load_model_shc(IGRF12)


class TestSIFM(TestCase, SHModelTestMixIn):
    reference_values = (
        5295.36, (30.0, 40.0, 8000.0),
        (15122.448070753977, 474.14615304317635, -14669.16289251053)
    )
    validity = decimal_year_to_mjd2000((2013.4976, 2015.4962))

    def load(self):
        return load_model_shc(SIFM)


class TestCHAOS5Static(TestCase, SHModelTestMixIn):
    reference_values = (
        0.0, (30.0, 40.0, 8000.0),
        (-0.019165363389425448, 0.017766807977599153, 0.007245125734944849)
    )
    validity = (-inf, inf)

    def load(self):
        return load_model_shc(CHAOS5_STATIC)


class TestCHAOS5Core(TestCase, SHModelTestMixIn):
    reference_values = (
        2192.51, (30.0, 40.0, 8000.0),
        (15126.610410635147, 302.75469100239826, -14477.55985460029)
    )
    validity = decimal_year_to_mjd2000((1997.0021, 2015.0007))

    def load(self):
        return load_model_shc(CHAOS5_CORE)


class TestCHAOS5CoreV4(TestCase, SHModelTestMixIn):
    reference_values = (
        2411.9, (30.0, 40.0, 8000.0),
        (15127.03768745214, 313.61814829613326, -14489.207459734534)
    )
    validity = decimal_year_to_mjd2000((1997.1020, 2016.1027))

    def load(self):
        return load_model_shc(CHAOS5_CORE_V4)


class TestCHAOS5Combined(TestCase, SHModelTestMixIn):
    reference_values = (
        2411.9, (30.0, 40.0, 8000.0),
        (15127.018522088763, 313.6359151041108, -14489.200214608843)
    )
    validity = decimal_year_to_mjd2000((1997.1020, 2016.1027))

    def load(self):
        return load_model_shc_combined(CHAOS5_STATIC, CHAOS5_CORE_V4)


class TestCHAOS6Static(TestCase, SHModelTestMixIn):
    reference_values = (
        0.0, (30.0, 40.0, 8000.0),
        (-0.006745769467490476, 0.00860457221837856, -0.010495388357779979)
    )
    validity = (-inf, inf)

    def load(self):
        return load_model_shc(CHAOS6_STATIC)


class TestCHAOS6Core(TestCase, SHModelTestMixIn):
    reference_values = (
        2503.33, (30.0, 40.0, 8000.0),
        (15127.146281343608, 318.51792709726175, -14493.952978715943)
    )
    validity = decimal_year_to_mjd2000((1997.102, 2016.6023))

    def load(self):
        return load_model_shc(CHAOS6_CORE)


class TestCHAOS6CoreX3(TestCase, SHModelTestMixIn):
    reference_values = (
        2685.9, (30.0, 40.0, 8000.0),
        (15127.196090133599, 328.5862052582883, -14503.664172833218)
    )
    validity = decimal_year_to_mjd2000((1997.102, 2017.6016))

    def load(self):
        return load_model_shc(CHAOS6_CORE_X3)


class TestCHAOS6Combined(TestCase, SHModelTestMixIn):
    reference_values = (
        2685.9, (30.0, 40.0, 8000.0),
        (15127.189344364127, 328.594809830505, -14503.674668221536)
    )
    validity = decimal_year_to_mjd2000((1997.102, 2017.6016))

    def load(self):
        return load_model_shc_combined(CHAOS6_CORE_X3, CHAOS6_STATIC)

#-------------------------------------------------------------------------------

class TestMMA2CSecondary(TestCase, DipoleSHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (1.7252467863888683, 0.27791273383414994, -0.12422361564742368)
    )
    validity = (6179.125, 6209.875)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mma_2c_internal(SWARM_MMA_SHA_2C_TEST_DATA)


class TestMMA2CPrimary(TestCase, DipoleSHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (-7.474051407972587, 3.531499380152684, -4.628812102394507)
    )
    validity = (6179.125, 6209.875)

    def load(self):
        return load_model_swarm_mma_2c_external(SWARM_MMA_SHA_2C_TEST_DATA)


class TestMMA2FGeoSecondary(TestCase, SHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (1.7678502698292433, 0.6267115585524842, 2.7484695371405405)
    )
    validity = (6179.03125, 6209.96875)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mma_2f_geo_internal(SWARM_MMA_SHA_2F_TEST_DATA)


class TestMMA2FGeoPrimary(TestCase, SHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (-9.114015792291584, 6.856282080637684, 3.208391426427198)
    )
    validity = (6179.03125, 6209.96875)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mma_2f_geo_external(SWARM_MMA_SHA_2F_TEST_DATA)


class TestMMA2FSMSecondary(TestCase, SHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (1.6186505587469782, 1.0283338998596887, 2.6779138076728497)
    )
    validity = (6179.03125, 6209.96875)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mma_2f_sm_internal(SWARM_MMA_SHA_2F_TEST_DATA)


class TestMMA2FSMPrimary(TestCase, SHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (-9.42096313502333, 5.586931375284516, 4.7449677343745975)
    )
    validity = (6179.03125, 6209.96875)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mma_2f_sm_external(SWARM_MMA_SHA_2F_TEST_DATA)

#-------------------------------------------------------------------------------

class TestMIOSecondary(TestCase, DipoleMIOSHModelTestMixIn):
    reference_values = (
        5661.87, (30.0, 40.0, 8000.0),
        (-0.5388282699123806, -0.17622120922727555, -1.6137152691151841)
    )
    validity = (-inf, inf)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mio_internal(SWARM_MIO_SHA_2_TEST_DATA)


class TestMIOPrimary(TestCase, DipoleMIOSHModelTestMixIn):
    model_class = DipoleMIOPrimaryGeomagneticModel
    reference_values = (
        5661.87,
        [
            (30.0, 40.0, 6400.0), # below ionosphere r < (a + h)
            (30.0, 40.0, 8000.0), # above ionosphere r > (a + h)
        ],
        [
            (-0.6061225119866813, -0.6088386296175435, -4.733769204526618),
            (0.2356719922628632, 0.19030444647263053, -1.9489199024730584),
        ]
    )
    validity = (-inf, inf)

    def load(self):
        return load_model_swarm_mio_external(SWARM_MIO_SHA_2_TEST_DATA)

    def _eval_reference(self, time, coords):
        height_radius = self.model.earth_radius + self.model.height

        if coords[2] <= height_radius:
            model = self.model.model_below_ionosphere
        else:
            model = self.model.model_above_ionosphere

        return self._eval_reference_mio(model, time, coords)

    def test_eval_single_reference_value_below_ionosphere(self):
        times, coords, results = self.reference_values
        assert_allclose(self.eval_model(times, coords[0]), results[0])

    def test_eval_single_reference_value_above_ionosphere(self):
        times, coords, results = self.reference_values
        assert_allclose(self.eval_model(times, coords[1]), results[1])


class TestMIOPrimaryAboveIonosphere(TestCase, DipoleMIOSHModelTestMixIn):
    reference_values = (
        5661.87,
        (30.0, 40.0, 8000.0),
        (0.2356719922628632, 0.19030444647263053, -1.9489199024730584)
    )
    validity = (-inf, inf)

    def load(self):
        return load_model_swarm_mio_external(
            SWARM_MIO_SHA_2_TEST_DATA, above_ionosphere=True
        )


class TestMIOPrimaryBelowIonosphere(TestCase, DipoleMIOSHModelTestMixIn):
    reference_values = (
        5661.87, (30.0, 40.0, 6400.0),
        (-0.6061225119866813, -0.6088386296175435, -4.733769204526618)
    )
    validity = (-inf, inf)

    def load(self):
        return load_model_swarm_mio_external(
            SWARM_MIO_SHA_2_TEST_DATA, above_ionosphere=False
        )

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
