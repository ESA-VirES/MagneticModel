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
# pylint: disable=missing-docstring,no-self-use,invalid-name

from __future__ import print_function
from unittest import TestCase, main
from itertools import product
from numpy import nan, inf, isinf, array, empty, full, nditer, asarray
from numpy.random import uniform
from numpy.testing import assert_allclose
from eoxmagmod.magnetic_time import mjd2000_to_magnetic_universal_time
from eoxmagmod.time_util import (
    decimal_year_to_mjd2000, decimal_year_to_mjd2000_simple,
)
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
    CHAOS6_CORE_LATEST, CHAOS6_STATIC,
    IGRF11, IGRF12, SIFM,
)
from eoxmagmod.magnetic_model.tests.data import (
    SWARM_MMA_SHA_2C_TEST_DATA,
    SWARM_MMA_SHA_2F_TEST_DATA,
    SWARM_MIO_SHA_2_TEST_DATA,
    CHAOS_MMA_TEST_DATA,
)
from eoxmagmod.magnetic_model.model import (
    SphericalHarmomicGeomagneticModel,
    DipoleSphericalHarmomicGeomagneticModel,
)
from eoxmagmod.magnetic_model.model_mio import (
    DipoleMIOPrimaryGeomagneticModel,
    DipoleMIOGeomagneticModel,
)
from eoxmagmod.magnetic_model.model_composed import (
    ComposedGeomagneticModel,
)
from eoxmagmod._pymm import (
    GEOCENTRIC_CARTESIAN, GEODETIC_ABOVE_WGS84, GEOCENTRIC_SPHERICAL, convert,
    GRADIENT, sheval,
)
from eoxmagmod.sheval_dipole import sheval_dipole


class SHModelTestMixIn(object):
    coord_type_in = GEOCENTRIC_SPHERICAL
    coord_type_out = GEOCENTRIC_SPHERICAL
    parameters = ("time", "location")
    scale = 1.0
    range_lat = range(-90, 91, 5)
    range_lon = range(-180, 181, 10)
    validity = None
    degree = None
    model_class = SphericalHarmomicGeomagneticModel
    options = {}

    @property
    def model(self):
        if not hasattr(self, "_model"):
            self._model = self.load()
        return self._model

    @property
    def coordinates(self):
        return convert(array([
            (lat, lon, 6371.2*uniform(1.0, 2.0)) for lat, lon
            in product(self.range_lat, self.range_lon)
        ]), GEOCENTRIC_SPHERICAL, self.coord_type_in)

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
        return self.model.eval(
            times, coords, input_coordinate_system=self.coord_type_in,
            output_coordinate_system=self.coord_type_out, **self.options
        )

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
            coord_type_in=self.coord_type_in,
            coord_type_out=self.coord_type_out,
            scale_gradient=-asarray(self.scale),
        )

    def test_class(self):
        self.assertIsInstance(self.model, self.model_class)

    def test_degree(self):
        if hasattr(self.model, 'degree'):
            self.assertEqual(self.model.degree, self.degree)

    def test_validity(self):
        assert_allclose(self.model.validity, self.validity)

    def test_parameters(self):
        self.assertEqual(set(self.model.parameters), set(self.parameters))

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
        coords = convert(
            array([(0, 0, 6371.2) for _ in times]),
            GEOCENTRIC_SPHERICAL, self.coord_type_in
        )
        assert_allclose(
            self.eval_model(times, coords),
            self.eval_reference(times, coords),
        )

    def test_eval_reference_values(self):
        times, src_coords, expected_result = self.reference_values
        coords = convert(src_coords, GEOCENTRIC_SPHERICAL, self.coord_type_in)
        model_result = self.eval_model(times, coords)
        try:
            assert_allclose(model_result, expected_result)
        except AssertionError:
            print(tuple(float(f) for f in model_result))
            raise

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
        lat_ngp, lon_ngp = self.model.north_pole #(time)
        coeff, degree = self.model.coefficients(time)
        return sheval_dipole(
            coords, degree, coeff[..., 0], coeff[..., 1], lat_ngp, lon_ngp,
            is_internal=is_internal, mode=GRADIENT,
            coord_type_in=self.coord_type_in,
            coord_type_out=self.coord_type_out,
            scale_gradient=-asarray(self.scale),
        )


class DipoleMIOSHModelTestMixIn(SHModelTestMixIn):
    parameters = ("time", "location", "f107", "subsolar_point")
    f107 = 70.0
    model_class = DipoleMIOGeomagneticModel

    def eval_model(self, times, coords):
        return self.model.eval(
            times, coords, f107=self.f107,
            input_coordinate_system=self.coord_type_in,
            output_coordinate_system=self.coord_type_out, **self.options
        )

    def _eval_reference(self, time, coords):
        return self._eval_reference_mio(self.model, time, coords)

    def _eval_reference_mio(self, model, time, coords):
        is_internal = model.coefficients.is_internal
        scale = -(
            1.0 + self.f107 * model.wolf_ratio
        ) * asarray(self.scale)
        lat_ngp, lon_ngp = model.north_pole #(time)
        coeff, degree = model.coefficients(
            time, mjd2000_to_magnetic_universal_time(time, lat_ngp, lon_ngp)
        )
        return sheval_dipole(
            coords, degree, coeff[..., 0], coeff[..., 1], lat_ngp, lon_ngp,
            is_internal=is_internal, mode=GRADIENT,
            coord_type_in=self.coord_type_in,
            coord_type_out=self.coord_type_out,
            scale_gradient=scale
        )

class ComposedModelTestMixIn(SHModelTestMixIn):
    components = None
    model_class = ComposedGeomagneticModel

    def load(self):
        composed_model = ComposedGeomagneticModel()
        for model, scale, parameters in self.components:
            composed_model.push(model, scale, **parameters)
        return composed_model

    def _eval_reference(self, time, coords):
        result = 0
        for model, scale, parameters in self.components:
            kwargs = {}
            kwargs.update(getattr(self, 'options', {}))
            kwargs.update(parameters)
            kwargs.update({
                "input_coordinate_system": self.coord_type_in,
                "output_coordinate_system": self.coord_type_out,
            })
            result += scale * model.eval(time, coords, **kwargs)
        return result

#-------------------------------------------------------------------------------

class TestComposedModelFull(TestCase, ComposedModelTestMixIn):
    parameters = ("time", "location", "f107", "subsolar_point")
    options = {"f107": 70, "scale": [1, 1, -1]}
    components = [
        (
            load_model_shc_combined(
                CHAOS6_STATIC, CHAOS6_CORE_LATEST,
                to_mjd2000=decimal_year_to_mjd2000_simple,
            ), 1.0, {}
        ),
        (load_model_swarm_mma_2c_internal(CHAOS_MMA_TEST_DATA), 1.0, {}),
        (load_model_swarm_mma_2c_external(CHAOS_MMA_TEST_DATA), 1.0, {}),
        (load_model_swarm_mio_internal(SWARM_MIO_SHA_2_TEST_DATA), 1.0, {}),
        (load_model_swarm_mio_external(SWARM_MIO_SHA_2_TEST_DATA), 1.0, {}),
    ]
    reference_values = (
        6201.125, (30.0, 40.0, 6400.0), # below ionosphere r < (a + h)
        (30289.26960669176, 2254.2152561714115, 31788.9069112966)
    )
    validity = (6179.00000, 6209.979167)


class TestComposedModelFullCartToWGS84(TestComposedModelFull):
    coord_type_in = GEOCENTRIC_CARTESIAN
    coord_type_out = GEODETIC_ABOVE_WGS84
    reference_values = (
        6201.125, (30.0, 40.0, 6400.0), # below ionosphere r < (a + h)
        (30381.359560356446, 2254.2152561714115, 31700.90609408948)
    )


class TestComposedModelFullWGS84ToCart(TestComposedModelFull):
    coord_type_in = GEODETIC_ABOVE_WGS84
    coord_type_out = GEOCENTRIC_CARTESIAN
    reference_values = (
        6201.125, (30.0, 40.0, 6400.0), # below ionosphere r < (a + h)
        (-34139.649212401426, -25703.898035517483, -10336.823485822737)
    )


class TestComposedModelDiffConstrained(TestCase, ComposedModelTestMixIn):
    parameters = ("time", "location")
    options = {"scale": [1, 1, -1]}
    components = [
        (
            load_model_shc_combined(
                CHAOS6_STATIC, CHAOS6_CORE_LATEST,
                to_mjd2000=decimal_year_to_mjd2000_simple
            ), 1.0, {"max_degree": 15}
        ),
        (
            load_model_shc(
                CHAOS6_CORE_LATEST,
                to_mjd2000=decimal_year_to_mjd2000_simple
            ), -1.0, {"min_degree": 5, "max_degree": 15}
        ),
    ]
    reference_values = (
        6201.125,
        (30.0, 40.0, 6400.0),
        (31513.116189124143, 2660.4562533211647, 30567.71751831368)
    )
    validity = (-1058.4945, 6976.49415)


class TestComposedModelDiffConstrainedCartToWGS84(TestComposedModelDiffConstrained):
    coord_type_in = GEOCENTRIC_CARTESIAN
    coord_type_out = GEODETIC_ABOVE_WGS84
    reference_values = (
        6201.125, (30.0, 40.0, 6400.0), # below ionosphere r < (a + h)
        (31601.658407935865, 2660.4562533211647, 30476.17154592741)
    )


class TestComposedModelDiffConstrainedWGS84ToCart(TestComposedModelDiffConstrained):
    coord_type_in = GEODETIC_ABOVE_WGS84
    coord_type_out = GEOCENTRIC_CARTESIAN
    reference_values = (
        6201.125, (30.0, 40.0, 6400.0), # below ionosphere r < (a + h)
        (-34059.382252106414, -25106.236099833302, -12007.300413035513)
    )


# TODO coordinate system conversions

#-------------------------------------------------------------------------------

class TestWMM2010(TestCase, SHModelTestMixIn):
    reference_values = (
        4566.0, (30.0, 40.0, 8000.0),
        (15123.605974201277, 431.1067254253052, 14617.02644010297)
    )
    degree = 12
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
    degree = 12
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
    degree = 739
    validity = decimal_year_to_mjd2000((2010.0, 2015.0))
    options = {"max_degree": 300}

    def load(self):
        return load_model_emm(EMM_2010_STATIC, EMM_2010_SECVAR)


class TestIGRF11(TestCase, SHModelTestMixIn):
    reference_values = (
        -14609.5, (30.0, 40.0, 8000.0),
        (15265.918081037888, -142.6442876878355, -14044.282413158882)
    )
    degree = 13
    validity = decimal_year_to_mjd2000((1900.0, 2015.0))

    def load(self):
        return load_model_igrf(IGRF11)


class TestIGRF12(TestCase, SHModelTestMixIn):
    reference_values = (
        -15522.5, (30.0, 40.0, 8000.0),
        (15259.57386772841, -159.00767967612023, -14015.952721753336)
    )
    degree = 13
    validity = decimal_year_to_mjd2000((1900.0, 2020.0))

    def load(self):
        return load_model_shc(IGRF12, interpolate_in_decimal_years=True)


class TestSIFM(TestCase, SHModelTestMixIn):
    reference_values = (
        5295.36, (30.0, 40.0, 8000.0),
        (15122.448070753977, 474.14615304317635, -14669.16289251053)
    )
    degree = 70
    validity = decimal_year_to_mjd2000((2013.4976, 2015.4962))

    def load(self):
        return load_model_shc(SIFM)


class TestCHAOS5Static(TestCase, SHModelTestMixIn):
    reference_values = (
        0.0, (30.0, 40.0, 8000.0),
        (-0.019165363389425448, 0.017766807977599153, 0.007245125734944849)
    )
    degree = 90
    validity = (-inf, inf)

    def load(self):
        return load_model_shc(CHAOS5_STATIC)


class TestCHAOS5Core(TestCase, SHModelTestMixIn):
    reference_values = (
        2192.51, (30.0, 40.0, 8000.0),
        (15126.611217467429, 302.7784261453687, -14477.586706907041)
    )
    degree = 20
    validity = decimal_year_to_mjd2000_simple((1997.0021, 2015.0007))

    def load(self):
        return load_model_shc(
            CHAOS5_CORE, to_mjd2000=decimal_year_to_mjd2000_simple,
        )


class TestCHAOS5CoreWithOverridenValidity(TestCHAOS5Core):
    validity = decimal_year_to_mjd2000_simple((1999.0, 2015.0))

    def load(self):
        return load_model_shc(
            CHAOS5_CORE,
            to_mjd2000=decimal_year_to_mjd2000_simple,
            validity_start=1999.0, validity_end=2015.0
        )


class TestCHAOS5CoreV4(TestCase, SHModelTestMixIn):
    reference_values = (
        2411.9, (30.0, 40.0, 8000.0),
        (15127.03768745214, 313.61814829613326, -14489.207459734534)
    )
    degree = 20
    validity = decimal_year_to_mjd2000((1997.1020, 2016.1027))

    def load(self):
        return load_model_shc(CHAOS5_CORE_V4)


class TestCHAOS5Combined(TestCase, SHModelTestMixIn):
    reference_values = (
        2411.9, (30.0, 40.0, 8000.0),
        (15127.018861793757, 313.6542615062537, -14489.218457022034)
    )
    degree = 90
    validity = decimal_year_to_mjd2000_simple((1997.1020, 2016.1027))

    def load(self):
        return load_model_shc_combined(
            CHAOS5_STATIC, CHAOS5_CORE_V4,
            to_mjd2000=decimal_year_to_mjd2000_simple,
        )


class TestCHAOS5CombinedOverridenValidity(TestCHAOS5Combined):
    validity = decimal_year_to_mjd2000_simple((1999.0, 2015.0))

    def load(self):
        return load_model_shc_combined(
            CHAOS5_STATIC, CHAOS5_CORE_V4,
            to_mjd2000=decimal_year_to_mjd2000_simple,
            validity_start=1999.0, validity_end=2015.0
        )


class TestCHAOS6Static(TestCase, SHModelTestMixIn):
    reference_values = (
        0.0, (30.0, 40.0, 8000.0),
        (-0.006745769467490476, 0.00860457221837856, -0.010495388357779979)
    )
    degree = 110
    validity = (-inf, inf)

    def load(self):
        return load_model_shc(CHAOS6_STATIC)


class TestCHAOS6Core(TestCase, SHModelTestMixIn):
    reference_values = (
        2503.33, (30.0, 40.0, 8000.0),
        (15127.112741712719, 318.52904189382974, -14493.826457120203)
    )
    degree = 20
    validity = decimal_year_to_mjd2000((1997.102, 2019.1006))

    def load(self):
        return load_model_shc(CHAOS6_CORE_LATEST)


class TestCHAOS6CoreWithOverridenValidity(TestCHAOS6Core):
    validity = decimal_year_to_mjd2000((2000.0, 2018.0))

    def load(self):
        return load_model_shc(
            CHAOS6_CORE_LATEST,
            validity_start=2000.0,
            validity_end=2018.0
        )


class TestCHAOS6Combined(TestCase, SHModelTestMixIn):
    reference_values = (
        2685.9, (30.0, 40.0, 8000.0),
        (15127.164745882454, 328.584740011225, -14503.5849807797)
    )
    degree = 110
    validity = decimal_year_to_mjd2000((1997.102, 2019.1006))

    def load(self):
        return load_model_shc_combined(CHAOS6_CORE_LATEST, CHAOS6_STATIC)


class TestCHAOS6CombinedOverridenValidity(TestCHAOS6Combined):
    validity = decimal_year_to_mjd2000((2000.0, 2018.0))

    def load(self):
        return load_model_shc_combined(
            CHAOS6_CORE_LATEST, CHAOS6_STATIC,
            validity_start=2000.0, validity_end=2018.0
        )

#-------------------------------------------------------------------------------

class TestMMA2CSecondary(TestCase, DipoleSHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (1.7252467863888683, 0.27791273383414994, -0.12422361564742368)
    )
    degree = 3
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
    degree = 2
    validity = (6179.125, 6209.875)

    def load(self):
        return load_model_swarm_mma_2c_external(SWARM_MMA_SHA_2C_TEST_DATA)


class TestChaosMMASecondary(TestCase, DipoleSHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (1.8492638163980442, 0.5125018012040559, 1.0821299594918217)
    )
    degree = 2
    validity = (6179.00000, 6209.979167)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mma_2c_internal(CHAOS_MMA_TEST_DATA)


class TestChaosMMAPrimary(TestCase, DipoleSHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (-8.667405753073385, 4.538967766836233, 6.576263454698334)
    )
    degree = 2
    validity = (6179.00000, 6209.979167)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mma_2c_external(CHAOS_MMA_TEST_DATA)


class TestMMA2FGeoSecondary(TestCase, SHModelTestMixIn):
    reference_values = (
        6194.5, (30.0, 40.0, 8000.0),
        (1.7678502698292433, 0.6267115585524842, 2.7484695371405405)
    )
    degree = 1
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
    degree = 1
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
    degree = 1
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
    degree = 1
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
    degree = 2
    validity = (-inf, inf)
    options = {"scale": [1, 1, -1]}
    scale = [1, 1, -1]

    def load(self):
        return load_model_swarm_mio_internal(SWARM_MIO_SHA_2_TEST_DATA)


class TestMIOPrimary(TestCase, DipoleMIOSHModelTestMixIn):
    model_class = DipoleMIOPrimaryGeomagneticModel
    degree = 2
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
    degree = 2
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
    degree = 2
    validity = (-inf, inf)

    def load(self):
        return load_model_swarm_mio_external(
            SWARM_MIO_SHA_2_TEST_DATA, above_ionosphere=False
        )

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
