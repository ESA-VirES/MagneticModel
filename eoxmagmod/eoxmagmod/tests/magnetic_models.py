#-------------------------------------------------------------------------------
#
#  World Magnetic Model - test
#
# Author: Martin Paces <martin.paces@eox.at>
#
#-------------------------------------------------------------------------------
# Copyright (C) 2014 EOX IT Services GmbH
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies of this Software or works derived from this Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#-------------------------------------------------------------------------------
#pylint: disable=missing-docstring

from unittest import TestCase, main
from itertools import product
from random import uniform
from numpy import array
from numpy.testing import assert_allclose
from eoxmagmod.wmm import read_model_wmm2010, read_model_wmm2015
from eoxmagmod.emm import read_model_emm2010
from eoxmagmod.igrf import read_model_igrf11
from eoxmagmod.shc import read_model_shc
from eoxmagmod.data import (
    CHAOS5_CORE, CHAOS5_CORE_V4, CHAOS5_STATIC,
    CHAOS6_CORE_LATEST, CHAOS6_STATIC,
    IGRF12, SIFM,
)
from eoxmagmod._pywmm import (
    GRADIENT, GEOCENTRIC_SPHERICAL, sheval,
)


class MagneticModelMixIn(object):
    is_internal = True
    validity = None
    model = None

    @property
    def coordinates(self):
        return array([
            (lat, lon, 6371.2*uniform(1.0, 2.0)) for lat, lon
            in product(range(-90, 91, 5), range(-180, 181, 10))
        ])

    @property
    def date(self):
        return uniform(*self.model.validity)

    def eval_reference(self, date, coords):
        degree = self.model.degree_static
        coef_g, coef_h = self.model.get_coef_static(date)
        return sheval(
            coords, mode=GRADIENT, is_internal=self.is_internal,
            degree=degree, coef_g=coef_g, coef_h=coef_h,
            coord_type_in=GEOCENTRIC_SPHERICAL,
            coord_type_out=GEOCENTRIC_SPHERICAL,
            scale_gradient=-1.0,
        )

    def eval_model(self, date, coords):
        return self.model.eval(
            coords, date,
            coord_type_in=GEOCENTRIC_SPHERICAL,
            coord_type_out=GEOCENTRIC_SPHERICAL,
        )

    def test_validity(self):
        self.assertEqual(self.model.validity, self.validity)

    def test_eval(self):
        date = self.date
        coords = self.coordinates
        assert_allclose(
            self.eval_model(date, coords),
            self.eval_reference(date, coords),
        )

#-------------------------------------------------------------------------------

class TestEMM2010(TestCase, MagneticModelMixIn):
    validity = (2010.0, 2015.0)
    model = read_model_emm2010()


class TestWMM2010(TestCase, MagneticModelMixIn):
    validity = (2010.0, 2015.0)
    model = read_model_wmm2010()


class TestWMM2015(TestCase, MagneticModelMixIn):
    validity = (2015.0, 2020.0)
    model = read_model_wmm2015()


class TestIGRF11(TestCase, MagneticModelMixIn):
    validity = (1900.0, 2015.0)
    model = read_model_igrf11()


class TestIGRF12(TestCase, MagneticModelMixIn):
    validity = (1900.0, 2020.0)
    model = read_model_shc(IGRF12)


class TestSIFM(TestCase, MagneticModelMixIn):
    validity = (2013.4976, 2015.4962)
    model = read_model_shc(SIFM)


class TestCHAOS5Static(TestCase, MagneticModelMixIn):
    validity = (2005.0021, 2005.0021)
    model = read_model_shc(CHAOS5_STATIC)


class TestCHAOS5Core(TestCase, MagneticModelMixIn):
    validity = (1997.0021, 2015.0007)
    model = read_model_shc(CHAOS5_CORE)


class TestCHAOS5CoreV4(TestCase, MagneticModelMixIn):
    validity = (1997.1020, 2016.1027)
    model = read_model_shc(CHAOS5_CORE_V4)


class TestCHAOS6Static(TestCase, MagneticModelMixIn):
    validity = (2005.0021, 2005.0021)
    model = read_model_shc(CHAOS6_STATIC)


class TestCHAOS6Core(TestCase, MagneticModelMixIn):
    validity = (1997.102, 2018.1013)
    model = read_model_shc(CHAOS6_CORE_LATEST)

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
