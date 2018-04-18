#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Coefficients.
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

from unittest import TestCase, main
from numpy import inf, array
from numpy.testing import assert_allclose
from eoxmagmod.magnetic_model.coefficients import (
    SparseSHCoefficients,
    SparseSHCoefficientsTimeDependent,
    SparseSHCoefficientsConstant,
    CombinedSHCoefficients,
)

class SHCoefficinetTestMixIn(object):

    def test_degree(self):
        self.assertEqual(self.coefficients.degree, self.degree)

    def test_validity(self):
        self.assertEqual(self.coefficients.validity, self.validity)

    def test_is_internal(self):
        self.assertEqual(self.coefficients.is_internal, self.is_internal)

#-------------------------------------------------------------------------------

class CombinedSHCoefficientsMixIn(object):
    is_internal = True
    times0 = array([2012.0, 2016.0, 2014.0])
    indices0 = array([(1, 0), (1, 1), (1, -1)])
    coeff0 = array([
        [1, 3, 2],
        [5, 15, 10],
        [10, 30, 20],
    ])
    indices1 = array([(2, 0), (2, 1), (2, -1), (2, 2), (2, -2)])
    coeff1 = array([1, 5, 10, 8, 12])
    options0 = {}
    options1 = {}
    degree = 2
    validity = (2012.0, 2016.0)

    @property
    def coefficients(self):
        return CombinedSHCoefficients(
            SparseSHCoefficientsTimeDependent(
                self.indices0, self.coeff0, self.times0, **self.options0
            ),
            SparseSHCoefficientsConstant(
                self.indices1, self.coeff1, **self.options1
            )
        )


class TestCombinedSHCoefficientsDefault(TestCase, SHCoefficinetTestMixIn, CombinedSHCoefficientsMixIn):
    def test_callable(self):
        coeff, degree = self.coefficients(2013.0)
        assert_allclose(coeff, [
            [0., 0.], [1.5, 0], [7.5, 15.0], [1, 0], [5, 10.0], [8., 12.],
        ])
        self.assertEqual(degree, self.degree)


class TestCombinedSHCoefficientsInternal(TestCombinedSHCoefficientsDefault):
    is_internal = True
    options0 = {"is_internal": True}
    options1 = {"is_internal": True}


class TestCombinedSHCoefficientsExternal(TestCombinedSHCoefficientsDefault):
    is_internal = False
    options0 = {"is_internal": False}
    options1 = {"is_internal": False}


class TestCombinedSHCoefficientsMixed(TestCase, CombinedSHCoefficientsMixIn):
    options0 = {"is_internal": True}
    options1 = {"is_internal": False}

    def test_mixed_type_failure(self):
        with self.assertRaises(ValueError):
            self.coefficients

#-------------------------------------------------------------------------------

class TestSparseSHCoefficientsConstantDefault(TestCase, SHCoefficinetTestMixIn):
    indices = array([(1, 0), (1, 1), (1, -1)])
    coeff = array([1, 5, 10])
    degree = 1
    is_internal = True
    validity = (-inf, inf)
    options = {}

    @property
    def coefficients(self):
        return SparseSHCoefficientsConstant(
            self.indices, self.coeff, **self.options
        )

    def test_callable(self):
        coeff, degree = self.coefficients(2013.0)
        assert_allclose(coeff, [[0., 0.], [1, 0], [5, 10.0]])
        self.assertEqual(degree, self.degree)


class TestSparseSHCoefficientsConstantInternal(TestSparseSHCoefficientsConstantDefault):
    is_internal = True
    options = {"is_internal": True}


class TestSparseSHCoefficientsConstantExternal(TestSparseSHCoefficientsConstantDefault):
    is_internal = False
    options = {"is_internal": False}


class TestSparseSHCoefficientsConstantExtendedValidity(TestSparseSHCoefficientsConstantDefault):
    validity = (2010.0, 2018.0)
    options = {"validity_start": 2010.0, "validity_end": 2018.0}

#-------------------------------------------------------------------------------

class TestSparseSHCoefficientsTimeDependentDefault(TestCase, SHCoefficinetTestMixIn):
    times = array([2012.0, 2016.0, 2014.0])
    indices = array([(1, 0), (1, 1), (1, -1)])
    coeff = array([
        [1, 3, 2],
        [5, 15, 10],
        [10, 30, 20],
    ])
    degree = 1
    is_internal = True
    validity = (2012.0, 2016.0)
    options = {}

    @property
    def coefficients(self):
        return SparseSHCoefficientsTimeDependent(
            self.indices, self.coeff, self.times, **self.options
        )

    def test_callable(self):
        coeff, degree = self.coefficients(2013.0)
        assert_allclose(coeff, [[0., 0.], [1.5, 0], [7.5, 15.0]])
        self.assertEqual(degree, self.degree)

    def test_callable_before_first_time(self):
        coeff, degree = self.coefficients(2011.0)
        assert_allclose(coeff, [[0., 0.], [0.5, 0], [2.5, 5.0]])
        self.assertEqual(degree, self.degree)


class TestSparseSHCoefficientsTimeDependentInternal(TestSparseSHCoefficientsTimeDependentDefault):
    is_internal = True
    options = {"is_internal": True}


class TestSparseSHCoefficientsTimeDependentExternal(TestSparseSHCoefficientsTimeDependentDefault):
    is_internal = False
    options = {"is_internal": False}


class TestSparseSHCoefficientsTimeDependentExtendedValidity(TestSparseSHCoefficientsTimeDependentDefault):
    validity = (2010.0, 2018.0)
    options = {"validity_start": 2010.0, "validity_end": 2018.0}

#-------------------------------------------------------------------------------

if __name__ == "__main__":
    main()
