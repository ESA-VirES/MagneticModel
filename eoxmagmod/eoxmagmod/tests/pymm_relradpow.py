#-------------------------------------------------------------------------------
#
#  Spherical Harmonic Expansion - Geomagnetic Model - tests
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
# pylint: disable=missing-docstring

from unittest import TestCase, main
from numpy import asarray, empty, linspace
from numpy.testing import assert_allclose
from eoxmagmod._pymm import relradpow


class TestRadialPowerSeries(TestCase):

    @classmethod
    def reference_internal(cls, rads, degree):
        return cls._reference(rads, degree, cls.reference_internal_scalar)

    @classmethod
    def reference_external(cls, rads, degree):
        return cls._reference(rads, degree, cls.reference_external_scalar)

    @staticmethod
    def _reference(rads, degree, method):
        """ Evaluate relative-radius power series. """
        rads = asarray(rads)
        size = rads.size
        shape = rads.shape

        nterms = degree + 1

        rads = rads.ravel()
        rrp = empty((size, nterms))

        for i in range(size):
            rrp[i, :] = method(rads[i], degree)

        return rrp.reshape((*shape, nterms))

    @staticmethod
    def reference_internal_scalar(value, degree):
        """ Reference implementation - internal field. """
        return asarray([value ** (-i - 2) for i in range(degree + 1)])

    @staticmethod
    def reference_external_scalar(value, degree):
        """ Reference implementation - external field. """
        return asarray([value ** (i - 1) for i in range(degree + 1)])

    def test_invalid_degree(self):
        self.assertRaises(ValueError, relradpow, 1.0, -1, 1.0)

    def test_invalid_reference_radius(self):
        self.assertRaises(ValueError, relradpow, 1.0, 0, -1.0)
        self.assertRaises(ValueError, relradpow, 1.0, 0, 0.0)

    def test_relradpow_zero_degree_external(self):
        assert_allclose(relradpow(2.0, 0, 1.0, is_internal=False), [0.5])

    def test_relradpow_zero_degree_internal(self):
        assert_allclose(relradpow(2.0, 0, 1.0, is_internal=True), [0.25])

    def test_relradpow_default(self):
        assert_allclose(
            relradpow(2*6371.2, 10),
            [
                0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625,
                0.001953125, 0.0009765625, 0.00048828125, 0.000244140625,
            ]
        )
        assert_allclose(
            relradpow(1.1*6371.2, 256),
            self.reference_internal(1.1, 256)
        )

    def test_relradpow_internal(self):
        assert_allclose(
            relradpow(2*6371.2, 10, is_internal=True),
            [
                0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625,
                0.001953125, 0.0009765625, 0.00048828125, 0.000244140625,
            ]
        )
        assert_allclose(
            relradpow(1.1*6371.2, 256, is_internal=True),
            self.reference_internal(1.1, 256)
        )

    def test_relradpow_external_with_custom_radius(self):
        assert_allclose(
            relradpow(4.0, 10, reference_radius=2.0, is_internal=False),
            [0.5, 1, 2, 4, 8, 16, 32, 64, 128, 256, 512]
        )
        assert_allclose(
            relradpow(1.1, 256, reference_radius=1.0, is_internal=False),
            self.reference_external(1.1, 256)
        )

    def test_relradpow_internal_with_custom_radius(self):
        assert_allclose(
            relradpow(4.0, 10, reference_radius=2.0, is_internal=True),
            [
                0.25, 0.125, 0.0625, 0.03125, 0.015625, 0.0078125, 0.00390625,
                0.001953125, 0.0009765625, 0.00048828125, 0.000244140625,
            ]
        )
        assert_allclose(
            relradpow(1.1, 256, reference_radius=1.0, is_internal=True),
            self.reference_internal(1.1, 256)
        )

    def test_relradpow_internal_array_input(self):
        rads = linspace(1.0, 2.4, 15)
        assert_allclose(
            relradpow(rads, 256, reference_radius=1.0, is_internal=True),
            self.reference_internal(rads, 256)
        )

    def test_relradpow_external_array_input(self):
        rads = linspace(1.0, 2.4, 15)
        assert_allclose(
            relradpow(rads, 256, reference_radius=1.0, is_internal=False),
            self.reference_external(rads, 256)
        )


if __name__ == "__main__":
    main()
