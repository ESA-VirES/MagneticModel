#-------------------------------------------------------------------------------
#
#  Magnetic Model
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
# all
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

import os.path
import datetime
import time
import numpy as np

# coordinate systems and their transformation
from _pywmm import (
    GEODETIC_ABOVE_WGS84, GEODETIC_ABOVE_EGM96,
    GEOCENTRIC_SPHERICAL, GEOCENTRIC_CARTESIAN,
    POTENTIAL, GRADIENT, POTENTIAL_AND_GRADIENT,
    convert, legendre, lonsincos, relradpow, spharpot, sphargrd,
    vrot_sph2geod, vrot_sph2cart, vrot_cart2sph, sheval,
)

COORD_TYPES = (
    (GEODETIC_ABOVE_WGS84, "GEODETIC_ABOVE_WGS84"),
    (GEODETIC_ABOVE_EGM96, "GEODETIC_ABOVE_EGM96"),
    (GEOCENTRIC_SPHERICAL, "GEOCENTRIC_SPHERICAL"),
    (GEOCENTRIC_CARTESIAN, "GEOCENTRIC_CARTESIAN"),
)

_GEODETIC_COORD_TYPES = (GEODETIC_ABOVE_WGS84, GEODETIC_ABOVE_EGM96)

EVAL_MODES = (
    (POTENTIAL, "POTENTIAL"),
    (GRADIENT, "GRADIENT"),
    (POTENTIAL_AND_GRADIENT, "POTENTIAL_AND_GRADIENT"),
)

def vrotate(arr, coord_in, coord_out, coord_type_in, coord_type_out):
    """ Rotate vectors from one coordinate system to another.
        Input:
            arr - array of the source vectors
            coord_in - source coordinates
            coord_out - destination coordinates
            coord_type_in - source coordinate system type
            coord_type_out - destination coordinate system type
        Output:
            arr_out - array of the rotated vectors
    """
    if coord_type_in == coord_type_out:
        return arr

    if coord_type_in in _GEODETIC_COORD_TYPES:
        if coord_type_out in _GEODETIC_COORD_TYPES:
            return arr
        elif coord_type_out == GEOCENTRIC_SPHERICAL:
            return vrot_sph2geod(arr, coord_out[..., 0] - coord_in[..., 0])
        elif coord_type_out == GEOCENTRIC_CARTESIAN:
            return vrot_sph2cart(arr, coord_in[..., 0], coord_in[..., 1])

    elif coord_type_in == GEOCENTRIC_SPHERICAL:
        if coord_type_out in _GEODETIC_COORD_TYPES:
            return vrot_sph2geod(arr, coord_out[..., 0] - coord_in[..., 0])
        elif coord_type_out == GEOCENTRIC_CARTESIAN:
            return vrot_sph2cart(arr, coord_in[..., 0], coord_in[..., 1])

    elif coord_type_in == GEOCENTRIC_CARTESIAN:
        if coord_type_out in _GEODETIC_COORD_TYPES:
            return vrot_cart2sph(arr, coord_out[..., 0], coord_out[..., 1])
        elif coord_type_out == GEOCENTRIC_SPHERICAL:
            return vrot_cart2sph(arr, coord_out[..., 0], coord_out[..., 1])

    raise ValueError("Unsupported coordinate type!")


def vnorm(arr):
    """Calculate norms for each vector form an input array of vectors."""
    return np.sqrt((arr*arr).sum(axis=arr.ndim-1))

def vincdecnorm(arr):
    """ Calculate vector inclinations (-90:90), declinations
    (-180,180), and the vector norms.
    """
    # equivalent to conversion of Cartesian to spherical coordinates
    tmp = convert(arr, GEOCENTRIC_CARTESIAN, GEOCENTRIC_SPHERICAL)
    return -tmp[...,0], tmp[...,1], tmp[...,2]


def to_year_fraction(date):
    """ Converts a Python `datetime.date` or `datetime.datetime` to a decimal
    year value.
    """

    def since_epoch(date):  # returns seconds since epoch
        return time.mktime(date.timetuple())

    year = date.year
    start_this_year = datetime.datetime(year=year, month=1, day=1)
    start_next_year = datetime.datetime(year=year+1, month=1, day=1)

    year_elapsed = since_epoch(date) - since_epoch(start_this_year)
    year_duration = since_epoch(start_next_year) - since_epoch(start_this_year)
    fraction = year_elapsed / year_duration

    return date.year + fraction


class MagneticModelBase(object):
    """ Base Magnetic model class """

    @property
    def validity(self):
        """Get interval of model validity."""
        raise NotImplementedError

    def is_valid_date(self, date):
        """Check whether the date is within the interval of validity."""
        vrange = self.validity
        return (date >= vrange[0]) and (date <= vrange[1])

    @property
    def degree_static(self):
        """Get the degree of the static model."""
        raise NotImplementedError

    @property
    def degree_secvar(self):
        """Get the degree of the secular variation model."""
        raise NotImplementedError

    def get_coef_static(self, date):
        """ Calculate model static coefficients for a date specified
        by a decimal year value.
        """
        raise NotImplementedError

    def get_coef_secvar(self, date):
        """Get secular variation coefficients."""
        raise NotImplementedError

    def print_info(self):
        """Print information about the model."""
        print self

    def eval(self, arr_in, date, coord_type_in=GEODETIC_ABOVE_WGS84,
                coord_type_out=None, secvar=False, mode=GRADIENT, maxdegree=-1,
                mindegree=-1, check_validity=True):
        """Evaluate spherical harmonic model for a given set of spatio-temporal
        coordinates.

        Input:
            arr_in - numpy array of (x, y, z) coordinates.
                     The type of the x, y, z, values depends on the selected
                     input coordinate system.

            date - The time as a decimal year value (e.g., 2015.23).
            coord_type_in - shall be set to one of the valid coordinate systems:
                        GEODETIC_ABOVE_WGS84 (default)
                        GEODETIC_ABOVE_EGM96
                        GEOCENTRIC_SPHERICAL
                        GEOCENTRIC_CARTESIAN

            coord_type_out - coordinate system of the output vector
                        GEODETIC_ABOVE_WGS84
                        GEODETIC_ABOVE_EGM96 (the same as WGS84)
                        GEOCENTRIC_SPHERICAL
                        GEOCENTRIC_CARTESIAN
                    Output coordinate system defaults to the input coordinate
                    system.

            secvar - if True secular variation of the magnetic field is
                     calculated rather than the magnetic field.

            mode - type of output to be produced. The possible values are:
                        POTENTIAL
                        GRADIENT (default)
                        POTENTIAL_AND_GRADIENT

            maxdegree - an optional maximum allowed modelled degree
                    (i.e., truncated evaluation). If set to -1 no limit
                    is imposed.

            mindegree - an optional min. allowed modelled degree
                    (i.e., truncated evaluation). When applied any coefficient
                    below this degree are set to zero. If set to -1 no limit
                    is imposed.

            check_validity - boolean flag controlling  whether the date
                    validity will is checked (True, by default) or not (False).

        Output:

            arr_out - output numpy array with the same shape as the input
                      array containing the calculated magnetic field parameters.
        """

        if isinstance(date, (datetime.date, datetime.datetime)):
            date = to_year_fraction(date)

        if mode not in dict(EVAL_MODES):
            raise ValueError("Invalid mode value!")

        if coord_type_in not in dict(COORD_TYPES):
            raise ValueError("Invalid input coordinate type!")

        if coord_type_out is None:
            coord_type_out = coord_type_in

        if coord_type_out not in dict(COORD_TYPES):
            raise ValueError("Invalid output coordinate type!")

        if check_validity and not self.is_valid_date(date):
            raise ValueError("The date is outside of the model validity range!")

        if secvar:
            degree = self.degree_secvar
            coef_g, coef_h = self.get_coef_secvar(date)
        else:
            degree = self.degree_static
            coef_g, coef_h = self.get_coef_static(date)

        if mindegree >= 0:
            mindegree = min(degree + 1, mindegree)
            idx = ((mindegree+1)*mindegree)/2
            coef_g[:idx] = 0
            coef_h[:idx] = 0

        if maxdegree >= 0:
            degree = min(maxdegree, degree)

        if mindegree > degree:
            # skip pointless evaluation when all coefficients are zero
            degree = 0

        return sheval(
            arr_in, degree, coef_g, coef_h, coord_type_in, coord_type_out, mode
        )


    def field_line(self, point, date, coord_type_in=GEODETIC_ABOVE_WGS84,
                coord_type_out=None, maxdegree=-1, check_validity=True,
                step=1e2, nstep=500):
        """Trace a field line passing trough given point in space.
        (Equivalent to ODE solution by means of the Euler method.)

        Input:
            point - as a numpy array of (x, y, z) coordinates.
                     The type of the x, y, z, values depends on the selected
                     input coordinate system.

            date - The time as a decimal year value (e.g., 2015.23).
            coord_type_in - shall be set to one of the valid coordinate systems:
                        GEODETIC_ABOVE_WGS84 (default)
                        GEODETIC_ABOVE_EGM96
                        GEOCENTRIC_SPHERICAL
                        GEOCENTRIC_CARTESIAN

            coord_type_out - coordinate system of the output vector
                        GEODETIC_ABOVE_WGS84
                        GEODETIC_ABOVE_EGM96 (the same as WGS84)
                        GEOCENTRIC_SPHERICAL
                        GEOCENTRIC_CARTESIAN
                    Output coordinate system defaults to the input coordinate
                    system.

            maxdegree - an optional maximum allowed modelled degree
                    (i.e., truncated evaluation). If set to -1 no limit
                    is imposed.

            check_validity - boolean flag controlling  whether the date
                    validity is checked (True, by default) or not (False).

        Output:
            arr_coord - array of point coordinates in the requested coordinate
                    system.
            arr_field - array of magnetic field vectors in the requested
                    coordinate system corresponding to the point coordinates.
        """
        if coord_type_out is None:
            coord_type_out = coord_type_in

        prm = {
            'maxdegree': maxdegree,
            'check_validity': check_validity,
        }

        def _field_line_(coord_cart0, date, step, nstep, **mprm):

            def feval(coord_cart):
                field_vector = self.eval(
                    coord_cart, date, GEOCENTRIC_CARTESIAN,
                    GEOCENTRIC_CARTESIAN, secvar=False, **mprm
                )
                relative_radius = vnorm(coord_cart) / 6571.0
                field_vector_norm = vnorm(field_vector)
                if field_vector_norm > 0:
                    normalised_field_vector = field_vector / field_vector_norm
                else:
                    normalised_field_vector = field_vector
                return relative_radius * normalised_field_vector, field_vector

            coord_cart = coord_cart0
            coords_cart = [coord_cart]
            field_vectors = []
            for _ in xrange(nstep):
                direction, field_vector = feval(coord_cart)
                field_vectors.append(field_vector)
                coord_cart = coord_cart + step * direction
                coord_gdt = convert(
                    coord_cart, GEOCENTRIC_CARTESIAN, GEODETIC_ABOVE_WGS84
                )
                coords_cart.append(coord_cart)
                if coord_gdt[2] < 0.0:
                    break
            _, field_vector = feval(coord_cart)
            field_vectors.append(field_vector)
            return coords_cart, field_vectors

        # initial coordinate conversion
        coord_cart = convert(point, coord_type_in, GEOCENTRIC_CARTESIAN)
        # trace the field-lines both to the north and south
        coords_n, field_n = _field_line_(coord_cart, date, +step, nstep, **prm)
        coords_s, field_s = _field_line_(coord_cart, date, -step, nstep, **prm)
        # join the north and reversed south segments
        coords_s.reverse()
        field_s.reverse()
        coords_cart = np.array(coords_s[:-1] + coords_n)
        field_cart = np.array(field_s[:-1] + field_n)
        # final coordinate system conversions
        coords_out = convert(coords_cart, GEOCENTRIC_CARTESIAN, coord_type_out)
        field_out = vrotate(
            field_cart, None, coords_out, GEOCENTRIC_CARTESIAN, coord_type_out
        )
        return (coords_out, field_out)


class MagneticModelComposed(MagneticModelBase):
    """ Composed Magnetic model."""

    def __init__(self, model0, model1, c0=1.0, c1=1.0):
        """ Model constructor """
        super(MagneticModelComposed, self).__init__()
        self.model0 = model0
        self.model1 = model1
        self.c0 = c0
        self.c1 = c1

    @property
    def validity(self):
        """Get interval of model validity."""
        v0 = self.model0.validity
        v1 = self.model1.validity
        return (max(v0[0], v1[0]), min(v0[1], v1[1]))

    @property
    def degree_static(self):
        """Get the degree of the static model."""
        return max(self.model0.degree_static, self.model1.degree_static)

    @property
    def degree_secvar(self):
        """Get the degree of the secular variation model."""
        return max(self.model0.degree_secvar, self.model1.degree_secvar)

    @staticmethod
    def _combine_coef(f0, (cg0, ch0), f1, (cg1, ch1)):
        cg = np.zeros(max(cg0.size, cg1.size))
        ch = np.zeros(max(ch0.size, ch1.size))
        cg[:cg0.size] = f0 * cg0
        ch[:ch0.size] = f0 * ch0
        cg[:cg1.size] += f1 * cg1
        ch[:ch1.size] += f1 * ch1
        return (cg, ch)

    def get_coef_static(self, date):
        """ Calculate model static coefficients for a date specified
        by a decimal year value.
        """
        return self._combine_coef(
            self.c0, self.model0.get_coef_static(date),
            self.c1, self.model1.get_coef_static(date),
        )

    def get_coef_secvar(self, date):
        """Get secular variation coefficients."""
        return self._combine_coef(
            self.c0, self.model0.get_coef_secvar(date),
            self.c1, self.model1.get_coef_secvar(date),
        )

    def print_info(self):
        """Print information about the model."""
        print "Composed Model (%g,%g) {"%(self.c0, self.c1)
        self.model0.print_info()
        self.model1.print_info()
        print "}"

    def __add__(self, other):
        return MagneticModelComposed(self, other, 1.0, 1.0)

    def __sub__(self, other):
        return MagneticModelComposed(self, other, 1.0, -1.0)


class MagneticModel(MagneticModelBase):
    """ Magnetic model class """

    def __init__(self, model_prm):
        """ Model constructor """
        super(MagneticModel, self).__init__()
        self.prm = model_prm

    def __add__(self, other):
        return MagneticModelComposed(self, other, 1.0, 1.0)

    def __sub__(self, other):
        return MagneticModelComposed(self, other, 1.0, -1.0)
