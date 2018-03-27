/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - spherical-harmonics
 * - evaluation of the potential gradient in the spherical coordinate system
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2014 EOX IT Services GmbH
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies of this Software or works derived from this Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *-----------------------------------------------------------------------------
*/

#ifndef PYWMM_SPHARGRD_H
#define PYWMM_SPHARGRD_H

#include "shc.h"
#include "pywmm_aux.h"
#include "pywmm_cconv.h"

/* python function definition */
#define DOC_SPHARGRD "\n"\
"  v_grad = sphargrd(latitude, degree, coef_g, coef_h, leg_p, leg_dp, rrp, lonsin, loncos, spherical=True)\n"\
"\n"\
"     Spherical harmonic evaluation of the gradient of the potential\n"\
"     (scalar) field in the (geocentric) spherical coordinates (latitude,\n"\
"     longitude, radius).\n"\
"     The input parameters are:\n"\
"       latitude - spherical (or geodetic) latitude in dg. at the evaluated\n"\
"                  location.\n"\
"       degree - degree of the spherical harmonic model.\n"\
"       coef_g - vector of spherical harmonic model coeficients.\n"\
"       coef_h - vector of spherical harmonic model coeficients.\n"\
"       leg_p - vector the Legendre polynomials.\n"\
"       leg_dp - vector the Legendre polynomials' derivations.\n"\
"       rrp - vector the relative radius powers.\n"\
"       lonsin - vector the the longitude cosines.\n"\
"       lonsin - vector the the longitude sines.\n"\
"       spherical - boolean flag indicating whether a geodentic spherical\n"\
"                   (default, True) or geodetic (WGS84, False) latitude is\n"\
"                   being used.\n"

static PyObject* sphargrd(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"latitude", "degree", "coef_g", "coef_h",
        "leg_p", "leg_dp", "rpp", "lonsin", "loncos", "spherical", NULL};

    int degree, nterm;
    double lat_sph, lat_in, tmp0, tmp1;
    int is_sph = 0;

    PyObject *obj_cg = NULL; // coef_g object
    PyObject *obj_ch = NULL; // coef_h object
    PyObject *obj_lp = NULL; // P object
    PyObject *obj_ldp = NULL; // dP object
    PyObject *obj_rrp = NULL; // rel.rad.pow. object
    PyObject *obj_lsin = NULL; // lonsin object
    PyObject *obj_lcos = NULL; // loncos object

    PyObject *arr_out = NULL; // output array
    PyObject *arr_cg = NULL; // coef_g array
    PyObject *arr_ch = NULL; // coef_h array
    PyObject *arr_lp = NULL; // P array
    PyObject *arr_ldp = NULL; // dP array
    PyObject *arr_rrp = NULL; // rel.rad.pow. array
    PyObject *arr_lsin = NULL; // lonsin array
    PyObject *arr_lcos = NULL; // loncos array

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "diOOOOOOO|i:sphargrd", keywords,
        &lat_in, &degree, &obj_cg, &obj_ch, &obj_lp, &obj_ldp, &obj_rrp,
        &obj_lsin, &obj_lcos, &is_sph
    ))
        goto exit;

    if (degree < 1)
    {
        PyErr_Format(PyExc_ValueError, "Invalid value %d of '%s'!", degree, keywords[1]);
        goto exit;
    }

    nterm = ((degree+1)*(degree+2))/2;

    // convert latitude to geocentric spherical latitude
    if (is_sph)
        lat_sph = lat_in;
    else
        conv_WGS84_to_sphECEF(&lat_sph, &tmp0, &tmp1, lat_in, 0.0, 0.0);

    // cast the objects to arrays
    if (NULL == (arr_cg=_get_as_double_array(obj_cg, 1, 1, NPY_IN_ARRAY, keywords[2])))
        goto exit;

    if (NULL == (arr_ch=_get_as_double_array(obj_ch, 1, 1, NPY_IN_ARRAY, keywords[3])))
        goto exit;

    if (NULL == (arr_lp=_get_as_double_array(obj_lp, 1, 1, NPY_IN_ARRAY, keywords[4])))
        goto exit;

    if (NULL == (arr_ldp=_get_as_double_array(obj_ldp, 1, 1, NPY_IN_ARRAY, keywords[5])))
        goto exit;

    if (NULL == (arr_rrp=_get_as_double_array(obj_rrp, 1, 1, NPY_IN_ARRAY, keywords[6])))
        goto exit;

    if (NULL == (arr_lsin=_get_as_double_array(obj_lsin, 1, 1, NPY_IN_ARRAY, keywords[7])))
        goto exit;

    if (NULL == (arr_lcos=_get_as_double_array(obj_lcos, 1, 1, NPY_IN_ARRAY, keywords[8])))
        goto exit;

    // check the arrays' dimensions
    if (_check_array_dim_le(arr_cg, 0, nterm, keywords[2]))
        goto exit;

    if (_check_array_dim_le(arr_ch, 0, nterm, keywords[3]))
        goto exit;

    if (_check_array_dim_le(arr_lp, 0, nterm, keywords[4]))
        goto exit;

    if (_check_array_dim_le(arr_ldp, 0, nterm, keywords[5]))
        goto exit;

    if (_check_array_dim_le(arr_rrp, 0, degree+1, keywords[6]))
        goto exit;

    if (_check_array_dim_le(arr_lsin, 0, degree+1, keywords[7]))
        goto exit;

    if (_check_array_dim_le(arr_lcos, 0, degree+1, keywords[8]))
        goto exit;

    // allocate the output array
    if (NULL == (arr_out = _get_new_double_array(1, NULL, 3)))
        goto exit;

    // the evaluation
    {
        double *out = (double*)PyArray_DATA(arr_out);
        shc_eval(NULL, out+0, out+1, out+2, degree, 0x2,
                DG2RAD*lat_sph, 0.0, PyArray_DATA(arr_cg),
                PyArray_DATA(arr_ch), PyArray_DATA(arr_lp),
                PyArray_DATA(arr_ldp), PyArray_DATA(arr_rrp),
                PyArray_DATA(arr_lsin), PyArray_DATA(arr_lcos));
    }

  exit:

    // decrease reference counters to the arrays
    if (arr_cg){Py_DECREF(arr_cg);}
    if (arr_ch){Py_DECREF(arr_ch);}
    if (arr_lp){Py_DECREF(arr_lp);}
    if (arr_ldp){Py_DECREF(arr_ldp);}
    if (arr_rrp){Py_DECREF(arr_rrp);}
    if (arr_lsin){Py_DECREF(arr_lsin);}
    if (arr_lcos){Py_DECREF(arr_lcos);}

    return arr_out;
 }

#endif  /* PYWMM_SPHARGRD_H */
