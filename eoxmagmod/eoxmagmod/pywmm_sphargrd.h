/*-----------------------------------------------------------------------------
 *
 * World Magnetic Model - C python bindings - spherical-harmonics
 * - evaluation of the potential gradient in the sherical coordinate system
 *
 * Project: World Magnetic Model - python interface
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

#include "pywmm_aux.h"
#include "pywmm_cconv.h"

/* python function definition */
#define DOC_SPHARGRD "\n"\
"  v_grad = sphargrd(latitude, degree, coef_g, coef_h, leg_p, leg_dp, rrp, lonsin, loncos, spherical=True)\n"\
"\n"\
"     Spherical harmonics evaluation of the gradient of the potential\n"\
"     (scalar) field in the (geocentric) spherical coordinates (latitude,\n"\
"     longitude, radius).\n"\
"     The input parameters are:\n"\
"       latitude - spherical (or geodetic) latitude in dg. of the evaluated\n"\
"                  location.\n"\
"       degree - degree of the spherical harmonic model.\n"\
"       coef_g - vector of spherical harmonic model coeficients.\n"\
"       coef_h - vector of spherical harmonic model coeficients.\n"\
"       leg_p - vector the Legendre polynomials.\n"\
"       leg_pp - vector the Legendre polynomials' derivations.\n"\
"       rrp - vector the relative radius powers.\n"\
"       lonsin - vector the the longitude cosines.\n"\
"       lonsin - vector the the longitude sines.\n"\
"       spherical - boolean flag indicating whether a geodentic spherical\n"\
"                   (default, True) or geodetic (WGS84, False) latitude is\n"\
"                   being used.\n"

static PyObject* sphargrd(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"latitude", "degree", "coef_g", "coef_h",
        "leg_p", "leg_p", "rpp", "lonsin", "loncos", "spherical", NULL};

    int degree, nterm;
    double lat_sph, lat_in, tmp0, tmp1;
    int is_sph = 0;

    PyObject *retval = NULL; // returned value
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
    if (!PyArg_ParseTupleAndKeywords(args, kwdict, "diOOOOOOO|i:", keywords,
            &lat_in, &degree, &obj_cg, &obj_ch, &obj_lp, &obj_ldp, &obj_rrp,
            &obj_lsin, &obj_lcos, &is_sph));

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
        const double *cg = (double*)PyArray_DATA(arr_cg);
        const double *ch = (double*)PyArray_DATA(arr_ch);
        const double *lp = (double*)PyArray_DATA(arr_lp);
        const double *ldp = (double*)PyArray_DATA(arr_ldp);
        const double *rrp = (double*)PyArray_DATA(arr_rrp);
        const double *lsin = (double*)PyArray_DATA(arr_lsin);
        const double *lcos = (double*)PyArray_DATA(arr_lcos);
        double *out = (double*)PyArray_DATA(arr_out);
        double cos_lat = cos(DG2RAD*lat_sph);
        double dv_lat = 0.0, dv_lon = 0.0, dv_rad = 0.0;
        int i, j;

        for (i = 1; i <= degree; ++i)
        {
            const int i_off = (i*(i+1))/2;

            for (j = 0; j <= i; ++j)
            {
                const int idx = i_off + j;
                const double tmp0 = cg[idx]*lcos[j] + ch[idx]*lsin[j];
                const double tmp1 = cg[idx]*lsin[j] - ch[idx]*lcos[j];

                dv_lat += tmp0 * rrp[i] * ldp[idx];
                dv_lon += tmp1 * rrp[i] * lp[idx] * j;
                dv_rad += tmp0 * rrp[i] * lp[idx] * (i+1);
            }
        }

        out[0] = -dv_lat;
        out[1] = +dv_lon/cos_lat;
        out[2] = -dv_rad;

        // handling of the geographic poles
        if (fabs(cos_lat) < 1e-10)
        {
            const double sin_lat = sin(DG2RAD*lat_sph);
            const double lsin1 = lsin[1], lcos1 = lcos[1];
            double sqn3, sqn1 = 1.0;
            double ps2, ps1 = 1.0, ps0 = 1.0,

            // i = 1
            dv_lon = (cg[2]*lsin1 - ch[2]*lcos1) * rrp[1];

            for (i = 2; i <= degree; ++i)
            {
                const int idx = 1 + (i*(i+1))/2;
                #define FDIV(a,b) ((double)(a)/(double)(b))
                const double tmp = FDIV((i-1)*(i-1)-1, (2*i-1)*(2*i-3));

                // evaluate ratio between the Gauss-normalised and Smidth
                // quasi-normalised associated Legendre functions.
                //  Equivalent to: sqrt((j==0?1:2)*(i-j)!/(i+j!))*(2i-1)!!/(i-j)!
                sqn1 = sqn1 * FDIV(2*i-1, i);
                sqn3 = sqn1 * sqrt(FDIV(i*2, i+1));
                #undef FDIV
                ps2 = ps1;
                ps1 = ps0;
                ps0 = sin_lat*ps1 - tmp*ps2;

                dv_lon += (cg[idx]*lsin1 - ch[idx]*lcos1) * rrp[i] * ps0 * sqn3;
            }

            out[1] = dv_lon;
        }
    }

    retval = arr_out;

  exit:

    // decrease reference counters to the arrays
    if (arr_cg){Py_DECREF(arr_cg);}
    if (arr_ch){Py_DECREF(arr_ch);}
    if (arr_lp){Py_DECREF(arr_lp);}
    if (arr_ldp){Py_DECREF(arr_ldp);}
    if (arr_rrp){Py_DECREF(arr_rrp);}
    if (arr_lsin){Py_DECREF(arr_lsin);}
    if (arr_lcos){Py_DECREF(arr_lcos);}
    if (!retval && arr_out){Py_DECREF(arr_out);}

    return retval;
 }

#endif  /* PYWMM_SPHARGRD_H */
