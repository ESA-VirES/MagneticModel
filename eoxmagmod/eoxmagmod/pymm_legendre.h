/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - associative Legendre functions evaluation
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2014-2022 EOX IT Services GmbH
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

#ifndef PYMM_LEGENDRE_H
#define PYMM_LEGENDRE_H

#include "spherical_harmonics.h"
#include "pymm_aux.h"
#include "pymm_cconv.h"

/* Python function definition */

static void _legendre(
    ARRAY_DATA *arrd_lat, ARRAY_DATA *arrd_p, ARRAY_DATA *arrd_dp,
    const int degree, const double *psqrt, const double *prsqrt
);


#define DOC_LEGENDRE "\n"\
"   p, dp = legendre(latitude, degree)\n"\
"\n"\
"     For given the latitude in degrees and the degree evaluate\n"\
"     the Schmidt semi-normalised associated Legendre functions and\n"\
"     their derivatives: \n"\
"         P_n^m(sin(latitude))  and  dP_n^m(sin(latitude))\n"\
"     where n = 0..degree and m = 0..n \n"\
"\n"\
"     The input parameters are:\n"\
"       latitude - spherical latitude in degrees (-90 <= latitude <= 90). It can be an array.\n"\
"                  It can be an array of values.\n"\
"       degree - degree of the spherical harmonic model.\n"\
"\n"\
"     Two arrays are returned, one for the values of the associated \n"\
"     Legendre functions and one for their derivatives.\n"\
"\n"


static PyObject* legendre(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"latitude", "degree", NULL};

    int degree;
    double *psqrt = NULL;
    double *prsqrt = NULL;
    PyArrayObject *arr_lat = NULL; // input array of latitudes
    PyArrayObject *arr_p = NULL; // P array
    PyArrayObject *arr_dp = NULL; // dP array
    PyObject *obj_lat = NULL; // input object
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "Oi:legendre", keywords, &obj_lat, &degree
    ))
        goto exit;

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "%s < 0", keywords[1]);
        goto exit;
    }

    // cast the input object to an array
    if (NULL == (arr_lat = _get_as_double_array(obj_lat, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    // create the output array
    {
        npy_intp ndim = PyArray_NDIM(arr_lat) + 1;
        npy_intp dims[ndim];
        npy_intp i;

        for (i = 0; i < ndim - 1; ++i)
        {
            dims[i] = PyArray_DIMS(arr_lat)[i];
        }

        dims[ndim-1] = (((npy_intp)degree + 1)*((npy_intp)degree + 2))/2;

        if (NULL == (arr_p = _get_new_array(ndim, dims, NPY_DOUBLE)))
            goto exit;

        if (NULL == (arr_dp = _get_new_array(ndim, dims, NPY_DOUBLE)))
            goto exit;
    }

    // allocate and fill the pre-calculated square-roots
    if (NULL == (psqrt = shc_presqrt(degree)))
    {
        PyErr_Format(PyExc_MemoryError, "Memory allocation error!");
        goto exit;
    }

    // allocate and fill the pre-calculated reciprocal square-roots
    if (NULL == (prsqrt = shc_prersqrt(degree)))
    {
        PyErr_Format(PyExc_MemoryError, "Memory allocation error!");
        goto exit;
    }


    if (NULL == (retval = Py_BuildValue("NN", arr_p, arr_dp)))
        goto exit;


    // evaluate associated Legendre functions and their derivatives
    {
        ARRAY_DATA arrd_lat = _array_to_arrd(arr_lat);
        ARRAY_DATA arrd_p = _array_to_arrd(arr_p);
        ARRAY_DATA arrd_dp = _array_to_arrd(arr_dp);

        _legendre(
            &arrd_lat,
            &arrd_p,
            &arrd_dp,
            degree,
            psqrt,
            prsqrt
        );
    }

  exit:

    // free the auxiliary arrays
    if (prsqrt) free(prsqrt);
    if (psqrt) free(psqrt);

    // decrease reference counters to the arrays
    if (arr_lat) Py_DECREF(arr_lat);
    if (!retval && arr_p) Py_DECREF(arr_p);
    if (!retval && arr_dp) Py_DECREF(arr_dp);

    return retval;
}

/*
 * Recursively iterate over the latitude array values and evaluate the Legendre
 * functions for each of them.
 */

static void _legendre(
    ARRAY_DATA *arrd_lat, ARRAY_DATA *arrd_p, ARRAY_DATA *arrd_dp,
    const int degree, const double *psqrt, const double *prsqrt
)
{
    if (arrd_lat->ndim > 0)
    {
        npy_intp i, n = arrd_lat->dim[0];

        for(i = 0; i < n; ++i)
        {
            ARRAY_DATA arrd_lat_item = _get_arrd_item_nocheck(arrd_lat, i);
            ARRAY_DATA arrd_p_item = _get_arrd_item_nocheck(arrd_p, i);
            ARRAY_DATA arrd_dp_item = _get_arrd_item_nocheck(arrd_dp, i);

            _legendre(
                &arrd_lat_item,
                &arrd_p_item,
                &arrd_dp_item,
                degree,
                psqrt,
                prsqrt
            );
        }
    }
    else
    {
        const double lat_sph = *((double*)arrd_lat->data) * DG2RAD;
        double *p = ((double*)arrd_p->data);
        double *dp = ((double*)arrd_dp->data);

        shc_legendre(p, dp, degree, lat_sph, psqrt, prsqrt);
    }
}

#endif  /* PYMM_LEGENDRE_H */
