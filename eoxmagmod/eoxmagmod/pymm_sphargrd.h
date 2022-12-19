/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - spherical-harmonics
 * - evaluation of the potential gradient in the spherical coordinate system
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

#ifndef PYMM_SPHARGRD_H
#define PYMM_SPHARGRD_H

#include "spherical_harmonics.h"
#include "pymm_aux.h"
#include "pymm_cconv.h"
#include "pymm_sphar_common.h"

void _sphargrd(
    ARRAY_DATA *arrd_grd, ARRAY_DATA *arrd_lat, ARRAY_DATA *arrd_coef,
    ARRAY_DATA *arrd_lp, ARRAY_DATA *arrd_ldp, ARRAY_DATA *arrd_rrp,
    ARRAY_DATA *arrd_lcs, const int degree, const int is_internal
);

/* Python function definition */

#define DOC_SPHARGRD "\n"\
"  v_grad = sphargrd(latitude, coef, leg_p, leg_dp, rrp, lcs, is_internal=True, degree=-1)\n"\
"\n"\
"     Spherical harmonic evaluation of the gradient of the potential\n"\
"     (scalar) field in the geocentric spherical coordinates (latitude,\n"\
"     longitude, radius).\n"\
"\n"\
"     The input parameters are:\n"\
"       latitude - spherical (or geodetic) latitude in degrees at the evaluated\n"\
"                  location.\n"\
"       coef - array of spherical harmonic model coefficients.\n"\
"       leg_p - array of Legendre polynomials.\n"\
"       leg_dp - array of Legendre polynomials' derivations.\n"\
"       rrp - array of relative radius powers.\n"\
"       lcs - array of longitude cosines and sines.\n"\
"       is_internal - boolean flag set to True by default. When set to False\n"\
"                     external field evaluation is used.\n"\
"       degree - degree of the spherical harmonic model. If not provided or \n"\
"                set to a negative number, then it is derived from the array \n"\
"                sizes.\n"\
"\n"\

static PyObject* sphargrd(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "latitude", "coef", "leg_p", "leg_dp",
        "rrp", "lcs", "is_internal", "degree", NULL
    };

    int degree = -1, is_internal;

    PyObject *obj_is_internal = NULL; // boolean flag
    PyObject *obj_lat = NULL; // latitude object
    PyObject *obj_coef = NULL; // coef_g object
    PyObject *obj_lp = NULL; // P object
    PyObject *obj_ldp = NULL; // dP object
    PyObject *obj_rrp = NULL; // RRP object
    PyObject *obj_lcs = NULL; // lon cos/sin object
    PyObject *retval = NULL; // output

    PyArrayObject *arr_grd = NULL; // output array
    PyArrayObject *arr_lat = NULL; // latitude array
    PyArrayObject *arr_coef = NULL; // coef_g array
    PyArrayObject *arr_lp = NULL; // P array
    PyArrayObject *arr_ldp = NULL; // dP array
    PyArrayObject *arr_rrp = NULL; // RRP array
    PyArrayObject *arr_lcs = NULL; // lon cos/sin array

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOOOOO|Oi:sphargrd", keywords,
        &obj_lat, &obj_coef, &obj_lp, &obj_ldp, &obj_rrp, &obj_lcs,
        &obj_is_internal, &degree
    ))
        goto exit;

    is_internal = (obj_is_internal == NULL) || PyObject_IsTrue(obj_is_internal);

    // cast the objects to arrays
    if (NULL == (arr_lat = _get_as_double_array(obj_lat, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_coef = _get_as_double_array(obj_coef, 2, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[1])))
        goto exit;

    if (NULL == (arr_lp = _get_as_double_array(obj_lp, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[3])))
        goto exit;

    if (NULL == (arr_ldp = _get_as_double_array(obj_ldp, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[4])))
        goto exit;

    if (NULL == (arr_rrp = _get_as_double_array(obj_rrp, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[5])))
        goto exit;

    if (NULL == (arr_lcs = _get_as_double_array(obj_lcs, 2, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[6])))
        goto exit;

    // extract degree from the array dimensions
    { 
        const int ndegrees = 5;
        npy_intp degrees[] = {
            _size_to_degree(PyArray_DIMS(arr_coef)[PyArray_NDIM(arr_coef)-2]),
            _size_to_degree(PyArray_DIMS(arr_lp)[PyArray_NDIM(arr_lp)-1]),
            _size_to_degree(PyArray_DIMS(arr_ldp)[PyArray_NDIM(arr_ldp)-1]),
            PyArray_DIMS(arr_rrp)[PyArray_NDIM(arr_rrp)-1] - 1,
            PyArray_DIMS(arr_lcs)[PyArray_NDIM(arr_lcs)-2] - 1,
        };

        npy_intp max_degree;
        int arg_idx;

         _get_max_degree(&arg_idx, &max_degree, degrees, ndegrees);

        if (max_degree < 0)
        {
            PyErr_Format(
                PyExc_ValueError,
                "Negative degree due to empty %s array!",
                keywords[arg_idx+1]
            );
            goto exit;
        }

        if ((max_degree > MAX_DEGREE)&&(degree < 0))
        {
            PyErr_Format(PyExc_ValueError, "Degree is larger than %g!", MAX_DEGREE);
            goto exit;
        }

        // set the applied degree
        if ((degree < 0)||(degree > max_degree))
            degree = max_degree;
    }

    // check array dimensions and allocate the output array
    {
        const int narr = 6;
        PyArrayObject *arr[] = {
            arr_lat, arr_coef, arr_lp, arr_ldp, arr_rrp, arr_lcs,
        };
        npy_intp arr_ndim[] = {0, 2, 1, 1, 1, 2};

        int arg_idx;
        npy_intp ndim, *dims;

        // extract maximum common shape
        _extract_common_shape(&arg_idx, &ndim, &dims, arr, arr_ndim, narr);

        // check if arrays match the common shape
        {
            int i;

            for (i = 0; i < narr; ++i)
            {
                if (i != arg_idx)
                {
                    if (_compare_dimensions(dims, PyArray_DIMS(arr[i]), PyArray_NDIM(arr[i]) - arr_ndim[i]))
                    {
                        PyErr_Format(
                            PyExc_ValueError,
                            "Shape of %s array does not match shape of %s!",
                            keywords[i], keywords[arg_idx]
                        );
                        goto exit;
                    }
                }

                // check last dimension of the pair-arrays
                if ((arr_ndim[i] == 2)&&(PyArray_DIMS(arr[i])[PyArray_NDIM(arr[i])-1] != 2))
                {
                    PyErr_Format(PyExc_ValueError, "Invalid shape of %s array!", keywords[i]);
                    goto exit;
                }
            }
        }

        // allocate the output array using the common shape
        {
            npy_intp j;
            npy_intp ndim_new = ndim + 1;
            npy_intp dims_new[ndim_new];

            for (j = 0; j < ndim; ++j)
                dims_new[j] = dims[j];
            dims_new[ndim] = 3;

            if (NULL == (arr_grd = _get_new_array(ndim_new, dims_new, NPY_DOUBLE)))
                goto exit;
        }
    }

    retval = (PyObject*) arr_grd;

    // evaluate gradient of the potential
    {
        ARRAY_DATA arrd_grd = _array_to_arrd(arr_grd);
        ARRAY_DATA arrd_lat = _array_to_arrd(arr_lat);
        ARRAY_DATA arrd_coef = _array_to_arrd(arr_coef);
        ARRAY_DATA arrd_lp = _array_to_arrd(arr_lp);
        ARRAY_DATA arrd_ldp = _array_to_arrd(arr_ldp);
        ARRAY_DATA arrd_rrp = _array_to_arrd(arr_rrp);
        ARRAY_DATA arrd_lcs = _array_to_arrd(arr_lcs);

        _sphargrd(
            &arrd_grd,
            &arrd_lat,
            &arrd_coef,
            &arrd_lp,
            &arrd_ldp,
            &arrd_rrp,
            &arrd_lcs,
            degree,
            is_internal
        );
    }

  exit:

    // decrease reference counters to the arrays
    if (arr_lat) {Py_DECREF(arr_lat);}
    if (arr_coef) {Py_DECREF(arr_coef);}
    if (arr_lp) {Py_DECREF(arr_lp);}
    if (arr_ldp) {Py_DECREF(arr_ldp);}
    if (arr_rrp) {Py_DECREF(arr_rrp);}
    if (arr_lcs) {Py_DECREF(arr_lcs);}
    if (!retval && arr_grd) {Py_DECREF(arr_grd);}

    return retval;
 }


/*
 * Recursively iterate over the input arrays and evaluate the SH series.
 */

void _sphargrd(
    ARRAY_DATA *arrd_grd, ARRAY_DATA *arrd_lat, ARRAY_DATA *arrd_coef,
    ARRAY_DATA *arrd_lp, ARRAY_DATA *arrd_ldp, ARRAY_DATA *arrd_rrp,
    ARRAY_DATA *arrd_lcs, const int degree, const int is_internal
)
{
    if (arrd_grd->ndim > 1)
    {
        npy_intp i, n = arrd_grd->dim[0];

        for(i = 0; i < n; ++i)
        {
            ARRAY_DATA arrd_grd_item = _get_arrd_item_nocheck(arrd_grd, i);
            ARRAY_DATA arrd_lat_item = _get_arrd_item(arrd_lat, i);
            ARRAY_DATA arrd_coef_item = _get_arrd_item_with_guard(arrd_coef, i, 2);
            ARRAY_DATA arrd_lp_item = _get_arrd_item_with_guard(arrd_lp, i, 1);
            ARRAY_DATA arrd_ldp_item = _get_arrd_item_with_guard(arrd_ldp, i, 1);
            ARRAY_DATA arrd_rrp_item = _get_arrd_item_with_guard(arrd_rrp, i, 1);
            ARRAY_DATA arrd_lcs_item = _get_arrd_item_with_guard(arrd_lcs, i, 2);

            _sphargrd(
                &arrd_grd_item,
                &arrd_lat_item,
                &arrd_coef_item,
                &arrd_lp_item,
                &arrd_ldp_item,
                &arrd_rrp_item,
                &arrd_lcs_item,
                degree,
                is_internal
            );
        }
    }
    else
    {
        double *grd = ((double*)arrd_grd->data);
        const double lat = *((double*)arrd_lat->data) * DG2RAD;
        const double (*coef)[2] = ((double(*)[2])arrd_coef->data);
        const double *lp = ((double*)arrd_lp->data);
        const double *ldp = ((double*)arrd_ldp->data);
        const double *rrp = ((double*)arrd_rrp->data);
        const double (*lcs)[2] = ((double(*)[2])arrd_lcs->data);

        shc_eval_dv(
            grd+0, grd+1, grd+2, degree, lat,
            coef, lcs, lp, ldp, rrp, is_internal
        );
    }
}

#endif  /* PYMM_SPHARGRD_H */
