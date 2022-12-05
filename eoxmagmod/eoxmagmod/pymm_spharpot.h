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

#ifndef PYMM_SPHARPOT_H
#define PYMM_SPHARPOT_H

#include "spherical_harmonics.h"
#include "pymm_aux.h"
#include "pymm_cconv.h"
#include "pymm_sphar_common.h"

void _spharpot(
    ARRAY_DATA *arrd_pot, ARRAY_DATA *arrd_rad, ARRAY_DATA *arrd_cg,
    ARRAY_DATA *arrd_ch, ARRAY_DATA *arrd_lp, ARRAY_DATA *arrd_rrp,
    ARRAY_DATA *arrd_lsin, ARRAY_DATA *arrd_lcos, const int degree
);

/* Python function definition */

#define DOC_SPHARPOT "\n"\
"  v = spharpot(radius, coef_g, coef_h, leg_p, rrp, lonsin, loncos, degree=-1)\n"\
"\n"\
"     Spherical harmonic evaluation of the scalar potential field\n"\
"     in the (geocentric) spherical coordinates (latitude, longitude, radius).\n"\
"\n"\
"     The input parameters are:\n"\
"       radius - radius, i.e., distance from the earth centre, at the\n"\
"                evaluated location.\n"\
"       coef_g - vector of spherical harmonic model coefficients.\n"\
"       coef_h - vector of spherical harmonic model coefficients.\n"\
"       leg_p - vector the Legendre polynomials.\n"\
"       rrp - vector the relative radius powers.\n"\
"       lonsin - vector of the longitude cosines.\n"\
"       lonsin - vector of the longitude sines.\n"\
"       degree - degree of the spherical harmonic model. If not provided or \n"\
"                set to a negative number, then it is derived from the array \n"\
"                sizes.\n"\


static PyObject* spharpot(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "radius", "coef_g", "coef_h", "leg_p", "rpp", "lonsin", "loncos",
        "degree", NULL
    };

    int degree = -1;

    PyObject *obj_rad = NULL; // radius object
    PyObject *obj_cg = NULL; // coef_g object
    PyObject *obj_ch = NULL; // coef_h object
    PyObject *obj_lp = NULL; // P object
    PyObject *obj_rrp = NULL; // RRP object
    PyObject *obj_lsin = NULL; // lon-sin object
    PyObject *obj_lcos = NULL; // lon-cos object
    PyObject *retval = NULL; // output

    PyArrayObject *arr_pot = NULL; // output array
    PyArrayObject *arr_rad = NULL; // radius array
    PyArrayObject *arr_cg = NULL; // coef_g array
    PyArrayObject *arr_ch = NULL; // coef_h array
    PyArrayObject *arr_lp = NULL; // P array
    PyArrayObject *arr_rrp = NULL; // RRP array
    PyArrayObject *arr_lsin = NULL; // lon-sin array
    PyArrayObject *arr_lcos = NULL; // lon-cos array

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOOOOOO|i:spharpot", keywords,
        &obj_rad, &obj_cg, &obj_ch, &obj_lp, &obj_rrp, &obj_lsin, &obj_lcos,
        &degree
    ))
        goto exit;

    // cast the objects to arrays
    if (NULL == (arr_rad = _get_as_double_array(obj_rad, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_cg = _get_as_double_array(obj_cg, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[1])))
        goto exit;

    if (NULL == (arr_ch = _get_as_double_array(obj_ch, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[2])))
        goto exit;

    if (NULL == (arr_lp = _get_as_double_array(obj_lp, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[3])))
        goto exit;

    if (NULL == (arr_rrp = _get_as_double_array(obj_rrp, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[5])))
        goto exit;

    if (NULL == (arr_lsin = _get_as_double_array(obj_lsin, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[6])))
        goto exit;

    if (NULL == (arr_lcos = _get_as_double_array(obj_lcos, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[7])))
        goto exit;

    // extract degree from the array dimensions
    { 
        const int ndegrees = 6;
        npy_intp degrees[] = {
            _size_to_degree(PyArray_DIMS(arr_cg)[PyArray_NDIM(arr_cg)-1]),
            _size_to_degree(PyArray_DIMS(arr_ch)[PyArray_NDIM(arr_ch)-1]),
            _size_to_degree(PyArray_DIMS(arr_lp)[PyArray_NDIM(arr_lp)-1]),
            PyArray_DIMS(arr_rrp)[PyArray_NDIM(arr_rrp)-1] - 1,
            PyArray_DIMS(arr_lsin)[PyArray_NDIM(arr_lsin)-1] - 1,
            PyArray_DIMS(arr_lcos)[PyArray_NDIM(arr_lcos)-1] - 1
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
        const int narr = 7;
        PyArrayObject *arr[] = {
            arr_rad, arr_cg, arr_ch, arr_lp, arr_rrp, arr_lsin, arr_lcos,
        };
        npy_intp arr_ndim[] = {0, 1, 1, 1, 1, 1, 1};

        int arg_idx;
        npy_intp ndim, *dims;

        // extract maximum common shape
        _extract_common_shape(&arg_idx, &ndim, &dims, arr, arr_ndim, narr);

        // check if arrays match the common shape
        {
            int i;

            for (i = 0; i < 7; ++i)
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
                    }
                }
            }
        }

        // allocate the output array using the common shape
        if (NULL == (arr_pot = _get_new_array(ndim, dims, NPY_DOUBLE)))
            goto exit;
    }

    retval = (PyObject*) arr_pot;

    // evaluate potential
    {
        ARRAY_DATA arrd_pot = _array_to_arrd(arr_pot);
        ARRAY_DATA arrd_rad = _array_to_arrd(arr_rad);
        ARRAY_DATA arrd_cg = _array_to_arrd(arr_cg);
        ARRAY_DATA arrd_ch = _array_to_arrd(arr_ch);
        ARRAY_DATA arrd_lp = _array_to_arrd(arr_lp);
        ARRAY_DATA arrd_rrp = _array_to_arrd(arr_rrp);
        ARRAY_DATA arrd_lsin = _array_to_arrd(arr_lsin);
        ARRAY_DATA arrd_lcos = _array_to_arrd(arr_lcos);

        _spharpot(
            &arrd_pot,
            &arrd_rad,
            &arrd_cg,
            &arrd_ch,
            &arrd_lp,
            &arrd_rrp,
            &arrd_lsin,
            &arrd_lcos,
            degree
        );
    }

  exit:

    // decrease reference counters to the arrays
    if (arr_rad) {Py_DECREF(arr_rad);}
    if (arr_cg) {Py_DECREF(arr_cg);}
    if (arr_ch) {Py_DECREF(arr_ch);}
    if (arr_lp) {Py_DECREF(arr_lp);}
    if (arr_rrp) {Py_DECREF(arr_rrp);}
    if (arr_lsin) {Py_DECREF(arr_lsin);}
    if (arr_lcos) {Py_DECREF(arr_lcos);}
    if (!retval && arr_pot) {Py_DECREF(arr_pot);}

    return retval;
}

/*
 * Recursively iterate over the input arrays and evaluate the SH series.
 */

void _spharpot(
    ARRAY_DATA *arrd_pot, ARRAY_DATA *arrd_rad, ARRAY_DATA *arrd_cg,
    ARRAY_DATA *arrd_ch, ARRAY_DATA *arrd_lp, ARRAY_DATA *arrd_rrp,
    ARRAY_DATA *arrd_lsin, ARRAY_DATA *arrd_lcos, const int degree
)
{
    if (arrd_pot->ndim > 0)
    {
        npy_intp i, n = arrd_pot->dim[0];

        for(i = 0; i < n; ++i)
        {
            ARRAY_DATA arrd_pot_item = _get_arrd_item_nocheck(arrd_pot, i);
            ARRAY_DATA arrd_rad_item = _get_arrd_item(arrd_rad, i);
            ARRAY_DATA arrd_cg_item = _get_arrd_item_with_guard(arrd_cg, i, 1);
            ARRAY_DATA arrd_ch_item = _get_arrd_item_with_guard(arrd_ch, i, 1);
            ARRAY_DATA arrd_lp_item = _get_arrd_item_with_guard(arrd_lp, i, 1);
            ARRAY_DATA arrd_rrp_item = _get_arrd_item_with_guard(arrd_rrp, i, 1);
            ARRAY_DATA arrd_lsin_item = _get_arrd_item_with_guard(arrd_lsin, i, 1);
            ARRAY_DATA arrd_lcos_item = _get_arrd_item_with_guard(arrd_lcos, i, 1);

            _spharpot(
                &arrd_pot_item,
                &arrd_rad_item,
                &arrd_cg_item,
                &arrd_ch_item,
                &arrd_lp_item,
                &arrd_rrp_item,
                &arrd_lsin_item,
                &arrd_lcos_item,
                degree
            );
        }
    }
    else
    {
        double *pot = ((double*)arrd_pot->data);
        const double rad = *((double*)arrd_rad->data);
        const double *cg = ((double*)arrd_cg->data);
        const double *ch = ((double*)arrd_ch->data);
        const double *lp = ((double*)arrd_lp->data);
        const double *rrp = ((double*)arrd_rrp->data);
        const double *lsin = ((double*)arrd_lsin->data);
        const double *lcos = ((double*)arrd_lcos->data);

        shc_eval_v(pot, degree, rad, cg, ch, lp, rrp, lsin, lcos);
    }
}
#endif  /* PYMM_SPHARPOT_H */
