/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - coefficients time interpolation
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *-----------------------------------------------------------------------------
 * Copyright (C) 2022 EOX IT Services GmbH
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

#ifndef PYMM_INTERP_H
#define PYMM_INTERP_H

#include "math.h"
#include "interpolation.h"
#include "pymm_aux.h"

#define INTERP_C0 0x00
#define INTERP_C1 0x01
#define INTERP_C1D1 0x11

#ifndef NAN
#define NAN (0.0/0.0)
#endif

static void _interp(ARRAY_DATA *arrd_t, ARRAY_DATA *arrd_c,
                    ARRAY_DATA *arrd_t0, ARRAY_DATA *arrd_c0,
                    f_get_interp_basis get_interp_basis,
                    f_interp_eval interp_eval, int extrapolate);

/* Python function definition */

#define DOC_INTERP "\n"\
"   coeff = interp(time, time0, coeff0, kind=INTERP1, extrapolate=False)\n"\
"\n"\
"     Interpolate time series of SH coefficients.\n"\
"\n"\
"     The input parameters are:\n"\
"       time - time at which the coefficients are interpolated.\n"\
"       time0[N] - times of the interpolated series of coefficients.\n"\
"       coeff0[...,N] - series of the interpolated SH coefficients\n"\
"       kind - interpolation kind. Currently supported\n"\
"           INTERP_C0   ... piecewise constant\n"\
"           INTERP_C1   ... piecewise linear\n"\
"           INTERP_C1D1 ... piecewise linear - 1st derivative\n"\
"       extrapolate - if True the out of bound values will be extrapolated\n"\
"\n"\
"     Output:\n"\
"       coeff[...] - coefficients at the given time\n"\
"\n"

static PyObject* interp(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "time", "time0", "coeff0", "kind", "extrapolate", NULL
    };

    int extrapolate = 0;
    int kind = INTERP_C1;
    size_t min_nt = (size_t)-1;
    f_get_interp_basis get_interp_basis = NULL;
    f_interp_eval interp_eval = NULL;

    PyArrayObject *arr_t = NULL; // interpolation times array (input)
    PyArrayObject *arr_c = NULL; // calculated coefficients (ouput)
    PyArrayObject *arr_t0 = NULL; // interpolated times array (input)
    PyArrayObject *arr_c0 = NULL; // interpolated coefficients array (input)

    PyObject *obj_t = NULL; // input object
    PyObject *obj_t0 = NULL; // input object
    PyObject *obj_c0 = NULL; // input object
    PyObject *obj_extrapolate = NULL; // boolean flag
    PyObject *retval = NULL; // output value

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOO|iO:interp", keywords,
        &obj_t, &obj_t0, &obj_c0, &kind, &obj_extrapolate
    ))
        goto exit;

    if (obj_extrapolate != NULL) {
        extrapolate = PyObject_IsTrue(obj_extrapolate);
    }

    switch (kind) {
        case INTERP_C0:
            min_nt = 1;
            get_interp_basis = get_interp0_basis;
            interp_eval = interp0_eval;
            break;
        case INTERP_C1:
            min_nt = 2;
            get_interp_basis = get_interp1_basis;
            interp_eval = interp1_eval;
            break;
        case INTERP_C1D1:
            min_nt = 2;
            get_interp_basis = get_interp1d1_basis;
            interp_eval = interp1_eval;
            break;
        default:
            PyErr_Format(PyExc_ValueError, "Invalid kind value!", keywords[3]);
            goto exit;
    }

    // cast the input objects to arrays
    if (NULL == (arr_t = _get_as_double_array(obj_t, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_t0 = _get_as_double_array(obj_t0, 1, 1, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[1])))
        goto exit;

    if (NULL == (arr_c0 = _get_as_double_array(obj_c0, 1, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[2])))
        goto exit;

    // check the array dimensions
    if (PyArray_DIMS(arr_t0)[0] != PyArray_DIMS(arr_c0)[PyArray_NDIM(arr_c0)-1]) {
        PyErr_Format(PyExc_ValueError, "The coefficients series does no match the corresponding times!");
        goto exit;
    }
    // check the array dimensions
    if (PyArray_DIMS(arr_t0)[0] < min_nt) {
        PyErr_Format(PyExc_ValueError, "The time series is expected to have at least %ld node%s for the requested interpolation kind.", min_nt, min_nt == 1 ? "" : "s");
        goto exit;
    }

    // create the output array
    {
        npy_intp ndim = PyArray_NDIM(arr_t) + PyArray_NDIM(arr_c0) - 1;
        npy_intp dims[ndim];
        npy_intp i;

        for (i = 0; i < PyArray_NDIM(arr_t); ++i)
        {
            dims[i] = PyArray_DIMS(arr_t)[i];
        }

        for (i = 0; i < PyArray_NDIM(arr_c0) - 1; ++i)
        {
            dims[i + PyArray_NDIM(arr_t)] = PyArray_DIMS(arr_c0)[i];
        }

        if (NULL == (arr_c = _get_new_array(ndim, dims, NPY_DOUBLE))) {
            goto exit;
        }
    }

    // perform coefficients interpolation
    {
        ARRAY_DATA arrd_t = _array_to_arrd(arr_t);
        ARRAY_DATA arrd_c = _array_to_arrd(arr_c);
        ARRAY_DATA arrd_t0 = _array_to_arrd(arr_t0);
        ARRAY_DATA arrd_c0 = _array_to_arrd(arr_c0);

        _interp(
            &arrd_t,
            &arrd_c,
            &arrd_t0,
            &arrd_c0,
            get_interp_basis,
            interp_eval,
            extrapolate
        );
    }

    // get return value
    retval = (PyObject*) arr_c;

  exit:

    // decrease reference counters to the arrays
    if (arr_t) Py_DECREF(arr_t);
    if (arr_t0) Py_DECREF(arr_t0);
    if (arr_c0) Py_DECREF(arr_c0);
    if (!retval && arr_c) Py_DECREF(arr_c);

    return retval;
}


/*
 * high level nD-array recursive coefficient interpolation
 */

static void _fill_value(ARRAY_DATA *arrd, double value);
static void _interp_eval(ARRAY_DATA *arrd_c, ARRAY_DATA *arrd_c0, INTERP_BASIS *basis, f_interp_eval interp_eval);

static void _interp(ARRAY_DATA *arrd_t, ARRAY_DATA *arrd_c,
                    ARRAY_DATA *arrd_t0, ARRAY_DATA *arrd_c0,
                    f_get_interp_basis get_interp_basis,
                    f_interp_eval interp_eval, int extrapolate)
{
    if (arrd_t->ndim > 0)
    {
        npy_intp i, n = arrd_t->dim[0];

        for(i = 0; i < n; ++i)
        {

            ARRAY_DATA arrd_t_item = _get_arrd_item_nocheck(arrd_t, i);
            ARRAY_DATA arrd_c_item = _get_arrd_item_nocheck(arrd_c, i);
            _interp(
                &arrd_t_item,
                &arrd_c_item,
                arrd_t0,
                arrd_c0,
                get_interp_basis,
                interp_eval,
                extrapolate
            );
        }
    }
    else
    {
        const double t = *((double*)arrd_t->data);
        const double *t0 = ((double*)arrd_t0->data);
        const size_t nt = arrd_t0->dim[0];

        INTERP_BASIS basis = get_interp_basis(t, t0, nt);

        if (extrapolate||(basis.i0 == basis.i)||(t == t0[nt-1]))
        {
            _interp_eval(arrd_c, arrd_c0, &basis, interp_eval);
            return;
        }

        _fill_value(arrd_c, NAN);
    }
}


static void _interp_eval(ARRAY_DATA *arrd_c, ARRAY_DATA *arrd_c0, INTERP_BASIS *basis, f_interp_eval interp_eval)
{
    if (arrd_c->ndim > 0)
    {
        npy_intp i, n = arrd_c->dim[0];

        for(i = 0; i < n; ++i) {

            ARRAY_DATA arrd_c_item = _get_arrd_item_nocheck(arrd_c, i);
            ARRAY_DATA arrd_c0_item = _get_arrd_item_nocheck(arrd_c0, i);
            _interp_eval(
                &arrd_c_item,
                &arrd_c0_item,
                basis,
                interp_eval
            );
        }
    }
    else
    {
        double *c = ((double*)arrd_c->data);
        const double *c0 = ((double*)arrd_c0->data);

        *c = interp_eval(basis, c0);
    }
}


static void _fill_value(ARRAY_DATA *arrd, double value) {
    if (arrd->ndim > 0)
    {
        npy_intp i, n = arrd->dim[0];

        for(i = 0; i < n; ++i) {
            ARRAY_DATA arrd_item = _get_arrd_item_nocheck(arrd, i);
            _fill_value(&arrd_item, value);
        }
    }
    else
    {
        *((double*)arrd->data) = value;
    }
}

#endif /*PYMM_INTERP_H*/
