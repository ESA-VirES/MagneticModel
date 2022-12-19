/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - bisection interval search
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

#ifndef PYMM_BISECT_H
#define PYMM_BISECT_H

#include "bisect.h"
#include "pymm_aux.h"

#define BISECT_SIDE_LEFT 0
#define BISECT_SIDE_RIGHT 1

static void _bisect_left(ARRAY_DATA *arrd_in, ARRAY_DATA *arrd_out, ARRAY_DATA *arrd_nodes);
static void _bisect_right(ARRAY_DATA *arrd_in, ARRAY_DATA *arrd_out, ARRAY_DATA *arrd_nodes);

/* Python function definition */

#define DOC_BISECT "\n"\
"   idx = bisect(v, x, side=BISECT_SIDE_LEFT)\n"\
"\n"\
"     Find indices of intervals, defined by the given sorted array or nodes,\n"\
"     the values fit in.\n"\
"\n"\
"     The input parameters are:\n"\
"       v - sorted 1D array defining the searched intervals.\n"\
"       x - array of values for which the interval are searched\n"\
"       side - define the open side of the searched interval\n"\
"              BISECT_SIDE_LEFT  v[idx] < x <= v[idx+1]\n"\
"              BISECT_SIDE_RIGHT v[idx] <= x < v[idx+1]\n"\
"\n"


static PyObject* bisect(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"v", "x", "side", NULL};

    int side = BISECT_SIDE_LEFT;
    PyArrayObject *arr_x = NULL; // x array (input)
    PyArrayObject *arr_v = NULL; // v array (input)
    PyArrayObject *arr_i = NULL; // i array (output)
    PyObject *obj_in_x = NULL; // input object
    PyObject *obj_in_v = NULL; // input object
    PyObject *retval = NULL; // output value

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OO|i:bisect", keywords, &obj_in_v, &obj_in_x, &side
    ))
        goto exit;

    if ((side != BISECT_SIDE_LEFT)&&(side != BISECT_SIDE_RIGHT))
    {
        PyErr_Format(PyExc_ValueError, "Invalid side value!", keywords[2]);
        goto exit;
    }

    // cast the input objects to arrays
    if (NULL == (arr_v = _get_as_double_array(obj_in_v, 1, 1, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_x = _get_as_double_array(obj_in_x, 0, 0, NPY_ARRAY_ALIGNED, keywords[1])))
        goto exit;

    // create the output arrays
    if (NULL == (arr_i = _get_new_array(PyArray_NDIM(arr_x), PyArray_DIMS(arr_x), NPY_INTP))) {
        goto exit;
    }

    // perform the bisection interval searches
    {
        ARRAY_DATA arrd_x = _array_to_arrd(arr_x);
        ARRAY_DATA arrd_i = _array_to_arrd(arr_i);
        ARRAY_DATA arrd_v = _array_to_arrd(arr_v);

        void (*_bisect)(ARRAY_DATA *, ARRAY_DATA *, ARRAY_DATA *) = (
            side == BISECT_SIDE_RIGHT ? _bisect_right : _bisect_left
        );

        _bisect(&arrd_x, &arrd_i, &arrd_v);
    }

    // get return value
    retval = (PyObject*) arr_i;

  exit:

    // decrease reference counters to the arrays
    if (arr_x) Py_DECREF(arr_x);
    if (arr_v) Py_DECREF(arr_v);
    if (!retval && arr_i) Py_DECREF(arr_i);

    return retval;
}


/*
 * high level nD-array interval search
 */

static void _bisect_left(ARRAY_DATA *arrd_in, ARRAY_DATA *arrd_out, ARRAY_DATA *arrd_nodes)
{
    if (arrd_in->ndim > 0)
    {
        npy_intp i, n = arrd_in->dim[0];

        for(i = 0; i < n; ++i) {
            ARRAY_DATA arrd_in_item = _get_arrd_item_nocheck(arrd_in, i);
            ARRAY_DATA arrd_out_item = _get_arrd_item_nocheck(arrd_out, i);
            _bisect_left(
                &arrd_in_item,
                &arrd_out_item,
                arrd_nodes
            );
        }
    }
    else
    {
        const size_t n = arrd_nodes->dim[0];
        const double *v = (const double*)arrd_nodes->data;
        const double x = *((double*)arrd_in->data);
        npy_intp *idx = (npy_intp*)arrd_out->data;

        *idx = bisect_left(x, v, n);
    }
}


static void _bisect_right(ARRAY_DATA *arrd_in, ARRAY_DATA *arrd_out, ARRAY_DATA *arrd_nodes)
{
    if (arrd_in->ndim > 0)
    {
        npy_intp i, n = arrd_in->dim[0];

        for(i = 0; i < n; ++i) {
            ARRAY_DATA arrd_in_item = _get_arrd_item_nocheck(arrd_in, i);
            ARRAY_DATA arrd_out_item = _get_arrd_item_nocheck(arrd_out, i);
            _bisect_right(
                &arrd_in_item,
                &arrd_out_item,
                arrd_nodes
            );
        }
    }
    else
    {
        const size_t n = arrd_nodes->dim[0];
        const double *v = (const double*)arrd_nodes->data;
        const double x = *((double*)arrd_in->data);
        npy_intp *idx = (npy_intp*)arrd_out->data;

        *idx = bisect_right(x, v, n);
    }
}

#endif /*PYMM_BISECT_H*/
