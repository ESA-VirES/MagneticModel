/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - longitude cos/sin spherical terms' evaluation
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

#ifndef PYMM_LONCOSSIN_H
#define PYMM_LONCOSSIN_H

#include <math.h>
#include "fourier_series.h"
#include "pymm_aux.h"

static void _loncossin(
    ARRAY_DATA *arrd_lon, ARRAY_DATA *arrd_cos_sin,
    const int degree, const int fast_alg
);

/* Python function definition */

#define DOC_LONCOSSIN "\n"\
"  loncossin = loncossin(longitude, degree, fast_alg=True)\n"\
"\n"\
"     For given longitude and degree, evaluate\n"\
"     the cosine and sine series: \n"\
"        [cos(i*longitude), sin(i*longitude)] for i in range(0, degree+1)\n"\
"     The longitude has to be entered in degrees. Array input is accepted.\n"\
"     The 'fast_alg' boolean options forces the subroutine to use a faster\n"\
"     but slightly less precise evaluation algorithm.\n"

static PyObject* loncossin(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"longitude", "degree", "fast_alg", NULL};

    int degree;
    int fast_alg = 1;
    PyArrayObject *arr_lon = NULL; // input array of longitudes
    PyArrayObject *arr_cos_sin = NULL; // cos/sin array
    PyObject *obj_lon = NULL; // input object
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "Oi|i:loncossin", keywords, &obj_lon, &degree, &fast_alg
    ))
        goto exit;

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "%s < 0", keywords[1]);
        goto exit;
    }


    // cast the input object to an array
    if (NULL == (arr_lon = _get_as_double_array(obj_lon, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    // create the output arrays
    {
        const npy_intp ndim = PyArray_NDIM(arr_lon) + 2;
        npy_intp dims[ndim];
        npy_intp i;

        for (i = 0; i < ndim - 2; ++i)
        {
            dims[i] = PyArray_DIMS(arr_lon)[i];
        }

        dims[ndim-2] = (npy_intp)degree + 1;
        dims[ndim-1] = 2;

        if (NULL == (arr_cos_sin = _get_new_array(ndim, dims, NPY_DOUBLE)))
            goto exit;
    }

    retval = (PyObject*) arr_cos_sin;

    // evaluate sin/cos series
    {
        ARRAY_DATA arrd_lon = _array_to_arrd(arr_lon);
        ARRAY_DATA arrd_cos_sin = _array_to_arrd(arr_cos_sin);

        _loncossin(
            &arrd_lon,
            &arrd_cos_sin,
            degree,
            fast_alg
        );
    }

  exit:

    // decrease reference counters to the arrays
    if (arr_lon) Py_DECREF(arr_lon);
    if (!retval && arr_cos_sin) Py_DECREF(arr_cos_sin);

    return retval;
}


/*
 * Recursively iterate over the longitude array values and evaluate the sin/cos
 * series.
 */

static void _loncossin(
    ARRAY_DATA *arrd_lon, ARRAY_DATA *arrd_cos_sin,
    const int degree, const int fast_alg
)
{
    if (arrd_lon->ndim > 0)
    {
        npy_intp i, n = arrd_lon->dim[0];

        for(i = 0; i < n; ++i)
        {
            ARRAY_DATA arrd_lon_item = _get_arrd_item_nocheck(arrd_lon, i);
            ARRAY_DATA arrd_cos_sin_item = _get_arrd_item_nocheck(arrd_cos_sin, i);

            _loncossin(
                &arrd_lon_item,
                &arrd_cos_sin_item,
                degree,
                fast_alg
            );
        }
    }
    else
    {
        const double lon = *((double*)arrd_lon->data) * DG2RAD;
        double (*cos_sin)[2] = ((double(*)[2])arrd_cos_sin->data);

        (fast_alg ? fs_cos_sin: fs_cos_sin_ref)(cos_sin, degree, lon);
    }
}

#endif  /* PYMM_LONCOSSIN_H */
