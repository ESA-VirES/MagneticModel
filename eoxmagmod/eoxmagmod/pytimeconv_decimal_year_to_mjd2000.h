/*-----------------------------------------------------------------------------
 *
 * Decimal year to MJD2000 conversion
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2018-2022 EOX IT Services GmbH
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

#ifndef PYTIMECONV_DECIMAL_YEAR_TO_MJD2000_H
#define PYTIMECONV_DECIMAL_YEAR_TO_MJD2000_H

#include "pymm_aux.h"
#include "time_conversion.h"
#include <math.h>

/*
 *  nD-array recursive evaluation
 */

static void _year_to_mjd2k(ARRAY_DATA arrd_in, ARRAY_DATA arrd_out)
{
    if (arrd_in.ndim > 0)
    {
        npy_intp i;
        for(i = 0; i < arrd_in.dim[0]; ++i)
            _year_to_mjd2k(
                _get_arrd_item(&arrd_in, i),
                _get_arrd_item(&arrd_out, i)
            );
    }
    else
    {
        #define S(a) (*((double*)(a).data))
        S(arrd_out) = decimal_year_to_mjd2k(S(arrd_in));
        #undef S 
    }
}


/*
 *  Python function definition
 */

#define DOC_DECIMAL_YEAR_TO_MJD2000 "\n"\
"   time_mjd2k = decimal_year_to_mjd2000(decimal_year)\n\n"\
"     Output:\n"\
"       time_mjd2k - array of the calculated MJD2000 times.\n"\
"\n"\
"     Parameters:\n"\
"       decimal_year - array of decimal years (up to 15 dimensions).\n"\
"\n"

static PyObject* pytimeconv_decimal_year_to_mjd2000(
    PyObject *self, PyObject *args, PyObject *kwdict
)
{
    static char *keywords[] = {"time_mjd2k", NULL};
    PyObject *obj_in = NULL; // input object
    PyArrayObject *arr_in = NULL; // input array
    PyArrayObject *arr_out = NULL; // return array
    PyObject *retval = NULL; // return object

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "O|:decimal_year_to_mjd2000", keywords, &obj_in
    ))
        goto exit;

    // cast input objects to arrays
    if (NULL == (arr_in = _get_as_double_array(obj_in, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    // check maximum allowed input array dimension
    if (PyArray_NDIM(arr_in) > (MAX_OUT_ARRAY_NDIM-1))
    {
        PyErr_Format(PyExc_ValueError, "Array dimension of '%s'"\
            " %d exceeds the allowed maximum value %d!", keywords[0],
            PyArray_NDIM(arr_in), (MAX_OUT_ARRAY_NDIM-1));
        goto exit;
    }

    // create a new output array
    if (NULL == (arr_out = (PyArrayObject*) PyArray_EMPTY(PyArray_NDIM(arr_in), PyArray_DIMS(arr_in), NPY_DOUBLE, 0)))
        goto exit;

    // convert decimal years to MJD2000 times
    _year_to_mjd2k(
        _array_to_arrd(arr_in),
        _array_to_arrd(arr_out)
    );

    // assign return value
    retval = (PyObject*) arr_out;

  exit:

    // decrease reference counters to the arrays
    if (arr_in) Py_DECREF(arr_in);
    if (!retval && arr_out) Py_DECREF(arr_out);

    return retval;
}

#endif /* PYTIMECONV_DECIMAL_YEAR_TO_MJD2000_H */
