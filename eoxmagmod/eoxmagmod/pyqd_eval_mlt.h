/*-----------------------------------------------------------------------------
 *
 * Magnetic Quasi Dipole Coordinates - C python bindings
 * -  Magnetic Local Time evaluation
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2016-2024 EOX IT Services GmbH
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

#ifndef PYQD_EVAL_MLT_H
#define PYQD_EVAL_MLT_H

#include "pymm_aux.h"
#include "qdipole/cqdipole.h"

/* Python function definition */

#define DOC_EVAL_MLT "\n"\
"   mlt = eval_mlt(qdlon, time, fname)\n"\
"     Inputs:\n"\
"       qdlon - quasi-dipole longitudes(s).\n"\
"       time  - MJD2000 time(s)\n"\
"     Outputs:\n"\
"       mlt - magnetic local times (s).\n"\
""

static PyObject* eval_mlt(PyObject *self, PyObject *args, PyObject *kwdict)
{
    int status;
    static char *keywords[] = {"qdlon", "time", NULL};

    PyObject *obj_qdlon = NULL; // gclon object
    PyObject *obj_time = NULL; // time object
    PyArrayObject *arr_qdlon = NULL; // qdlon array
    PyArrayObject *arr_time = NULL; // time array
    PyArrayObject *arr_mlt = NULL; // mlt array
    PyObject *retval = NULL;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OO:eval_mlt", keywords, &obj_qdlon, &obj_time
    ))
        goto exit;

    #define NPY_REQ (NPY_ARRAY_ALIGNED|NPY_ARRAY_C_CONTIGUOUS)

    // cast the objects to arrays
    if (NULL == (arr_qdlon = _get_as_double_array(obj_qdlon, 0, 1, NPY_REQ, keywords[1])))
        goto exit;

    if (NULL == (arr_time = _get_as_double_array(obj_time, 0, 1, NPY_REQ, keywords[3])))
        goto exit;

    // check the dimensions
    npy_intp ndim = PyArray_NDIM(arr_qdlon);
    npy_intp *dims = PyArray_DIMS(arr_qdlon);

    if(_check_arr_dims_all_eq(arr_time, ndim, dims, keywords[3]))
        goto exit;

    // create the output arrays
    if (NULL == (arr_mlt = (PyArrayObject*) PyArray_EMPTY(ndim, dims, NPY_DOUBLE, 0)))
        goto exit;

    // evaluate the output values
    status = c_eval_mlt(
        (double*) PyArray_DATA(arr_mlt),
        (double*) PyArray_DATA(arr_qdlon),
        (double*) PyArray_DATA(arr_time),
        ndim == 0 ? 1 : dims[0]
    );

    if (status) {
        PyErr_Format(
            PyExc_RuntimeError,
            "Call to c_eval_mlt() failed with an error! error_code = %d", status
        );
        goto exit;
    }

    retval = (PyObject*) arr_mlt;

  exit:

    // decrease reference counters to the arrays
    if (arr_qdlon) Py_DECREF(arr_qdlon);
    if (arr_time) Py_DECREF(arr_time);
    if (!retval && arr_mlt) Py_DECREF(arr_mlt);

    return retval;
}

#endif  /* PYQD_EVAL_MLT_H */
