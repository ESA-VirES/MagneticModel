/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - evaluation of the relative radius power series
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

#ifndef PYMM_RELRADPOW_H
#define PYMM_RELRADPOW_H

#include "spherical_harmonics.h"
#include "pymm_aux.h"

/* Earth radius in km */
#define RADIUS  6371.2

#ifndef STR
#define SRINGIFY(x) #x
#define STR(x) SRINGIFY(x)
#endif

static void _relradpow(
    ARRAY_DATA *arrd_rad, ARRAY_DATA *arrd_rrp,
    const int degree, const int is_internal, const double reference_radius
);

/* Python function definition */

#define DOC_RELRADPOW "\n"\
"   rrp = relradpow(radius, degree, reference_radius="STR(RADIUS)", is_internal=True)\n"\
"\n"\
"    Calculate relative radius series for the given radius or array\n"\
"    of radii.\n"\
"\n"\
"    By default when the 'is_internal' flag is set to True, evaluate\n"\
"    for the given 'radius' and 'reference_radius' relative radius power\n"\
"    series:\n"\
"      (reference_radius/radius)**(i+2) for i in range(0, degree+1) .\n"\
"\n"\
"    When the 'is_internal' flag is set to False, evaluate\n"\
"    for the given 'radius' and 'reference_radius' relative radius power\n"\
"    series:\n"\
"      (radius/reference_radius)**(i+1) for i in range(0, degree+1) .\n"\
"\n"

static PyObject* relradpow(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "radius", "degree", "reference_radius", "is_internal", NULL
    };

    int degree;
    int is_internal;
    double reference_radius = RADIUS; // reference radius
    PyArrayObject *arr_rad = NULL; // input array of radii
    PyArrayObject *arr_rrp = NULL; // radial series array
    PyObject *obj_is_internal = NULL; // boolean flag
    PyObject *obj_rad = NULL; // input object
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "Oi|dO:relradpow", keywords,
        &obj_rad, &degree, &reference_radius, &obj_is_internal
    ))
        goto exit;

    is_internal = (obj_is_internal == NULL) || PyObject_IsTrue(obj_is_internal);

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "%s < 0", keywords[1]);
        goto exit;
    }

    if (reference_radius <= 0.0)
    {
        PyErr_Format(PyExc_ValueError, "%s <= 0", keywords[2]);
        goto exit;
    }

    // cast the input object to an array
    if (NULL == (arr_rad = _get_as_double_array(obj_rad, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    // create the output array
    {
        npy_intp ndim = PyArray_NDIM(arr_rad) + 1;
        npy_intp dims[ndim];
        npy_intp i;

        for (i = 0; i < ndim - 1; ++i)
        {
            dims[i] = PyArray_DIMS(arr_rad)[i];
        }

        dims[ndim-1] = (npy_intp)degree + 1;

        if (NULL == (arr_rrp = _get_new_array(ndim, dims, NPY_DOUBLE)))
            goto exit;
    }

    retval = (PyObject*) arr_rrp;

    // evaluate relative radial power series
    {
        ARRAY_DATA arrd_rad = _array_to_arrd(arr_rad);
        ARRAY_DATA arrd_rrp = _array_to_arrd(arr_rrp);

        _relradpow(
            &arrd_rad,
            &arrd_rrp,
            degree,
            is_internal,
            reference_radius
        );
    }

  exit:

    // decrease reference counters to the arrays
    if (arr_rad) Py_DECREF(arr_rad);
    if (!retval && arr_rrp) Py_DECREF(arr_rrp);

    return retval;
}


/*
 * Recursively iterate over the radius array values and evaluate the radial
 * functions for each of them.
 */

static void _relradpow(
    ARRAY_DATA *arrd_rad, ARRAY_DATA *arrd_rrp,
    const int degree, const int is_internal, const double reference_radius
)
{
    if (arrd_rad->ndim > 0)
    {
        npy_intp i, n = arrd_rad->dim[0];

        for(i = 0; i < n; ++i)
        {
            ARRAY_DATA arrd_rad_item = _get_arrd_item_nocheck(arrd_rad, i);
            ARRAY_DATA arrd_rrp_item = _get_arrd_item_nocheck(arrd_rrp, i);

            _relradpow(
                &arrd_rad_item,
                &arrd_rrp_item,
                degree,
                is_internal,
                reference_radius
            );
        }
    }
    else
    {
        const double rel_rad = *((double*)arrd_rad->data) / reference_radius;
        double *rrp = ((double*)arrd_rrp->data);

        (is_internal ? shc_relradpow_internal : shc_relradpow_external)(rrp, degree, rel_rad);
    }
}

#endif  /* PYMM_RELRADPOW_H */
