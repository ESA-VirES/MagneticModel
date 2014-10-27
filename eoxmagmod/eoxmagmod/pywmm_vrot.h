/*-----------------------------------------------------------------------------
 *
 * World Magnetic Model - C python bindings - vector rotation
 *  (i.e., vector coordinate system transformation)
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

#ifndef PYWMM_VROT_H
#define PYWMM_VROT_H

#include "math_aux.h"
#include "geo_conv.h"
#include "pywmm_aux.h"

#include <stdio.h>

/* recursive vector rotation */

static void _vrot_sph2geod(ARRAY_DATA arrd_in, ARRAY_DATA arrd_dlat, ARRAY_DATA arrd_out)
{
    if (arrd_in.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < arrd_in.dim[0]; ++i)
            _vrot_sph2geod(_get_arrd_item(&arrd_in, i),
                _get_arrd_item(&arrd_dlat, i), _get_arrd_item(&arrd_out, i));
    }
    else
    {
        #define SV(a) (*((double*)((a).data)))
        #define P(a,i) ((double*)((a).data+(i)*(a).stride[0]))
        #define V(a,i) (*P(a,i))

        const double a = DG2RAD*SV(arrd_dlat);
        const double sin_a = sin(a);
        const double cos_a = cos(a);

        V(arrd_out, 1) = V(arrd_in, 1);
        rot2d( P(arrd_out, 2), P(arrd_out, 0), V(arrd_in, 2), V(arrd_in, 0), sin_a, cos_a);

        #undef V
        #undef P
    }
}

/* python function definition */

#define DOC_VROT_SPH2GEOD "\n"\
"   arr_out = vrot_sph2geod(arr_in, arr_lat_sph, arr_lat_geod)\n"\
"\n"\
"     Rotate vectors from the spherical geocentric to the geodetic\n"\
"     coordinate system (or back) for a given difference of latitudes.\n"\
"     The inputs are:\n"\
"         arr_in - array of the input vectors\n"\
"         arr_dlat - difference of the latitudes in dg..\n"

static PyObject* vrot_sph2geod(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"arr_in", "arr_dlat"};
    PyObject *obj_in = NULL; // input object
    PyObject *obj_dlat = NULL; // input object
    PyObject *arr_in = NULL; // input array
    PyObject *arr_dlat = NULL; // input array
    PyObject *arr_out = NULL; // output array
    PyObject *retval = NULL;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwdict,
            "OO|:vrot_sph2geod", keywords, &obj_in, &obj_dlat));

    // cast the objects to arrays
    if (NULL == (arr_in=_get_as_double_array(obj_in, 1, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_dlat=_get_as_double_array(obj_dlat, 0, 0, NPY_ALIGNED, keywords[1])))
        goto exit;

    // check maximum allowed input array dimension
    if (PyArray_NDIM(arr_in) > MAX_OUT_ARRAY_NDIM)
    {
        PyErr_Format(PyExc_ValueError, "The input array dimension of '%s'"\
            " %d exceeds the allowed maximum value %d!", keywords[0],
            PyArray_NDIM(arr_in), MAX_OUT_ARRAY_NDIM);
        goto exit;
    }

    // check the dimensions
    if (_check_array_dim_eq(arr_in, -1, 3, keywords[0]))
        goto exit;

    if ((PyArray_NDIM(arr_in) > 1)||(PyArray_NDIM(arr_dlat) > 0))
    {
        int d;
        if (_check_array_dim_eq(arr_dlat, -1, 1, keywords[1]))
            goto exit;

        if (PyArray_NDIM(arr_in) != PyArray_NDIM(arr_dlat))
        {
            PyErr_Format(PyExc_ValueError, "Shape mismatch between '%s' and "
                "'%s'!", keywords[0], keywords[1]);
            goto exit;
        }

        for (d = 0 ; d < (PyArray_NDIM(arr_in)-1); ++d)
        {
            if (PyArray_DIM(arr_in, d) != PyArray_DIM(arr_dlat, d))
            {
                PyErr_Format(PyExc_ValueError, "Shape mismatch between '%s' "
                    "and '%s'!", keywords[0], keywords[1]);
                goto exit;
            }
        }
    }

    // create the output array
    if (NULL == (arr_out = _get_new_double_array(PyArray_NDIM(arr_in), PyArray_DIMS(arr_in), 3)))
        goto exit;

    // rorate the vector(s)
    _vrot_sph2geod(_array_to_arrd(arr_in), _array_to_arrd(arr_dlat), _array_to_arrd(arr_out));

    retval = arr_out;

  exit:

    // decrease reference counters to the arrays
    if (arr_in){Py_DECREF(arr_in);}
    if (arr_dlat){Py_DECREF(arr_dlat);}
    if (!retval && arr_out){Py_DECREF(arr_out);}

    return retval;
}

#endif  /* PYWMM_VROT_H */
