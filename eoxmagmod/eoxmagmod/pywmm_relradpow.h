/*-----------------------------------------------------------------------------
 *
 * World Magnetic Model - C python bindings
 * - evaluation of the relative radius power series
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

#ifndef PYWMM_RELRADPOW_H
#define PYWMM_RELRADPOW_H

#include "pywmm_aux.h"
#include "sph_harm.h"

/* Earth radius in km */
#define RADIUS  6371.2

#ifndef STR
#define SRINGIFY(x) #x
#define STR(x) SRINGIFY(x)
#endif

/* python function definition */

#define DOC_RELRADPOW "\n"\
"   rrp = relradpow(radius, degree, earth_radius="STR(RADIUS)")\n"\
"\n"\
"     For given 'radius' (geocentric spherical), evaluate relative radius power\n"\
"     series:\n"\
"       (earth_radius/radius)**(i+2) for i in range(0, degree+1) .\n"\


static PyObject* relradpow(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"latitude", "degree", "spherical", NULL};

    int degree;
    double rad, rad0 = RADIUS; // radius and reference radius
    PyObject *arr_rrp = NULL; // P array
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwdict,
            "di|d:relradpow", keywords, &rad, &degree, &rad0));

    if (degree < 1)
    {
        PyErr_Format(PyExc_ValueError, "Invalid value %d of '%s'!", degree, keywords[1]);
        goto exit;
    }

    if (rad < 0.0)
    {
        PyErr_Format(PyExc_ValueError, "Invalid value %g of '%s'!", rad, keywords[0]);
        goto exit;
    }

    if (rad0 <= 0.0)
    {
        PyErr_Format(PyExc_ValueError, "Invalid value %g of '%s'!", rad0, keywords[2]);
        goto exit;
    }

    // create the output array
    if (NULL == (arr_rrp = _get_new_double_array(1, NULL, degree+1)))
        goto exit;

    // evaluate the relative radius power
    rel_rad_pow(PyArray_DATA(arr_rrp), degree, rad/rad0);

    retval = arr_rrp;

  exit:

    // decrease reference counters to the arrays
    if (!retval && arr_rrp){Py_DECREF(arr_rrp);}

    return retval;
}

#endif  /* PYWMM_RELRADPOW_H */
