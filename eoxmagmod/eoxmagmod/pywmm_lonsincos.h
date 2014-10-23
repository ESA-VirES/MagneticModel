/*-----------------------------------------------------------------------------
 *
 * World Magnetic Model - C python bindings
 * - asociative Legendre functions evaluation
 *
 * Project: World Magnetic Model - python interface
 * Author: Martin Paces <martin.paces@eox.at>
 *
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

#ifndef PYWMM_LONSINCOS_H
#define PYWMM_LONSINCOS_H

#include "GeomagnetismHeader.h"
#include "pywmm_coord.h"
#include "pywmm_cconv.h"


/* python function definition */

#define DOC_LONSINCOS "\n"\
"  lonsin, loncos = lonsincos(longitude, degree, fast_alg=True)\n"\
"\n"\
"     For given 'longitude' (geocentric spherical) and 'degree', evaluate\n"\
"     the following sine and cosine series: \n"\
"        cos(i*longitude) for i in range(0, degree+1)\n"\
"        sin(i*longitude) for i in range(0, degree+1)\n"\
"     The longitude has to be entered in dg..\n"\
"     The 'fast_alg' boolean options forces the subroutine to use a faster\n"\
"     but sligtly less precise evaluation algorithm.\n"

static PyObject* lonsincos(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"longitude", "degree", "fast_alg", NULL};

    int degree;
    int fast_alg = 1;
    double lon_dg;
    PyObject *arr_lonsin = NULL; // lonsin array
    PyObject *arr_loncos = NULL; // loncos array
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwdict,
            "di|i:lonsincos", keywords, &lon_dg, &degree, &fast_alg));

    if (degree < 1)
    {
        PyErr_Format(PyExc_ValueError, "Invalid value %d of '%s'!", degree, keywords[1]);
        goto exit;
    }

    // create the output arrays
    if (NULL == (arr_lonsin = _get_new_double_array(1, NULL, degree+1)))
        goto exit;

    if (NULL == (arr_loncos = _get_new_double_array(1, NULL, degree+1)))
        goto exit;

    if (NULL == (retval = Py_BuildValue("NN", arr_lonsin, arr_loncos)))
        goto exit;

    {
        int i;
        double *lonsin = (double*)PyArray_DATA(arr_lonsin);
        double *loncos = (double*)PyArray_DATA(arr_loncos);
        double lon_rad = DG2RAD*lon_dg;
        double sin_lon = sin(lon_rad);
        double cos_lon = cos(lon_rad);
        double sl, sl_last, cl, cl_last;

        lonsin[0] = 0.0;
        loncos[0] = 1.0;
        lonsin[1] = sl_last = sin_lon;
        loncos[1] = cl_last = cos_lon;

        if (fast_alg)
        {
            // Faster evaluation based on pure recurrent
            // addition/substration and multiplication:
            //  sin(a + b) = cos(a)*sin(b) + sin(a)*cos(b)
            //  cos(a + b) = cos(a)*cos(b) - sin(a)*sin(b)
            for (i = 2; i <= degree; ++i)
            {
                lonsin[i] = sl = cl_last*sin_lon + sl_last*cos_lon;
                loncos[i] = cl = cl_last*cos_lon - sl_last*sin_lon;
                sl_last = sl;
                cl_last = cl;
            }
        }
        else
        {
            // Slower evaluation calling sin/cos for each term.
            for (i = 2; i <= degree; ++i)
            {
                lonsin[i] = sin(i*lon_rad);
                loncos[i] = cos(i*lon_rad);
            }
        }
    }

  exit:

    // decrease reference counters to the arrays
    if (!retval && arr_lonsin){Py_DECREF(arr_lonsin);}
    if (!retval && arr_loncos){Py_DECREF(arr_loncos);}

    return retval;
}

#endif  /* PYWMM_LONSINCOS_H */
