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

#ifndef PYWMM_LEGENDRE_H
#define PYWMM_LEGENDRE_H

#include "GeomagnetismHeader.h"
#include "pywmm_aux.h"
#include "pywmm_cconv.h"


/* python function definition */

#define DOC_LEGENDRE "\n"\
"   p, dp = legendre(latitude, degree, spherical=True)\n"\
"\n"\
"     For given 'latitude' and model's 'degree', evaluate asociative Legendre\n"\
"     functions. The input parameters are:\n"\
"       latitude - spherical (or geodetic) latitude in dg. of the evaluated\n"\
"                  location.\n"\
"       degree - degree of the spherical harmonic model.\n"\
"       spherical - boolean flag indicating whether a geodentic spherical\n"\
"                   (default, True) or geodetic (WGS84, False) latitude is\n"\
"                   being used.\n"


static PyObject* legendre(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"latitude", "degree", "spherical", NULL};

    int degree, nterm;
    double lat_sph, lat_in, tmp0, tmp1;
    int is_sph = 1;
    PyObject *arr_p = NULL; // P array
    PyObject *arr_dp = NULL; // dP array
    PyObject *retval = NULL; // output tuple

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwdict,
            "di|i:legendre", keywords, &lat_in, &degree, &is_sph));

    if (degree < 1)
    {
        PyErr_Format(PyExc_ValueError, "Invalid value %d of '%s'!", degree, keywords[1]);
        goto exit;
    }

    nterm = ((degree+1)*(degree+2))/2;

    // convert latitude to geocentric spherical latitude
    if (is_sph)
        lat_sph = lat_in;
    else
        conv_WGS84_to_sphECEF(&lat_sph, &tmp0, &tmp1, lat_in, 0.0, 0.0);

    // create the output arrays
    if (NULL == (arr_p = _get_new_double_array(1, NULL, nterm)))
        goto exit;

    if (NULL == (arr_dp = _get_new_double_array(1, NULL, nterm)))
        goto exit;

    if (NULL == (retval = Py_BuildValue("NN", arr_p, arr_dp)))
        goto exit;

    {
        MAGtype_LegendreFunction legfcn;
        MAGtype_CoordSpherical scoord;
        legfcn.Pcup = (double*) PyArray_DATA(arr_p);
        legfcn.dPcup = (double*) PyArray_DATA(arr_dp);
        scoord.phig = lat_sph;

        MAG_AssociatedLegendreFunction(scoord, degree, &legfcn);
    }

  exit:

    // decrease reference counters to the arrays
    if (!retval && arr_p){Py_DECREF(arr_p);}
    if (!retval && arr_dp){Py_DECREF(arr_dp);}

    return retval;
}

#endif  /* PYWMM_LEGENDRE_H */

