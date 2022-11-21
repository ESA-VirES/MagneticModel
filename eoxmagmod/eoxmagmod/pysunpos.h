/*-----------------------------------------------------------------------------
 *
 * Solar position
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2017-2022 EOX IT Services GmbH
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
#ifndef PYSUNPOS_H
#define PYSUNPOS_H

#include "pymm_aux.h"
#include "sun_ephemeris.h"
#include <math.h>

/*
 *  nD-array recursive evaluation
 */

#ifndef RAD2DEG
#define RAD2DEG (180.0/M_PI)
#endif

#ifndef DEG2RAD
#define DEG2RAD (M_PI/180.0)
#endif

static void _sunpos(
    ARRAY_DATA arrd_dim,
    ARRAY_DATA arrd_mjd,
    ARRAY_DATA arrd_lat,
    ARRAY_DATA arrd_lon,
    ARRAY_DATA arrd_rad,
    ARRAY_DATA arrd_dtt,
    ARRAY_DATA arrd_out
)
{
    if (arrd_dim.ndim > 0)
    {
        npy_intp i;
        for(i = 0; i < arrd_dim.dim[0]; ++i)
            _sunpos(
                _get_arrd_item(&arrd_dim, i),
                _get_arrd_item(&arrd_mjd, i),
                _get_arrd_item(&arrd_lat, i),
                _get_arrd_item(&arrd_lon, i),
                _get_arrd_item(&arrd_rad, i),
                _get_arrd_item(&arrd_dtt, i),
                _get_arrd_item(&arrd_out, i)
            );
    }
    else
    {
        #define P(a,i) ((double*)((a).data+(i)*(a).stride[0]))
        #define V(a,i) (*P(a,i))

        double rasc, decl, hang, zenith, azimuth;
        double mjd2k, dtt, lat, lon, rad, pres, temp;

        mjd2k = V(arrd_mjd,0); // decimal days since 2000-01-01T00:00:00
        dtt = V(arrd_dtt,0); // seconds (TT to UT offset)
        lat = V(arrd_lat,0) * DEG2RAD; // deg --> rad
        lon = V(arrd_lon,0) * DEG2RAD; // deg --> rad
        rad = V(arrd_rad,0); // km
        pres = 0; // atm
        temp = 0; // dgC

        sunpos5equat(&decl, &rasc, &hang, mjd2k, dtt, lon);
        sunpos5eq2hor(&azimuth, &zenith, decl, hang, lat, rad, pres, temp);

        V(arrd_out,0) = RAD2DEG * decl; // rad --> deg
        V(arrd_out,1) = RAD2DEG * rasc; // rad --> deg
        V(arrd_out,2) = RAD2DEG * hang; // rad --> deg
        V(arrd_out,3) = RAD2DEG * azimuth; // rad --> deg
        V(arrd_out,4) = RAD2DEG * zenith; // rad --> deg

        #undef V
        #undef P
    }
}


/*
 * Python function definition
 */

#define DOC_SUNPOS "\n"\
"   arr_out = sunpos(time_mjd2k, lat, lon, rad, dtt)\n\n"\
"     Output:\n"\
"       arr_out - array of the Sun equatorial and horizontal coordinates:\n"\
"                 - declination\n"\
"                 - right ascension\n"\
"                 - hour angle\n"\
"                 - azimuth\n"\
"                 - zenith\n"\
"                 All angles are in deg.\n\n"\
"     Parameters:\n"\
"       time_mjd2k - array of MJD2000 times (up to 15 dimensions).\n"\
"       lat - array of latitudes in deg.\n"\
"       lon - array of longitudes in deg.\n"\
"       rad - array of radii in km (0 means no parallax correction).\n"\
"       dtt - array of offsets to TT in sec. \n"

static PyObject* pysunpos_sunpos(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"time_mjd2k", "lat", "lon", "rad", "dtt", NULL};
    PyObject *obj_mjd = NULL;
    PyObject *obj_lat = NULL;
    PyObject *obj_lon = NULL;
    PyObject *obj_rad = NULL;
    PyObject *obj_dtt = NULL;
    PyArrayObject *arr_mjd = NULL;
    PyArrayObject *arr_lat = NULL;
    PyArrayObject *arr_lon = NULL;
    PyArrayObject *arr_rad = NULL;
    PyArrayObject *arr_dtt = NULL;
    PyArrayObject *arr_out = NULL; // output array
    PyArrayObject *arr = NULL;
    PyObject *retval = NULL; // return object
    int idx = 0;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOOOO|:sunpos", keywords,
        &obj_mjd, &obj_lat, &obj_lon, &obj_rad, &obj_dtt
    ))
        goto exit;

    // cast the objects to arrays
    if (NULL == (arr_mjd = _get_as_double_array(obj_mjd, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_lat = _get_as_double_array(obj_lat, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_lon = _get_as_double_array(obj_lon, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_rad = _get_as_double_array(obj_rad, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_dtt = _get_as_double_array(obj_dtt, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    // check the array dimensions and allocate the output array
    // get array with the maximum dimension
    arr = arr_mjd;
    idx = 0;
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_lat)) {idx = 1; arr = arr_lat;}
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_lon)) {idx = 2; arr = arr_lon;}
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_rad)) {idx = 3; arr = arr_rad;}
    if (PyArray_NDIM(arr) < PyArray_NDIM(arr_dtt)) {idx = 4; arr = arr_dtt;}

    // check maximum allowed input array dimension
    if (PyArray_NDIM(arr) > (MAX_OUT_ARRAY_NDIM-1))
    {
        PyErr_Format(PyExc_ValueError, "Array dimension of '%s'"\
            " %d exceeds the allowed maximum value %d!", keywords[idx],
            PyArray_NDIM(arr), (MAX_OUT_ARRAY_NDIM-1));
        goto exit;
    }

    // check array dimensions
    if (PyArray_NDIM(arr_mjd) > 0)
        if (_check_equal_shape(arr_mjd, arr, keywords[0], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_lat) > 0)
        if (_check_equal_shape(arr_lat, arr, keywords[1], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_lon) > 0)
        if (_check_equal_shape(arr_lon, arr, keywords[2], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_rad) > 0)
        if (_check_equal_shape(arr_rad, arr, keywords[3], keywords[idx]))
            goto exit;

    if (PyArray_NDIM(arr_dtt) > 0)
        if (_check_equal_shape(arr_dtt, arr, keywords[4], keywords[idx]))
            goto exit;

    // create a new output array (adding one extra dimension)
    {
        npy_intp *dims_src = PyArray_DIMS(arr);
        npy_intp dims_dst[MAX_OUT_ARRAY_NDIM];
        int i, ndim_src = PyArray_NDIM(arr);

        for (i = 0; i < ndim_src; ++i) dims_dst[i] = dims_src[i];
        dims_dst[ndim_src] = 5;

        if (NULL == (arr_out = (PyArrayObject*) PyArray_EMPTY(ndim_src + 1, dims_dst, NPY_DOUBLE, 0)))
            goto exit;
    }

    // evaluate solar positions
    _sunpos(
        _array_to_arrd(arr),
        _array_to_arrd(arr_mjd),
        _array_to_arrd(arr_lat),
        _array_to_arrd(arr_lon),
        _array_to_arrd(arr_rad),
        _array_to_arrd(arr_dtt),
        _array_to_arrd(arr_out)
    );

    // assign return value
    retval = (PyObject*) arr_out;

  exit:

    // decrease reference counters to the arrays
    if (arr_mjd) Py_DECREF(arr_mjd);
    if (arr_lat) Py_DECREF(arr_lat);
    if (arr_lon) Py_DECREF(arr_lon);
    if (arr_rad) Py_DECREF(arr_rad);
    if (arr_dtt) Py_DECREF(arr_dtt);
    if (!retval && arr_out) Py_DECREF(arr_out);

    return retval;
}

#endif /* PYSUNPOS_H */
