/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - coordinate conversions
 *
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

#ifndef PYWMM_CCONV_H
#define PYWMM_CCONV_H

#include "wmm/GeomagnetismHeader.h"
#include "wmm/EGM9615.h"
#include "geo_conv.h"
#include "pywmm_aux.h"
#include "pywmm_coord.h"

#ifndef NAN
#define NAN (0.0/0.0)
#endif

/* NOTE: Geomagnetism library is used to get the geoid height */
static double delta_EGM96_to_WGS84(double lat, double lon)
{
    double delta_h = NAN;
    MAGtype_Geoid Geoid;

    /* Sets EGM-96 model file parameters (see also  EGM9615.h)*/
    Geoid.GeoidHeightBuffer = GeoidHeightBuffer; /* see EGM9615.h */
    Geoid.NumbGeoidCols = 1441; /* 360 degrees of longitude at 15 minute spacing */
    Geoid.NumbGeoidRows = 721; /* 180 degrees of latitude  at 15 minute spacing */
    Geoid.NumbHeaderItems = 6; /* min, max lat, min, max long, lat, long spacing*/
    Geoid.ScaleFactor = 4; /* 4 grid cells per degree at 15 minute spacing  */
    Geoid.NumbGeoidElevs = Geoid.NumbGeoidCols * Geoid.NumbGeoidRows;
    Geoid.Geoid_Initialized = 1;
    Geoid.UseGeoid = 1;

    MAG_GetGeoidHeight(lat, lon, &delta_h, &Geoid);

    return delta_h * 1e-3; /* scale to km */
}

/* low-level conversion handlers */

typedef void (*f_conv)(double*, double*, double*, double, double, double);

static void conv_identity(double *x1, double *y1, double *z1,
                        double x0, double y0, double z0)
{
    *x1 = x0; *y1 = y0; *z1 = z0;
}

static void conv_EGM96_to_WGS84(double *lat1, double *lon1, double *h1,
                        double lat0, double lon0, double h0)
{
    *lat1 = lat0; *lon1 = lon0; *h1 = h0 + delta_EGM96_to_WGS84(lat0, lon0);
}

static void conv_WGS84_to_EGM96(double *lat1, double *lon1, double *h1,
                        double lat0, double lon0, double h0)
{
    *lat1 = lat0; *lon1 = lon0; *h1 = h0 - delta_EGM96_to_WGS84(lat0, lon0);
}

static void conv_WGS84_to_cartECEF(double *x, double *y, double *z,
                        double lat, double lon, double h)
{
    geodetic2geocentric_cart(x, y, z, lat, lon, h, WGS84_A, WGS84_EPS2);
}

static void conv_WGS84_to_sphECEF(double *lat1, double *lon1, double *r1,
                        double lat0, double lon0, double h0)
{
    geodetic2geocentric_sph(r1, lat1, lon1, lat0, lon0, h0, WGS84_A, WGS84_EPS2);
    *lat1 *= RAD2DG;
    *lon1 *= RAD2DG;
}

static void conv_EGM96_to_cartECEF(double *x, double *y, double *z,
                        double lat, double lon, double h)
{
    double tmp_lat, tmp_lon, tmp_h;
    conv_EGM96_to_WGS84(&tmp_lat, &tmp_lon, &tmp_h, lat, lon, h);
    conv_WGS84_to_cartECEF(x, y, z, lat, lon, tmp_h);
}

static void conv_EGM96_to_sphECEF(double *lat1, double *lon1, double *r1,
                        double lat0, double lon0, double h0)
{
    double tmp_lat, tmp_lon, tmp_h;
    conv_EGM96_to_WGS84(&tmp_lat, &tmp_lon, &tmp_h, lat0, lon0, h0);
    conv_WGS84_to_sphECEF(lat1, lon1, r1, lat0, lon0, tmp_h);
}

static void conv_cartECEF_to_sphECEF(double *lat, double *lon, double *r,
                         double x, double y, double z)
{
    cart2sph(r, lat, lon, x, y, z);
    *lat *= RAD2DG;
    *lon *= RAD2DG;
}

static void conv_cartECEF_to_WGS84(double *lat, double *lon, double *h,
                         double x, double y, double z)
{
    geocentric_cart2geodetic(lat, lon, h, x, y, z, WGS84_A, WGS84_EPS2);
}


static void conv_cartECEF_to_EGM96(double *lat, double *lon, double *h,
                         double x, double y, double z)
{
    double tmp_lat, tmp_lon, tmp_h;
    conv_cartECEF_to_WGS84(&tmp_lat, &tmp_lon, &tmp_h, x, y, z);
    conv_WGS84_to_EGM96(lat, lon, h, tmp_lat, tmp_lon, tmp_h);

}

static void conv_sphECEF_to_cartECEF(double *x, double *y, double *z,
                         double lat, double lon, double r)
{
    sph2cart(x, y, z, r, DG2RAD*lat, DG2RAD*lon);
}

static void conv_sphECEF_to_WGS84(double *lat1, double *lon1, double *h1,
                           double lat0, double lon0, double r0)
{
    geocentric_sph2geodetic(lat1, lon1, h1, r0, DG2RAD*lat0, DG2RAD*lon0, WGS84_A, WGS84_EPS2);
}

static void conv_sphECEF_to_EGM96(double *lat1, double *lon1, double *h1,
                           double lat0, double lon0, double r0)
{
    double tmp_lat, tmp_lon, tmp_h;
    conv_sphECEF_to_WGS84(&tmp_lat, &tmp_lon, &tmp_h, lat0, lon0, r0);
    conv_WGS84_to_EGM96(lat1, lon1, h1, tmp_lat, tmp_lon, tmp_h);
}

static f_conv _get_fconv(int ct_in, int ct_out)
{
    f_conv mtx[4][4] = {
        {conv_identity, conv_WGS84_to_EGM96, conv_WGS84_to_sphECEF, conv_WGS84_to_cartECEF},
        {conv_EGM96_to_WGS84, conv_identity, conv_EGM96_to_sphECEF, conv_EGM96_to_cartECEF},
        {conv_sphECEF_to_WGS84, conv_sphECEF_to_EGM96, conv_identity, conv_sphECEF_to_cartECEF},
        {conv_cartECEF_to_WGS84, conv_cartECEF_to_EGM96, conv_cartECEF_to_sphECEF, conv_identity}
    };

    return mtx[ct_in][ct_out];
}

/* recursive coordinate conversion */

static void _convert(ARRAY_DATA arrd_in, ARRAY_DATA arrd_out, f_conv fconv)
{
    if (arrd_in.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < arrd_in.dim[0]; ++i)
            _convert(_get_arrd_item(&arrd_in, i), _get_arrd_item(&arrd_out, i), fconv);
    }
    else
    {
        #define P(a,i) ((double*)((a).data+(i)*(a).stride[0]))
        #define V(a,i) (*P(a,i))

        fconv(P(arrd_out, 0), P(arrd_out, 1), P(arrd_out, 2),
              V(arrd_in, 0), V(arrd_in, 1), V(arrd_in, 2));

        #undef V
        #undef P
    }
}


/* python function definition */

#define DOC_CONVERT "\n"\
"   arr_out = convert(arr_in, coord_type_in=GEODETIC_ABOVE_WGS84, coord_type_out=GEODETIC_ABOVE_WGS84)\n"\
"\n     Convert coordinates between different coordinate systems.\n"


static PyObject* convert(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"arr_in", "coord_type_in", "coord_type_out",NULL};
    int ct_in = CT_GEODETIC_ABOVE_WGS84;
    int ct_out = CT_GEODETIC_ABOVE_WGS84;
    PyObject *obj_in = NULL; // input object
    PyObject *arr_in = NULL; // input array
    PyObject *arr_out = NULL; // output array
    PyObject *retval = NULL;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "O|ii:convert", keywords, &obj_in, &ct_in, &ct_out
    ))
        goto exit;

    // check the type of the coordinate transformation
    if (CT_INVALID == _check_coord_type(ct_in, keywords[0])) goto exit;
    if (CT_INVALID == _check_coord_type(ct_out, keywords[1])) goto exit;

    // cast the object to an array
    if (NULL == (arr_in=_get_as_double_array(obj_in, 1, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    // check maximum allowed input array dimension
    if (PyArray_NDIM(arr_in) > MAX_OUT_ARRAY_NDIM)
    {
        PyErr_Format(PyExc_ValueError, "The input array dimension of '%s'"\
            " %d exceeds the allowed maximum value %d!", keywords[0],
            PyArray_NDIM(arr_in), MAX_OUT_ARRAY_NDIM);
        goto exit;
    }

    // check the last dimension (required length of the coordinates vectors)
    if (_check_array_dim_eq(arr_in, -1, 3, keywords[0]))
        goto exit;

    // fast-track identity
    if (ct_in == ct_out) return arr_in;

    // create the output array
    if (NULL == (arr_out = _get_new_double_array(PyArray_NDIM(arr_in), PyArray_DIMS(arr_in), 3)))
        goto exit;

    // process
    _convert(_array_to_arrd(arr_in), _array_to_arrd(arr_out), _get_fconv(ct_in, ct_out));

    //retval = Py_None; Py_INCREF(Py_None);
    retval = arr_out;

  exit:

    // decrease reference counters to the arrays
    if (arr_in){Py_DECREF(arr_in);}
    if (!retval && arr_out){Py_DECREF(arr_out);}

    return retval;
}

#endif  /* PYWMM_CCONV_H */
