/*-----------------------------------------------------------------------------
 *
 * World Magnetic Model - C python bindings - magnetic model evaluation
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

#ifndef PYWMM_GEOMAG_H
#define PYWMM_GEOMAG_H

#include "GeomagnetismHeader.h"
#include "pywmm_coord.h"
#include "pywmm_cconv.h"

typedef struct {
    int degree;
    int nterm;
    int coord;
    double *coef_g;
    double *coef_h;
} MODEL;

/* single point evaluation */

static void _geomag_eval(double *fx, double *fy, double *fz,
                double x, double y, double z, const MODEL *model)
{
    double glat, glon, gh;
    double clat, clon, cr;

    // convert coordinates
    switch(model->coord)
    {
        case CT_GEODETIC_ABOVE_WGS84:
            glat = x; glon = y; gh = z;
            conv_WGS84_to_sphECEF(&clat, &clon, &cr, glat, glon, gh);
            break;
        case CT_GEODETIC_ABOVE_EGM96:
            conv_EGM96_to_WGS84(&glat, &glon, &gh, x, y, z);
            conv_WGS84_to_sphECEF(&clat, &clon, &cr, glat, glon, gh);
            break;
        case CT_GEOCENTRIC_SPHERICAL:
            clat = x; clon = y; cr = z;
            conv_sphECEF_to_WGS84(&glat, &glon, &gh, clat, clon, cr);
            break;
        case CT_GEOCENTRIC_CARTESIAN:
            conv_cartECEF_to_sphECEF(&clat, &clon, &cr, x,  y, z);
            conv_sphECEF_to_WGS84(&glat, &glon, &gh, clat, clon, cr);
        default:
            return;
    }

    // model evaluation
    {
        MAGtype_LegendreFunction *LegendreFunction = MAG_AllocateLegendreFunctionMemory(model->nterm);
        MAGtype_SphericalHarmonicVariables *SphVariables = MAG_AllocateSphVarMemory(model->degree);

        MAGtype_Ellipsoid Ellip;
        MAGtype_Geoid Geoid; // not used
        MAG_SetDefaults(&Ellip, &Geoid);

        MAGtype_MagneticModel MagneticModel;

        MAGtype_CoordGeodetic CoordGeodetic;
        MAGtype_CoordSpherical CoordSpherical;

        MAGtype_MagneticResults MagneticResultsSph, MagneticResultsGeo;

        MagneticModel.nMax = model->degree;
        MagneticModel.Main_Field_Coeff_G = model->coef_g;
        MagneticModel.Main_Field_Coeff_H = model->coef_h;

        CoordGeodetic.lambda = glon;
        CoordGeodetic.phi = glat;
        CoordGeodetic.HeightAboveEllipsoid = gh;

        CoordSpherical.lambda = clon;
        CoordSpherical.phig = clat;
        CoordSpherical.r = cr;

        MAG_ComputeSphericalHarmonicVariables(Ellip, CoordSpherical, model->degree, SphVariables);
        MAG_AssociatedLegendreFunction(CoordSpherical, model->degree, LegendreFunction);

        // Accumulate the spherical harmonic coefficients
        MAG_Summation(LegendreFunction, &MagneticModel, *SphVariables, CoordSpherical, &MagneticResultsSph);

        // Map the computed Magnetic fields to Geodeitic coordinates
        MAG_RotateMagneticVector(CoordSpherical, CoordGeodetic, MagneticResultsSph, &MagneticResultsGeo);

        MAG_FreeLegendreMemory(LegendreFunction);
        MAG_FreeSphVarMemory(SphVariables);

        *fx = MagneticResultsGeo.Bx;
        *fy = MagneticResultsGeo.By;
        *fz = MagneticResultsGeo.Bz;
    }
}

/* recursive model_evaluation */

static void _geomag(ARRAY_DATA arrd_in, ARRAY_DATA arrd_out, MODEL *model)
{
    if (arrd_in.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < arrd_in.dim[0]; ++i)
            _geomag(_get_arrd_item(&arrd_in, i), _get_arrd_item(&arrd_out, i), model);
        return;
    }

    #define P(a,i) ((double*)((a).data+(i)*(a).stride[0]))
    #define V(a,i) (*P(a,i))

    _geomag_eval(P(arrd_out, 0), P(arrd_out, 1), P(arrd_out, 2),
              V(arrd_in, 0), V(arrd_in, 1), V(arrd_in, 2), model);
/*
    fconv(P(arrd_out, 0), P(arrd_out, 1), P(arrd_out, 2),
          V(arrd_in, 0), V(arrd_in, 1), V(arrd_in, 2));
*/
    #undef V
    #undef P
}

/* python function definition */

#define DOC_GEOMAG "\n"\
"   arr_out = geomag_static(arr_in, degree, coef_g, coef_h, coord_type=GEODETIC_ABOVE_WGS84)\n"\
"\n     Evaluate WMM model for a given array of 3D location vectors.\n"\
"     Parameters:\n"\
"        arr_in - array of 3D coordinates (up to 16 dimensions).\n"\
"        degree - degree of the spherical harmonic model.\n"\
"        coef_g - vector of spherical harmonic model coeficients.\n"\
"        coef_h - vector of spherical harmonic model coeficients.\n"\
"        cood_type - type of the input coordinates.\n"


static PyObject* geomag(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {"arr_in", "degree", "coef_g", "coef_h",
                               "coord_type", NULL };
    int ct_in = CT_GEODETIC_ABOVE_WGS84;
    int nterm, degree = 0;
    PyObject *obj_in = NULL; // input object
    PyObject *obj_cg = NULL; // coef_g object
    PyObject *obj_ch = NULL; // coef_h object
    PyObject *arr_in = NULL; // input array
    PyObject *arr_cg = NULL; // coef_g array
    PyObject *arr_ch = NULL; // coef_h array
    PyObject *arr_out = NULL; // output array
    PyObject *retval = NULL;

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(args, kwdict,
            "OiOO|i:geomag", keywords, &obj_in, &degree, &obj_cg, &obj_ch, &ct_in))
        goto exit;

    // check the type of the coordinate transformation
    if (CT_INVALID == _check_coord_type(ct_in, keywords[4])) goto exit;

    // cast the objects to arrays
    if (NULL == (arr_in=_get_as_double_array(obj_in, 1, 0, NPY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_cg=_get_as_double_array(obj_cg, 1, 1, NPY_IN_ARRAY, keywords[2])))
        goto exit;

    if (NULL == (arr_ch=_get_as_double_array(obj_ch, 1, 1, NPY_IN_ARRAY, keywords[3])))
        goto exit;

    if (degree < 1)
    {
        PyErr_Format(PyExc_ValueError, "Invalid value %d of '%s'!", degree, keywords[1]);
        goto exit;
    }

    // check dimensions of the coeficient arrays
    nterm = ((degree+1)*(degree+2))/2;

    if (PyArray_DIM(arr_cg, 0) != nterm)
    {
        PyErr_Format(PyExc_ValueError, "The size of '%s'"\
            " %d does not match the required value %d!", keywords[2],
            (int)PyArray_DIM(arr_cg, 0), nterm);
        goto exit;
    }

    if (PyArray_DIM(arr_ch, 0) != nterm)
    {
        PyErr_Format(PyExc_ValueError, "The size of '%s'"\
            " %d does not match the required value %d!", keywords[3],
            (int)PyArray_DIM(arr_ch, 0), nterm);
        goto exit;
    }

    // check maximum allowed input array dimension
    if (PyArray_NDIM(arr_in) > MAX_OUT_ARRAY_NDIM)
    {
        PyErr_Format(PyExc_ValueError, "The input array dimension of '%s'"\
            " %d exceeds the allowed maximum value %d!", keywords[0],
            PyArray_NDIM(arr_in), MAX_OUT_ARRAY_NDIM);
        goto exit;
    }

    // check the last dimension (required length of the coordinates vectors)
    if (_check_last_array_dimension(arr_in, 3, keywords[0]))
        goto exit;

    // create the output array
    if (NULL == (arr_out = _get_new_double_array(PyArray_NDIM(arr_in), PyArray_DIMS(arr_in), 3)))
        goto exit;

    // process
    //_convert(_array_to_arrd(arr_in), _array_to_arrd(arr_out), _get_fconv(ct_in, ct_out));
    {
        MODEL model;
        model.degree = degree;
        model.nterm = nterm;
        model.coord = ct_in;
        model.coef_g = (double*) PyArray_DATA(arr_cg);
        model.coef_h = (double*) PyArray_DATA(arr_ch);

        _geomag(_array_to_arrd(arr_in), _array_to_arrd(arr_out), &model);
    }
    //retval = Py_None; Py_INCREF(Py_None);
    retval = arr_out;

  exit:

    // decrease reference counters to the arrays
    if (arr_in){Py_DECREF(arr_in);}
    if (arr_cg){Py_DECREF(arr_cg);}
    if (arr_ch){Py_DECREF(arr_ch);}
    if (!retval && arr_out){Py_DECREF(arr_out);}

    return retval;
}

#endif  /* PYWMM_GEOMAG_H */
