/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - spherical harmonic model evaluation
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

#ifndef PYMM_SHEVAL_H
#define PYMM_SHEVAL_H

#include <math.h>
#include <string.h>
#include <stdlib.h>

#include "math_aux.h"
#include "geo_conversion.h"
#include "fourier_series.h"
#include "spherical_harmonics.h"
#include "pymm_aux.h"
#include "pymm_coord.h"
#include "pymm_cconv.h"
#include "pymm_sphar_common.h"

#ifndef NAN
#define NAN (0.0/0.0)
#endif

/* Earth radius in km */
#ifndef RADIUS
#define RADIUS  6371.2
#endif

#ifndef STR
#define SRINGIFY(x) #x
#define STR(x) SRINGIFY(x)
#endif

/* modes of evaluation - enumerate */
typedef enum {
    SM_INVALID = 0x0,
    SM_POTENTIAL = 0x1,
    SM_GRADIENT = 0x2,
    SM_POTENTIAL_AND_GRADIENT = 0x3,
} SHEVAL_MODE;

static SHEVAL_MODE _check_sheval_mode(int mode, const char *label);


/* magnetic model - auxiliary structure */
typedef struct {
    int is_internal;
    int degree;
    int nterm;
    int coord_in;
    int coord_out;
    double elps_a;
    double elps_eps2;
    double crad_ref;
    double clat_last;
    double clon_last;
    double crad_last;
    double scale_potential;
    double scale_gradient0;
    double scale_gradient1;
    double scale_gradient2;
    double *lp;
    double *ldp;
    double (*lcs)[2];
    double *rrp;
    double *psqrt;
    double *prsqrt;
} MODEL;

static void _model_reset(MODEL *model);
static void _model_destroy(MODEL *model);
static int _model_init(
    MODEL *model, int is_internal, int degree, int coord_in, int coord_out,
    const double scale_potential, const double *scale_gradient
);

static void _sheval_pot(
    ARRAY_DATA arrd_pot,
    ARRAY_DATA arrd_x, ARRAY_DATA arrd_coef, MODEL *model
);
static void _sheval_grd(
    ARRAY_DATA arrd_grd,
    ARRAY_DATA arrd_x, ARRAY_DATA arrd_coef, MODEL *model
);
static void _sheval_pot_and_grd(
    ARRAY_DATA arrd_pot, ARRAY_DATA arrd_grd,
    ARRAY_DATA arrd_x, ARRAY_DATA arrd_coef, MODEL *model
);

/* Python function definition */

#define DOC_SHEVAL "\n"\
"   arr_out = sheval(arr_x, coef, coord_type_in=GEODETIC_ABOVE_WGS84,\n"\
"                    coord_type_out=GEODETIC_ABOVE_WGS84, mode=GRADIENT,\n"\
"                    is_internal=True, degree=-1,\n"\
"                    scale_potential=1.0, scale_gradient=1.0)\n"\
"\n"\
"     Parameters:\n"\
"       arr_x  - array of 3D coordinates.\n"\
"       coef - array of spherical harmonic model coefficients.\n"\
"       coord_type_in - type of the input coordinates.\n"\
"       coord_type_out - type of the output coordinates frame.\n"\
"       mode - quantity to be evaluated:\n"\
"                  POTENTIAL\n"\
"                  GRADIENT (default)\n"\
"                  POTENTIAL_AND_GRADIENT\n"\
"       is_internal - boolean flag set to True by default. When set to False \n"\
"                     external field evaluation is used.\n"\
"       degree - degree of the spherical harmonic model. If not provided or \n"\
"                set to a negative number, then it is derived from the array \n"\
"                sizes.\n"\
"       scale_potential - optional scalar multiplied with the result potentials.\n"\
"       scale_gradient - optional scalar or 3 element array multiplied with\n"\
"                        the result gradient components.\n"\
"\n"


static PyObject* sheval(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "arr_x", "coef", "coord_type_in", "coord_type_out",
        "mode", "is_internal", "degree", "scale_potential", "scale_gradient",
        NULL
    };
    int ct_in = CT_GEODETIC_ABOVE_WGS84;
    int ct_out = CT_GEODETIC_ABOVE_WGS84;
    int degree = -1;
    int mode = SM_GRADIENT;
    int is_internal;
    double scale_potential = 1.0;
    double scale_gradient[3] = {1.0, 1.0, 1.0};
    PyObject *obj_is_internal = NULL; // boolean flag
    PyObject *obj_x = NULL; // input object
    PyObject *obj_coef = NULL; // coefficients object
    PyObject *obj_scale = NULL; // gradient scale object
    PyArrayObject *arr_x = NULL; // coordinates array
    PyArrayObject *arr_coef = NULL; // coefficients array
    PyArrayObject *arr_pot = NULL; // output array
    PyArrayObject *arr_grd = NULL; // output array
    PyArrayObject *arr_scale = NULL; // gradient scale array
    PyObject *retval = NULL;
    MODEL model = {0};

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OO|iiiOidO:sheval", keywords,
        &obj_x, &obj_coef, &ct_in, &ct_out, &mode,
        &obj_is_internal, &degree, &scale_potential, &obj_scale
    ))
        goto exit;

    is_internal = (obj_is_internal == NULL) || PyObject_IsTrue(obj_is_internal);

    // check the type of the coordinate transformation
    if (CT_INVALID == _check_coord_type(ct_in, keywords[4]))
        goto exit;

    if (CT_INVALID == _check_coord_type(ct_out, keywords[5]))
        goto exit;

    // check the operation mode
    if (SM_INVALID == _check_sheval_mode(mode, keywords[6]))
        goto exit;

    // cast the objects to arrays
    if (NULL == (arr_x = _get_as_double_array(obj_x, 1, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_coef = _get_as_double_array(obj_coef, 2, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_IN_ARRAY, keywords[2])))
        goto exit;

    // extract degree from the array dimensions
    { 
        const int ndegrees = 1;
        npy_intp degrees[] = {
            _size_to_degree(PyArray_DIMS(arr_coef)[PyArray_NDIM(arr_coef)-2]),
        };

        npy_intp max_degree;
        int arg_idx;

         _get_max_degree(&arg_idx, &max_degree, degrees, ndegrees);

        if (max_degree < 0)
        {
            PyErr_Format(
                PyExc_ValueError,
                "Negative degree due to empty %s array!",
                keywords[arg_idx+1]
            );
            goto exit;
        }

        if ((max_degree > MAX_DEGREE)&&(degree < 0))
        {
            PyErr_Format(PyExc_ValueError, "Degree is larger than %g!", MAX_DEGREE);
            goto exit;
        }

        // set the applied degree
        if ((degree < 0)||(degree > max_degree))
            degree = max_degree;
    }

    // check array dimensions and allocate the output array
    {
        const int narr = 2;
        PyArrayObject *arr[] = {arr_x, arr_coef};
        npy_intp arr_ndim[] = {1, 2};

        int arg_idx;
        npy_intp ndim, *dims;

        // extract maximum common shape
        _extract_common_shape(&arg_idx, &ndim, &dims, arr, arr_ndim, narr);

        // check if arrays match the common shape
        {
            int i;

            for (i = 0; i < narr; ++i)
            {
                if (i != arg_idx)
                {
                    if (_compare_dimensions(dims, PyArray_DIMS(arr[i]), PyArray_NDIM(arr[i]) - arr_ndim[i]))
                    {
                        PyErr_Format(
                            PyExc_ValueError,
                            "Shape of %s array does not match shape of %s!",
                            keywords[i], keywords[arg_idx]
                        );
                        goto exit;
                    }
                }

                // check last dimension of the pair-arrays
                if ((arr_ndim[i] == 2)&&(PyArray_DIMS(arr[i])[PyArray_NDIM(arr[i])-1] != 2))
                {
                    PyErr_Format(PyExc_ValueError, "Invalid shape of %s array!", keywords[i]);
                    goto exit;
                }
            }
        }

        // allocate the output arrays using the common shape

        if (mode & SM_POTENTIAL)
        {
            if (NULL == (arr_pot = _get_new_array(ndim, dims, NPY_DOUBLE)))
                goto exit;
        }

        if (mode & SM_GRADIENT)
        {
            npy_intp j;
            npy_intp ndim_new = ndim + 1;
            npy_intp dims_new[ndim_new];

            for (j = 0; j < ndim; ++j)
                dims_new[j] = dims[j];
            dims_new[ndim] = 3;

            if (NULL == (arr_grd = _get_new_array(ndim_new, dims_new, NPY_DOUBLE)))
                goto exit;
        }

    }

    // handle the gradient scale factors
    if (NULL != obj_scale)
    {
        if (NULL == (arr_scale = _get_as_double_array(obj_scale, 0, 1, NPY_ARRAY_IN_ARRAY, keywords[9])))
            goto exit;

        if (_extract_1d_double_array(scale_gradient, 3, arr_scale, keywords[9]))
            goto exit;
    }

    // evaluate the spherical harmonics

    if(_model_init(
        &model, is_internal, degree, ct_in, ct_out,
        scale_potential, scale_gradient
    ))
        goto exit;

    switch(mode)
    {
        case SM_POTENTIAL:
            _sheval_pot(
                _array_to_arrd(arr_pot),
                _array_to_arrd(arr_x),
                _array_to_arrd(arr_coef),
                &model
            );
            retval = (PyObject*) arr_pot;
            break;

        case SM_GRADIENT:
            _sheval_grd(
                _array_to_arrd(arr_grd),
                _array_to_arrd(arr_x),
                _array_to_arrd(arr_coef),
                &model
            );
            retval = (PyObject*) arr_grd;
            break;

        case SM_POTENTIAL_AND_GRADIENT:
            _sheval_pot_and_grd(
                _array_to_arrd(arr_pot),
                _array_to_arrd(arr_grd),
                _array_to_arrd(arr_x),
                _array_to_arrd(arr_coef),
                &model
            );
            if (NULL == (retval = Py_BuildValue("NN", (PyObject*) arr_pot, (PyObject*) arr_grd)))
                goto exit;
            break;
    }

  exit:

    _model_destroy(&model);

    // decrease reference counters to the arrays
    if (arr_x) Py_DECREF(arr_x);
    if (arr_coef) Py_DECREF(arr_coef);
    if (arr_scale) Py_DECREF(arr_scale);
    if (!retval && arr_grd) Py_DECREF(arr_grd);
    if (!retval && arr_pot) Py_DECREF(arr_pot);

    return retval;
}


/* checking the mode of the model evaluation */

static SHEVAL_MODE _check_sheval_mode(int mode, const char *label)
{
    switch (mode)
    {
        case SM_POTENTIAL:
            return SM_POTENTIAL;
        case SM_GRADIENT:
            return SM_GRADIENT;
        case SM_POTENTIAL_AND_GRADIENT:
            return SM_POTENTIAL_AND_GRADIENT;
    }

    PyErr_Format(PyExc_ValueError, "Invalid mode value!");
    return SM_INVALID;
}


/* high-level nD-array recursive batch model_evaluation */

static void _model_eval(
    double *fpot, double *fx, double *fy, double *fz,
    const double x, const double y, const double z,
    const double (*coeff)[2], MODEL *model, const int mode
);


static void _sheval_pot(
    ARRAY_DATA arrd_pot,
    ARRAY_DATA arrd_x,
    ARRAY_DATA arrd_coef,
    MODEL *model
)
{
    if (arrd_pot.ndim > 0)
    {
        npy_intp i;
        for(i = 0; i < arrd_pot.dim[0]; ++i)
        {
            _sheval_pot(
                _get_arrd_item_nocheck(&arrd_pot, i),
                _get_arrd_item_with_guard(&arrd_x, i, 1),
                _get_arrd_item_with_guard(&arrd_coef, i, 2),
                model
            );
        }
        return;
    }
    else
    {
        double *pot = ((double*)arrd_pot.data);
        const double (*coeff)[2] = ((double(*)[2])arrd_coef.data);
        const double *x = ((double*)arrd_x.data);

        _model_eval(
            pot, NULL, NULL, NULL, x[0], x[1], x[2],
            coeff, model, SM_POTENTIAL
        );
    }
}


static void _sheval_grd(
    ARRAY_DATA arrd_grd,
    ARRAY_DATA arrd_x,
    ARRAY_DATA arrd_coef,
    MODEL *model
)
{
    if (arrd_grd.ndim > 1)
    {
        npy_intp i;
        for(i = 0; i < arrd_grd.dim[0]; ++i)
        {
            _sheval_grd(
                _get_arrd_item_nocheck(&arrd_grd, i),
                _get_arrd_item_with_guard(&arrd_x, i, 1),
                _get_arrd_item_with_guard(&arrd_coef, i, 2),
                model
            );
        }
        return;
    }
    else
    {
        double *grd = ((double*)arrd_grd.data);
        const double *x = ((double*)arrd_x.data);
        const double (*coeff)[2] = ((double(*)[2])arrd_coef.data);

        _model_eval(
            NULL, grd+0, grd+1, grd+2, x[0], x[1], x[2],
            coeff, model, SM_GRADIENT
        );
    }
}


static void _sheval_pot_and_grd(
    ARRAY_DATA arrd_pot,
    ARRAY_DATA arrd_grd,
    ARRAY_DATA arrd_x,
    ARRAY_DATA arrd_coef,
    MODEL *model
)
{
    if (arrd_pot.ndim > 0)
    {
        npy_intp i;
        for(i = 0; i < arrd_pot.dim[0]; ++i)
        {
            _sheval_pot_and_grd(
                _get_arrd_item_nocheck(&arrd_pot, i),
                _get_arrd_item_nocheck(&arrd_grd, i),
                _get_arrd_item_with_guard(&arrd_x, i, 1),
                _get_arrd_item_with_guard(&arrd_coef, i, 2),
                model
            );
        }
        return;
    }
    else
    {
        double *pot = ((double*)arrd_pot.data);
        double *grd = ((double*)arrd_grd.data);
        const double *x = ((double*)arrd_x.data);
        const double (*coeff)[2] = ((double(*)[2])arrd_coef.data);

        _model_eval(
            pot, grd+0, grd+1, grd+2, x[0], x[1], x[2],
            coeff, model, SM_POTENTIAL_AND_GRADIENT
        );
    }

}

/* single point evaluation */

static void _model_eval(
    double *fpot, double *fx, double *fy, double *fz,
    const double x, const double y, const double z,
    const double (*coeff)[2], MODEL *model, const int mode
)
{
    double glat, glon, ghgt;
    double clat, clon, crad;
    double flat, flon, frad;
    double tmp;

    void (*shc_relradpow)(double *rrp, int degree, double relrad) = (
        model->is_internal ? shc_relradpow_internal : shc_relradpow_external
    );

    // convert the input coordinates
    switch(model->coord_in)
    {
        case CT_GEODETIC_ABOVE_WGS84:
            glat = x; glon = y; ghgt = z;
            geodetic2geocentric_sph(&crad, &clat, &clon, glat, glon, ghgt, model->elps_a, model->elps_eps2);
            break;
        case CT_GEOCENTRIC_SPHERICAL:
            clat = DG2RAD*x; clon = DG2RAD*y; crad = z;
            if (model->coord_out == CT_GEODETIC_ABOVE_WGS84) {
                geocentric_sph2geodetic(&glat, &glon, &ghgt, crad, clat, clon, model->elps_a, model->elps_eps2);
            }
            break;
        case CT_GEOCENTRIC_CARTESIAN:
            cart2sph(&crad, &clat, &clon, x, y, z);
            if (model->coord_out == CT_GEODETIC_ABOVE_WGS84) {
                geocentric_sph2geodetic(&glat, &glon, &ghgt, crad, clat, clon, model->elps_a, model->elps_eps2);
            }
            break;
        default:
            return;
    }

    // associative Legendre functions
    if (model->clat_last != clat)
        shc_legendre(model->lp, model->ldp, model->degree, clat, model->psqrt, model->prsqrt);

    // longitude sines/cosines series
    if (model->clon_last != clon)
        fs_cos_sin(model->lcs, model->degree, clon);

    // relative radial powers
    if (model->crad_last != crad)
        shc_relradpow(model->rrp, model->degree, crad/model->crad_ref);

    // save the last evaluated coordinate
    model->clat_last = clat;
    model->clon_last = clon;
    model->crad_last = crad;

    // evaluate the model
    shc_eval(
        fpot, &flat, &flon, &frad, model->degree, mode,
        clat, crad, coeff, model->lcs, model->lp, model->ldp,
        model->rrp, model->is_internal
    );

    if (mode & SM_GRADIENT)
    {
        // project the produced vectors to the desired coordinate system
        switch(model->coord_out)
        {
            case CT_GEOCENTRIC_SPHERICAL:
                *fx = flat;
                *fy = flon;
                *fz = frad;
                break;

            case CT_GEODETIC_ABOVE_WGS84:
                tmp = clat - DG2RAD*glat;
                rot2d(fz, fx, frad, flat, sin(tmp), cos(tmp));
                *fy = flon;
                break;

            case CT_GEOCENTRIC_CARTESIAN:
                rot2d(&tmp, fz, frad, flat, sin(clat), cos(clat));
                rot2d(fx, fy, tmp, flon, sin(clon), cos(clon));
                break;
        }

        // final scaling
        *fx *= model->scale_gradient0;
        *fy *= model->scale_gradient1;
        *fz *= model->scale_gradient2;
    }

    if (mode & SM_POTENTIAL)
    {
        // final scaling
        *fpot *= model->scale_potential;
    }
}


/* model structure - reset */

static void _model_reset(MODEL *model)
{
    memset(model, 0, sizeof(MODEL));
}


/* model structure - destruction */

static void _model_destroy(MODEL *model)
{
    if(model->lp) free(model->lp);
    if(model->ldp) free(model->ldp);
    if(model->lcs) free(model->lcs);
    if(model->rrp) free(model->rrp);
    if(model->psqrt) free(model->psqrt);
    if(model->prsqrt) free(model->prsqrt);

    _model_reset(model);
}


/* model structure - initialization */

static int _model_init(
    MODEL *model, int is_internal, int degree, int coord_in, int coord_out,
    const double scale_potential, const double *scale_gradient)
{
    const int nterm = ((degree+1)*(degree+2))/2;

    _model_reset(model);

    model->is_internal = is_internal;
    model->degree = degree;
    model->nterm = nterm;
    model->coord_in = coord_in;
    model->coord_out = coord_out;
    model->elps_a = WGS84_A;
    model->elps_eps2 = WGS84_EPS2;
    model->crad_ref = RADIUS;
    model->clat_last = NAN;
    model->clon_last = NAN;
    model->crad_last = NAN;
    model->scale_potential = scale_potential;
    model->scale_gradient0 = scale_gradient[0];
    model->scale_gradient1 = scale_gradient[1];
    model->scale_gradient2 = scale_gradient[2];

    if (NULL == (model->lp = (double*)malloc(nterm*sizeof(double))))
    {
        PyErr_Format(PyExc_MemoryError, "_model_ts_init: lp");
        goto error;
    }

    if (NULL == (model->ldp = (double*)malloc(nterm*sizeof(double))))
    {
        PyErr_Format(PyExc_MemoryError, "_model_ts_init: ldp");
        goto error;
    }

    if (NULL == (model->lcs = (double(*)[2])malloc((degree+1)*sizeof(double[2]))))
    {
        PyErr_Format(PyExc_MemoryError, "_model_ts_init: lcs");
        goto error;
    }

    if (NULL == (model->rrp = (double*)malloc((degree+1)*sizeof(double))))
    {
        PyErr_Format(PyExc_MemoryError, "_model_ts_init: rrp");
        goto error;
    }

    if (NULL == (model->psqrt = shc_presqrt(degree)))
    {
        PyErr_Format(PyExc_MemoryError, "_model_ts_init: psqrt");
        goto error;
    }

    if (NULL == (model->prsqrt = shc_prersqrt(degree)))
    {
        PyErr_Format(PyExc_MemoryError, "_model_ts_init: prsqrt");
        goto error;
    }

    return 0;

  error:

    _model_destroy(model);

    return 1;
}

#endif  /* PYMM_SHEVAL_H */
