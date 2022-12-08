/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - spherical harmonic model evaluation with time-dependent coefficients
 *   evaluated by means of the 2D Fourier series (MIO)
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *-----------------------------------------------------------------------------
 * Copyright (C) 2022 EOX IT Services GmbH
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

#ifndef PYMM_SHEVAL2DFS_H
#define PYMM_SHEVAL2DFS_H

#include "fourier_series.h"
#include "pymm_aux.h"
#include "pymm_sheval.h"
#include "pymm_fourier2d.h"

/* coefficient set - auxiliary structure */
typedef struct Fourier2DCoefSet {
    PyArrayObject *arr_ab;
    PyArrayObject *arr_cm;
    const double (*ab)[2];
    const int *cm;
    ptrdiff_t *offset;
    size_t ncoef;
    int degree;
    FOURIER_2D f2d;
} FOURIER_2D_COEF_SET;


static int _f2d_coefset_parse(
    FOURIER_2D_COEF_SET *coefset,
    const char *label,
    PyObject *obj_coef_set
);
static int _f2d_coefset_init(
    FOURIER_2D_COEF_SET *coefset,
    PyArrayObject *arr_ab,
    PyArrayObject *arr_cm,
    const ptrdiff_t min_degree1,
    const ptrdiff_t min_degree2,
    const double scale1,
    const double scale2
);
static void _f2d_coefset_destroy(FOURIER_2D_COEF_SET *coefset);
static void _f2d_coefset_reset(FOURIER_2D_COEF_SET *coefset);

/* 2D Fourier series magnetic model - auxiliary structure */
typedef struct Model2DFS {
    MODEL sh_model;
    double (*coef)[2];
    FOURIER_2D_COEF_SET *coefset;
    double time1_last;
    double time2_last;
} MODEL_2DFS;

static void _model_2dfs_reset(MODEL_2DFS *model);
static void _model_2dfs_destroy(MODEL_2DFS *model);
static int _model_2dfs_init(
    MODEL_2DFS *model, const int is_internal, const int degree,
    const int coord_in, const int coord_out,
    FOURIER_2D_COEF_SET *coefset,
    const double scale_potential, const double *scale_gradient
);

static void _sheval2dfs_pot(
    ARRAY_DATA arrd_pot,
    ARRAY_DATA arrd_t1,
    ARRAY_DATA arrd_t2,
    ARRAY_DATA arrd_x,
    MODEL_2DFS *model
);
static void _sheval2dfs_grd(
    ARRAY_DATA arrd_grd,
    ARRAY_DATA arrd_t1,
    ARRAY_DATA arrd_t2,
    ARRAY_DATA arrd_x,
    MODEL_2DFS *model
);
static void _sheval2dfs_pot_and_grd(
    ARRAY_DATA arrd_pot,
    ARRAY_DATA arrd_grd,
    ARRAY_DATA arrd_t1,
    ARRAY_DATA arrd_t2,
    ARRAY_DATA arrd_x,
    MODEL_2DFS *model
);


/* Python function definition */

#define DOC_SHEVAL2DFS "\n"\
"   arr_out = sheval2dfs(\n"\
"       arr_t1, arr_t2, arr_x, coef_set,\n"\
"       coord_type_in=GEODETIC_ABOVE_WGS84,\n"\
"       coord_type_out=GEODETIC_ABOVE_WGS84,\n"\
"       mode=GRADIENT, is_internal=True,\n"\
"       scale_potential=1.0, scale_gradient=1.0)\n"\
"   )\n"\
"\n"\
"     Parameters:\n"\
"       arr_t1 - array of seasonal times.\n"\
"       arr_t2 - array of diurnal times.\n"\
"       arr_x  - array of 3D coordinates.\n"\
"       coef_set = (coef_ab, coef_nm, min_degree1, min_degree2, scale1)\n"\
"              - a sets of 2D Fourier series coefficients. \n"\
"       coef_ab[c,n1,n2,2] - series of the 2D Fourier series coefficients\n"\
"       coef_nm - [c,2] mapping of the spherical harmonic model coefficients.\n"\
"       min_degree1 - seasonal time minimal Fourier series min. degree\n"\
"               (max. degree is determined from the coef_abdimension n1)\n"\
"       min_degree2 - diurnal time minimal Fourier series min. degree\n"\
"               (max. degree is determined from the coef_ab dimension n2)\n"\
"       scale1 - seasonal time scaling factor\n"\
"       scale2 - diurnal time scaling factor\n"\
"       coord_type_in - type of the input coordinates.\n"\
"       coord_type_our - type of the output coordinates frame.\n"\
"       mode - quantity to be evaluated:\n"\
"                  POTENTIAL\n"\
"                  GRADIENT (default)\n"\
"                  POTENTIAL_AND_GRADIENT\n"\
"       is_internal - boolean flag set to True by default. When set to False \n"\
"                     external field evaluation is used.\n"\
"       scale_potential - scalar value multiplied with the result potentials.\n"\
"       scale_gradient - scalar or 3 element array multiplied with the result\n"\
"                        gradient components.\n"\
"\n"


static PyObject* sheval2dfs(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "arr_t1", "arr_t2", "arr_x", "coef_set",
        "coord_type_in", "coord_type_out", "mode",
        "is_internal", "scale_potential", "scale_gradient", NULL
    };

    int degree = 0, mode = 0x2, is_internal;
    int ct_in = CT_GEODETIC_ABOVE_WGS84;
    int ct_out = CT_GEODETIC_ABOVE_WGS84;
    double scale_potential = 1.0;
    double scale_gradient[3] = {1.0, 1.0, 1.0};
    PyObject *obj_is_internal = NULL; // boolean flag
    PyObject *obj_t1 = NULL; // times
    PyObject *obj_t2 = NULL; // times
    PyObject *obj_x = NULL; // spatial coordinates
    PyObject *obj_coef_set = NULL; // coefficient set
    PyObject *obj_scale = NULL; // gradient scale object
    PyObject *retval = NULL;
    PyArrayObject *arr_t1 = NULL; // times array
    PyArrayObject *arr_t2 = NULL; // times array
    PyArrayObject *arr_x = NULL; // spatial coordinates array
    PyArrayObject *arr_pot = NULL; // output array
    PyArrayObject *arr_grd = NULL; // output array
    PyArrayObject *arr_scale = NULL; // gradient scale array
    FOURIER_2D_COEF_SET coefset = {0};
    MODEL_2DFS model = {0};

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOOO|iiiOdO:sheval2dfs", keywords,
        &obj_t1, &obj_t2, &obj_x, &obj_coef_set, &ct_in, &ct_out, &mode,
        &obj_is_internal, &scale_potential, &obj_scale
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
    if (NULL == (arr_t1 = _get_as_double_array(obj_t1, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_t2 = _get_as_double_array(obj_t2, 0, 0, NPY_ARRAY_ALIGNED, keywords[1])))
        goto exit;

    if (NULL == (arr_x = _get_as_double_array(obj_x, 1, 0, NPY_ARRAY_ALIGNED, keywords[2])))
        goto exit;

    // check the arrays' dimensions
    {
        npy_intp i, ndim = MIN(PyArray_NDIM(arr_t1), PyArray_NDIM(arr_t2));
        for (i = 0; i < ndim; ++i) {
            if (PyArray_DIMS(arr_t1)[i] != PyArray_DIMS(arr_t2)[i]) {
                PyErr_Format(PyExc_ValueError, "Shape of %s array "\
                    " does not match shape of %s!", keywords[0], keywords[1]);
                goto exit;
            }
        }
    }
    {
        npy_intp i, ndim = MIN(PyArray_NDIM(arr_t1), PyArray_NDIM(arr_x) - 1);
        for (i = 0; i < ndim; ++i) {
            if (PyArray_DIMS(arr_t1)[i] != PyArray_DIMS(arr_x)[i]) {
                PyErr_Format(PyExc_ValueError, "Shape of %s array "\
                    " does not match shape of %s!", keywords[0], keywords[2]);
                goto exit;
            }
        }
    }
    {
        npy_intp i, ndim = MIN(PyArray_NDIM(arr_t2), PyArray_NDIM(arr_x) - 1);
        for (i = 0; i < ndim; ++i) {
            if (PyArray_DIMS(arr_t2)[i] != PyArray_DIMS(arr_x)[i]) {
                PyErr_Format(PyExc_ValueError, "Shape of %s array "\
                    " does not match shape of %s!", keywords[1], keywords[2]);
                goto exit;
            }
        }
    }

    // check the last dimension of the position array
    if (_check_array_dim_eq(arr_x, -1, 3, keywords[0]))
        goto exit;


    // parse coefficient set
    _f2d_coefset_parse(&coefset, keywords[3], obj_coef_set);

    degree = coefset.degree;

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "Invalid model degree %d!", degree);
        goto exit;
    }

    // handle the gradient scale factors
    if (NULL != obj_scale)
    {
        if (NULL == (arr_scale = _get_as_double_array(obj_scale, 0, 1, NPY_ARRAY_IN_ARRAY, keywords[9])))
            goto exit;

        if (_extract_1d_double_array(scale_gradient, 3, arr_scale, keywords[9]))
            goto exit;
    }

    // create new output arrays
    {
        npy_intp ndim = 0;
        npy_intp *dims_src = NULL;

        ndim = MAX(ndim, PyArray_NDIM(arr_t1));
        ndim = MAX(ndim, PyArray_NDIM(arr_t2));
        ndim = MAX(ndim, PyArray_NDIM(arr_x) - 1);

        if (PyArray_NDIM(arr_t1) == ndim)
            dims_src = PyArray_DIMS(arr_t1);
        else if (PyArray_NDIM(arr_t2) == ndim)
            dims_src = PyArray_DIMS(arr_t2);
        else
            dims_src = PyArray_DIMS(arr_x);

        {
            npy_intp i;
            npy_intp dims[ndim+1];

            for (i = 0; i < ndim; ++i)
                dims[i] = dims_src[i];
            dims[ndim] = 3;

            if (mode & SM_POTENTIAL)
                if (NULL == (arr_pot = _get_new_array(ndim, dims, NPY_DOUBLE)))
                    goto exit;

            if (mode & SM_GRADIENT)
                if (NULL == (arr_grd = _get_new_array(ndim + 1, dims, NPY_DOUBLE)))
                    goto exit;
        }
    }

    // evaluate model

    if (_model_2dfs_init(
        &model, is_internal, degree, ct_in, ct_out,
        &coefset, scale_potential, scale_gradient
    ))
        goto exit;


    switch(mode)
    {
        case SM_POTENTIAL:
            _sheval2dfs_pot(
                _array_to_arrd(arr_pot),
                _array_to_arrd(arr_t1),
                _array_to_arrd(arr_t2),
                _array_to_arrd(arr_x),
                &model
            );
            retval = (PyObject*) arr_pot;
            break;

        case SM_GRADIENT:
            _sheval2dfs_grd(
                _array_to_arrd(arr_grd),
                _array_to_arrd(arr_t1),
                _array_to_arrd(arr_t2),
                _array_to_arrd(arr_x),
                 &model
            );
            retval = (PyObject*) arr_grd;
            break;

        case SM_POTENTIAL_AND_GRADIENT:
            _sheval2dfs_pot_and_grd(
                _array_to_arrd(arr_pot),
                _array_to_arrd(arr_grd),
                _array_to_arrd(arr_t1),
                _array_to_arrd(arr_t2),
                _array_to_arrd(arr_x),
                 &model
            );
            if (NULL == (retval = Py_BuildValue("NN", (PyObject*) arr_pot, (PyObject*) arr_grd)))
                goto exit;
            break;
    }

  exit:

    _model_2dfs_destroy(&model);
    _f2d_coefset_destroy(&coefset);

    // decrease reference counters to the arrays
    if (arr_t1) Py_DECREF(arr_t1);
    if (arr_t2) Py_DECREF(arr_t2);
    if (arr_x) Py_DECREF(arr_x);
    if (arr_scale) Py_DECREF(arr_scale);
    if (!retval && arr_grd) Py_DECREF(arr_grd);
    if (!retval && arr_pot) Py_DECREF(arr_pot);

    return retval;
}


/* high-level nD-array recursive batch model_evaluation */

static void _sheval2dfs_eval_coeff(
    const double time1,
    const double time2,
    MODEL_2DFS *model
);

static void _sheval2dfs_pot(
    ARRAY_DATA arrd_pot,
    ARRAY_DATA arrd_t1,
    ARRAY_DATA arrd_t2,
    ARRAY_DATA arrd_x,
    MODEL_2DFS *model
)
{
    if ((arrd_t1.ndim > 0) || (arrd_t2.ndim > 0))
    {
        npy_intp i, n = arrd_t1.ndim > 0 ? arrd_t1.dim[0] : arrd_t2.dim[0];
        for(i = 0; i < n; ++i)
        {
            _sheval2dfs_pot(
                _get_arrd_item(&arrd_pot, i),
                _get_arrd_item(&arrd_t1, i),
                _get_arrd_item(&arrd_t2, i),
                _get_arrd_vector_item(&arrd_x, i),
                model
            );
        }
        return;
    }
    else
    {
        const double time1 = *((double*)arrd_t1.data);
        const double time2 = *((double*)arrd_t2.data);
        ARRAY_DATA arrd_coef = {
            .data=model->coef, .ndim=0, .dim=NULL, .stride=NULL
        };

        if ((time1 != model->time1_last)||(time2 != model->time2_last))
            _sheval2dfs_eval_coeff(time1, time2, model);

        _sheval_pot(arrd_pot, arrd_x, arrd_coef, &(model->sh_model));
    }
}

static void _sheval2dfs_grd(
    ARRAY_DATA arrd_grd,
    ARRAY_DATA arrd_t1,
    ARRAY_DATA arrd_t2,
    ARRAY_DATA arrd_x,
    MODEL_2DFS *model
)
{
    if ((arrd_t1.ndim > 0) || (arrd_t2.ndim > 0))
    {
        npy_intp i, n = arrd_t1.ndim > 0 ? arrd_t1.dim[0] : arrd_t2.dim[0];
        for(i = 0; i < n; ++i)
        {
            _sheval2dfs_grd(
                _get_arrd_item(&arrd_grd, i),
                _get_arrd_item(&arrd_t1, i),
                _get_arrd_item(&arrd_t2, i),
                _get_arrd_vector_item(&arrd_x, i),
                model
            );
        }
        return;
    }
    else
    {
        const double time1 = *((double*)arrd_t1.data);
        const double time2 = *((double*)arrd_t2.data);
        ARRAY_DATA arrd_coef = {
            .data=model->coef, .ndim=0, .dim=NULL, .stride=NULL
        };

        if ((time1 != model->time1_last)||(time2 != model->time2_last))
            _sheval2dfs_eval_coeff(time1, time2, model);

        _sheval_grd(arrd_grd, arrd_x, arrd_coef, &(model->sh_model));
    }
}

static void _sheval2dfs_pot_and_grd(
    ARRAY_DATA arrd_pot,
    ARRAY_DATA arrd_grd,
    ARRAY_DATA arrd_t1,
    ARRAY_DATA arrd_t2,
    ARRAY_DATA arrd_x,
    MODEL_2DFS *model
)
{
    if ((arrd_t1.ndim > 0) || (arrd_t2.ndim > 0))
    {
        npy_intp i, n = arrd_t1.ndim > 0 ? arrd_t1.dim[0] : arrd_t2.dim[0];
        for(i = 0; i < n; ++i)
        {
            _sheval2dfs_pot_and_grd(
                _get_arrd_item(&arrd_pot, i),
                _get_arrd_item(&arrd_grd, i),
                _get_arrd_item(&arrd_t1, i),
                _get_arrd_item(&arrd_t2, i),
                _get_arrd_vector_item(&arrd_x, i),
                model
            );
        }
        return;
    }
    else
    {
        const double time1 = *((double*)arrd_t1.data);
        const double time2 = *((double*)arrd_t2.data);
        ARRAY_DATA arrd_coef = {
            .data=model->coef, .ndim=0, .dim=NULL, .stride=NULL
        };

        if ((time1 != model->time1_last)||(time2 != model->time2_last))
            _sheval2dfs_eval_coeff(time1, time2, model);

        _sheval_pot_and_grd(arrd_pot, arrd_grd, arrd_x, arrd_coef, &(model->sh_model));
    }
}

static void _sheval2dfs_eval_coeff(
    const double time1,
    const double time2,
    MODEL_2DFS *model
)
{
    FOURIER_2D_COEF_SET *coefset = model->coefset;
    FOURIER_2D *f2d = &coefset->f2d;

    double (*coef)[2] = model->coef;

    const ptrdiff_t size1 = f2d->max_degree1 - f2d->min_degree1 + 1;
    const ptrdiff_t size2 = f2d->max_degree2 - f2d->min_degree2 + 1;
    const ptrdiff_t ab_stride = size1 * size2;

    const ptrdiff_t n = coefset->ncoef;
    ptrdiff_t i;

    _fourier2d_cossin(f2d, time1, time2);

    for (i = 0; i < n; ++i)
    {
        size_t offset = abs(coefset->offset[i]);
        size_t kind = coefset->offset[i] < 0;
        const double (*ab)[2] = coefset->ab + i*ab_stride ;

        coef[offset][kind] =_fourier2d_eval(f2d, ab);
    }

    model->time1_last = time1;
    model->time2_last = time2;
}

/* model structure - zero reset */

static void _model_2dfs_reset(MODEL_2DFS *model) {
    memset(model, 0, sizeof(MODEL_2DFS));
}


/* model structure - destruction */

static void _model_2dfs_destroy(MODEL_2DFS *model)
{
    if (model->coef) free((double*)model->coef);

    _model_destroy(&(model->sh_model));

    _model_2dfs_reset(model);
}


/* model structure - initialization */

static int _model_2dfs_init(
    MODEL_2DFS *model, const int is_internal, const int degree,
    const int coord_in, const int coord_out,
    FOURIER_2D_COEF_SET *coefset,
    const double scale_potential, const double *scale_gradient
)
{
    _model_2dfs_reset(model);

    // initialize nested single-time model
    if (_model_init(
        &(model->sh_model), is_internal, degree, coord_in, coord_out,
        scale_potential, scale_gradient
    ))
        goto error;

    // allocate memory for the single time coefficients
    if (NULL == (model->coef = (double(*)[2])calloc(model->sh_model.nterm, sizeof(double[2]))))
    {
        PyErr_Format(PyExc_MemoryError, "_model_2dfs_init: coef");
        goto error;
    }

    // initialize coefficients time-series
    model->coefset = coefset;
    model->time1_last = NAN;
    model->time2_last = NAN;

    return 0;

  error:

    _model_2dfs_destroy(model);

    return 1;
}


/* 2D Fourier series coefficients set - reset */

static void _f2d_coefset_reset(FOURIER_2D_COEF_SET *coefset)
{
    memset(coefset, 0, sizeof(FOURIER_2D_COEF_SET));
}


/* 2D Fourier series coefficients set - destroy */

static void _f2d_coefset_destroy(FOURIER_2D_COEF_SET *coefset)
{
    if (coefset->arr_ab) Py_DECREF(coefset->arr_ab);
    if (coefset->arr_cm) Py_DECREF(coefset->arr_cm);

    _fourier2d_destroy(&coefset->f2d);

    _f2d_coefset_reset(coefset);
}


/* 2D Fourier series coefficients set - initialize */

static int _f2d_coefset_init(
    FOURIER_2D_COEF_SET *coefset,
    PyArrayObject *arr_ab,
    PyArrayObject *arr_cm,
    const ptrdiff_t min_degree1,
    const ptrdiff_t min_degree2,
    const double scale1,
    const double scale2
)
{
    const ptrdiff_t size1 = PyArray_DIMS(arr_ab)[1];
    const ptrdiff_t size2 = PyArray_DIMS(arr_ab)[2];

    coefset->arr_ab = arr_ab;
    coefset->arr_cm = arr_cm;

    coefset->ab = (double(*)[2])PyArray_DATA(arr_ab),
    coefset->cm = (int*) PyArray_DATA(arr_cm),
    coefset->ncoef = PyArray_DIMS(arr_cm)[0];

    if (NULL == (coefset->offset = (ptrdiff_t*)malloc(coefset->ncoef*sizeof(ptrdiff_t))))
    {
        PyErr_Format(PyExc_MemoryError, "_f2d_coefset_init: offset");
        goto error;
    }

    {
        size_t i;

        for (i = 0; i < coefset->ncoef; ++i) {
            int degree = coefset->cm[i*2];
            int order = coefset->cm[i*2+1];
            ptrdiff_t offset = ((ptrdiff_t)degree)*((ptrdiff_t)degree+1)/2 + ABS(order);
            ptrdiff_t sign = order < 0 ? -1 : 1;

            coefset->offset[i] = sign * offset;
            coefset->degree = MAX(coefset->degree, degree);
        }
    }

    if (_fourier2d_init(&coefset->f2d, min_degree1, size1, min_degree2, size2, scale1, scale2))
        goto error;

    return 0;

  error:

    _f2d_coefset_destroy(coefset);

    return 1;
}


/* 2D Fourier series coefficients set - parse inputs */

static int _f2d_coefset_parse(
    FOURIER_2D_COEF_SET *coefset,
    const char *label,
    PyObject *obj_coef_set
)
{
    static char *keywords[] = {
        "coef_ab", "coef_cm", "min_degree1", "min_degree2", "scale1", "scale2", NULL
    };

    PyObject *obj_ab = NULL;
    PyObject *obj_cm = NULL;
    PyArrayObject *arr_ab = NULL;
    PyArrayObject *arr_cm = NULL;

    ptrdiff_t min_degree1;
    ptrdiff_t min_degree2;
    double scale1;
    double scale2;

    const size_t buffer_size = 128;
    char buffer[buffer_size];

    // check tuple item
    if (!PyTuple_Check(obj_coef_set))
    {
        PyErr_Format(PyExc_ValueError, "The '%s' is expected to be a tuple!", label);
        goto error;
    }

    if (PyTuple_Size(obj_coef_set) != 6)
    {
        PyErr_Format(PyExc_ValueError, "The '%s' tuple is expected to have 6 items!", label);
        goto error;
    }

    // extract and parse set values

    if (NULL == (obj_ab = PyTuple_GetItem(obj_coef_set, 0)))
        goto error;

    if (NULL == (obj_cm = PyTuple_GetItem(obj_coef_set, 1)))
        goto error;

    snprintf(buffer, buffer_size, "%s.%s", label, keywords[2]);
    if (_parse_long_value(&min_degree1, PyTuple_GetItem(obj_coef_set, 2), buffer))
        goto error;

    snprintf(buffer, buffer_size, "%s.%s", label, keywords[3]);
    if (_parse_long_value(&min_degree2, PyTuple_GetItem(obj_coef_set, 3), buffer))
        goto error;

    snprintf(buffer, buffer_size, "%s.%s", label, keywords[4]);
    if (_parse_double_value(&scale1, PyTuple_GetItem(obj_coef_set, 4), buffer))
        goto error;

    snprintf(buffer, buffer_size, "%s.%s", label, keywords[5]);
    if (_parse_double_value(&scale2, PyTuple_GetItem(obj_coef_set, 5), buffer))
        goto error;

    // extract array objects

    snprintf(buffer, buffer_size, "%s.%s", label, keywords[0]);
    if (NULL == (arr_ab = _get_as_double_array(obj_ab, 4, 4, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_IN_ARRAY, buffer)))
        goto error;

    snprintf(buffer, buffer_size, "%s.%s", label, keywords[1]);
    if (NULL == (arr_cm = _get_as_int_array(obj_cm, 2, 2, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_IN_ARRAY, buffer)))
        goto error;

    // check dimensions

    if (PyArray_DIMS(arr_ab)[3] != 2) {
        PyErr_Format(PyExc_ValueError, "Incorrect shape of %s.%s!", label, keywords[1]);
        goto error;
    }

    if (PyArray_DIMS(arr_cm)[1] != 2) {
        PyErr_Format(PyExc_ValueError, "Incorrect shape of %s.%s!", label, keywords[1]);
        goto error;
    }

    if (PyArray_DIMS(arr_ab)[0] != PyArray_DIMS(arr_cm)[0]) {
        PyErr_Format(PyExc_ValueError, "Shape mismatch between %s.%s and %s.%s!", label, keywords[0], label, keywords[1]);
        goto error;
    }

    return _f2d_coefset_init(
        coefset, arr_ab, arr_cm,
        min_degree1, min_degree2, scale1, scale2
    );

  error:

    if (arr_ab) Py_DECREF(arr_ab);
    if (arr_cm) Py_DECREF(arr_cm);

    return 1;
}

#endif  /* PYMM_SHEVAL2DFS_H */
