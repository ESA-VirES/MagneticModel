/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - spherical harmonic model evaluation with temporal coefficient dependency
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

#ifndef PYMM_SHEVALTEMP_H
#define PYMM_SHEVALTEMP_H

#include "interpolation.h"
#include "pymm_sheval.h"

#define SHEVALTEMP_MIN_SPLINE_ORDER 1
#define SHEVALTEMP_MAX_SPLINE_ORDER 2
#define SHEVALTEMP_DEFAULT_SPLINE_ORDER 2

#define SHEVALTEMP_MAX_COEF_LIST_SIZE 16

#ifndef NAN
#define NAN (0.0/0.0)
#endif

#ifndef MIN
#define MIN(a,b) ((b)<(a)?(b):(a))
#endif

#ifndef MAX
#define MAX(a,b) ((b)>(a)?(b):(a))
#endif

#ifndef ABS
#define ABS(a) ((a)<0?-(a):(a))
#endif

typedef struct ModelTs MODEL_TS;
typedef struct CoefSet COEF_SET;

typedef void (*f_coeff_set_interp)(double time, MODEL_TS *model, const COEF_SET *coefset);
static void _coeff_set_interp0(double time, MODEL_TS *model, const COEF_SET *coefset);
static void _coeff_set_interp1(double time, MODEL_TS *model, const COEF_SET *coefset);

/* deficient set - auxiliary structure */
typedef struct CoefSet {
    PyArrayObject *arr_ct;
    PyArrayObject *arr_cv;
    PyArrayObject *arr_cm;
    const double *ct;
    const double *cv;
    const int *cm;
    ptrdiff_t *offset;
    f_coeff_set_interp coeff_set_interp;
    size_t ntime;
    size_t ncoef;
    int spline_order;
    int degree;
    ptrdiff_t idx_last;
} COEF_SET;


/* coefficient set - auxiliary structure */
static void _coefset_reset(COEF_SET *coefset) {
    memset(coefset, 0, sizeof(COEF_SET));
}

/* model structure - destruction */
static void _coefset_destroy(COEF_SET *coefset)
{
    if(NULL != coefset->offset)
        free(coefset->offset);

    if (coefset->arr_ct) Py_DECREF(coefset->arr_ct);
    if (coefset->arr_cv) Py_DECREF(coefset->arr_cv);
    if (coefset->arr_cm) Py_DECREF(coefset->arr_cm);

    _coefset_reset(coefset);
}

/* model structure - initialization */

static int _coefset_init(
    COEF_SET *coefset,
    PyArrayObject *arr_ct,
    PyArrayObject *arr_cv,
    PyArrayObject *arr_cm,
    int spline_order)
{
    _coefset_reset(coefset);

    coefset->arr_ct = arr_ct;
    coefset->arr_cv = arr_cv;
    coefset->arr_cm = arr_cm;
    coefset->ct = PyArray_DATA(arr_ct),
    coefset->cv = PyArray_DATA(arr_cv),
    coefset->cm = PyArray_DATA(arr_cm),
    coefset->ntime = PyArray_DIMS(arr_ct)[0];
    coefset->ncoef = PyArray_DIMS(arr_cm)[0];
    coefset->spline_order = spline_order;
    coefset->degree = 0;
    coefset->idx_last = -1;

    if (NULL == (coefset->offset = (ptrdiff_t*)calloc(coefset->ncoef, sizeof(ptrdiff_t))))
    {
        PyErr_Format(PyExc_MemoryError, "_coefset_init: coefset");
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

    switch (spline_order) {
        case 1:
            coefset->coeff_set_interp = _coeff_set_interp0;
            break;
        case 2:
            coefset->coeff_set_interp = _coeff_set_interp1;
            break;
        default:
            PyErr_Format(PyExc_ValueError, "Invalid spline order %d!", spline_order);
            goto error;
    }

    return 0;

  error:
    _coefset_destroy(coefset);
    return 1;
}

/* magnetic model - auxiliary structure */
typedef struct ModelTs {
    MODEL sh_model;
    const COEF_SET *coefsets;
    double time_last;
    size_t ncoefset;
} MODEL_TS;

/* model structure - zero reset */
static void _shevaltemp_model_reset(MODEL_TS *model) {
    memset(model, 0, sizeof(MODEL_TS));
}

/* model structure - destruction */
static void _shevaltemp_model_destroy(MODEL_TS *model)
{
    if(NULL != model->sh_model.cg)
        free((double*)model->sh_model.cg);

    if(NULL != model->sh_model.ch)
        free((double*)model->sh_model.ch);

    _sheval_model_destroy(&(model->sh_model));

    _shevaltemp_model_reset(model);
}

/* model structure - initialization */
static int _shevaltemp_model_init(
    MODEL_TS *model, const int is_internal, const int degree,
    const int coord_in, const int coord_out,
    const COEF_SET *coefsets, size_t ncoefset,
    const double scale_potential, const double *scale_gradient)
{
    _shevaltemp_model_reset(model);

    // initialize nested single-time model
    if (_sheval_model_init(
        &(model->sh_model), is_internal, degree, coord_in, coord_out,
        NULL, NULL, scale_potential, scale_gradient
    ))
        goto error;

    // allocate memory for the single time coefficients
    if (NULL == (model->sh_model.cg = (double*)calloc(model->sh_model.nterm, sizeof(double))))
    {
        PyErr_Format(PyExc_MemoryError, "_shevaltemp_model_init: cg");
        goto error;
    }

    if (NULL == (model->sh_model.ch = (double*)calloc(model->sh_model.nterm, sizeof(double))))
    {
        PyErr_Format(PyExc_MemoryError, "_shevaltemp_model_init: ch");
        goto error;
    }

    // initialize coefficients time-series
    model->coefsets = coefsets;
    model->ncoefset = ncoefset;
    model->time_last = NAN;

    return 0;

  error:
    _shevaltemp_model_destroy(model);
    return 1;
}

/* coefficients time-series interpolation */

static void _coeff_set_interp0(double time, MODEL_TS *model, const COEF_SET *coefset) {
    INTERP_BASIS basis = get_interp0_basis(time, coefset->ct, coefset->ntime);

    if (coefset->idx_last == basis.i) {
        return;
    }

    double *cg = (double*)model->sh_model.cg;
    double *ch = (double*)model->sh_model.ch;

    size_t i, n = coefset->ncoef;

    for (i = 0; i < n; ++i)
    {
        double *target = (coefset->offset[i] >= 0 ? cg : ch);
        ptrdiff_t offset = abs(coefset->offset[i]);
        target[offset] = interp0_eval(&basis, coefset->cv + i*coefset->ntime);
    }

     ((COEF_SET*)coefset)->idx_last = basis.i;
}

static void _coeff_set_interp1(double time, MODEL_TS *model, const COEF_SET *coefset) {
    INTERP_BASIS basis = get_interp1_basis(time, coefset->ct, coefset->ntime);

    double *cg = (double*)model->sh_model.cg;
    double *ch = (double*)model->sh_model.ch;

    size_t i, n = coefset->ncoef;

    for (i = 0; i < n; ++i)
    {
        double *target = (coefset->offset[i] >= 0 ? cg : ch);
        ptrdiff_t offset = abs(coefset->offset[i]);
        target[offset] = interp1_eval(&basis, coefset->cv + i*coefset->ntime);
    }
}

static void _coeff_interp(double time, MODEL_TS *model) {

    size_t i, n = model->ncoefset;

    for (i = 0; i < n; ++i) {
        const COEF_SET *coefset = model->coefsets + i;
        coefset->coeff_set_interp(time, model, coefset);
    }

    model->time_last = time;
}



/* high-level nD-array recursive batch model_evaluation */

static void _shevaltemp1(ARRAY_DATA arrd_t, ARRAY_DATA arrd_x, ARRAY_DATA arrd_pot, MODEL_TS *model)
{
    if (arrd_t.ndim > 0)
    {
        npy_intp i, n = arrd_t.dim[0];
        for(i = 0; i < n; ++i)
        {
            _shevaltemp1(
                _get_arrd_item(&arrd_t, i),
                _get_arrd_vector_item(&arrd_x, i),
                _get_arrd_item(&arrd_pot, i),
                model
            );
        }
        return;
    }
    else
    {
        const double time = *((double*)arrd_t.data);
        if (time != model->time_last)
            _coeff_interp(time, model);
        _sheval1(arrd_x, arrd_pot, &(model->sh_model));
    }
}

static void _shevaltemp2(ARRAY_DATA arrd_t, ARRAY_DATA arrd_x, ARRAY_DATA arrd_grd, MODEL_TS *model)
{
    if (arrd_t.ndim > 0)
    {
        npy_intp i, n = arrd_t.dim[0];
        for(i = 0; i < n; ++i)
        {
            _shevaltemp2(
                _get_arrd_item(&arrd_t, i),
                _get_arrd_vector_item(&arrd_x, i),
                _get_arrd_item(&arrd_grd, i),
                model
            );
        }
        return;
    }
    else
    {
        const double time = *((double*)arrd_t.data);
        if (time != model->time_last)
            _coeff_interp(time, model);
        _sheval2(arrd_x, arrd_grd, &(model->sh_model));
    }
}


static void _shevaltemp3(ARRAY_DATA arrd_t, ARRAY_DATA arrd_x, ARRAY_DATA arrd_pot, ARRAY_DATA arrd_grd, MODEL_TS *model)
{
    if (arrd_t.ndim > 0)
    {
        npy_intp i, n = arrd_t.dim[0];
        for(i = 0; i < n; ++i)
        {
            _shevaltemp3(
                _get_arrd_item(&arrd_t, i),
                _get_arrd_vector_item(&arrd_x, i),
                _get_arrd_item(&arrd_pot, i),
                _get_arrd_item(&arrd_grd, i),
                model
            );
        }
        return;
    }
    else
    {
        const double time = *((double*)arrd_t.data);
        if (time != model->time_last)
            _coeff_interp(time, model);
        _sheval3(arrd_x, arrd_pot, arrd_grd, &(model->sh_model));
    }
}

/* Python function definition */

#define DOC_SHEVALTEMP "\n"\
"   arr_out = shevaltemp(\n"\
"       arr_t, arr_x, degree,\n"\
"       coef_set_list,\n"\
"       coord_type_in=GEODETIC_ABOVE_WGS84,\n"\
"       coord_type_out=GEODETIC_ABOVE_WGS84,\n"\
"       mode=GRADIENT, is_internal=True,\n"\
"       scale_potential=1.0, scale_gradient=1.0)\n"\
"   )\n"\
"\n"\
"     Parameters:\n"\
"       arr_t  - array of times.\n"\
"       arr_x  - array of 3D coordinates.\n"\
"       degree - degree of the spherical harmonic model.\n"\
"       coef_set_list = [(coef_t, coef_gh, coef_nm, spline_order), ...]\n"\
"              - a list of sets of interpolated cooeficients. \n"\
"       coef_t - times of spherical harmonic model coefficients.\n"\
"       coef_gh - time-series of spherical harmonic model coefficients.\n"\
"       coef_nm - N,M mapping of the spherical harmonic model coefficients.\n"\
"       spline_order - spline order of the temporal extrapolation.\n"\
"                      (currently, the only allowed value is 2)\n"\
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

static int parse_coefset(COEF_SET *coefset, size_t idx, PyObject *obj_coef_list_item);

static PyObject* shevaltemp(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "arr_t", "arr_x", "degree", "coef_set_list",
        "coord_type_in", "coord_type_out", "mode",
        "is_internal", "scale_potential", "scale_gradient", NULL
    };
    int ct_in = CT_GEODETIC_ABOVE_WGS84;
    int ct_out = CT_GEODETIC_ABOVE_WGS84;
    int degree = 0, mode = 0x2, is_internal;
    size_t ncoefset = 0;
    double scale_potential = 1.0;
    double scale_gradient[3] = {1.0, 1.0, 1.0};
    PyObject *obj_is_internal = NULL; // boolean flag
    PyObject *obj_t = NULL; // times
    PyObject *obj_x = NULL; // spatial coordinates
    PyObject *obj_coef_list = NULL; // list of coefficient objects
    PyObject *obj_scale = NULL; // gradient scale object
    PyArrayObject *arr_t = NULL; // times array
    PyArrayObject *arr_x = NULL; // spatial coordinates array
    PyArrayObject *arr_pot = NULL; // output array
    PyArrayObject *arr_grd = NULL; // output array
    PyArrayObject *arr_scale = NULL; // gradient scale array
    PyObject *retval = NULL;
    MODEL_TS model;
    COEF_SET coefsets[SHEVALTEMP_MAX_COEF_LIST_SIZE];

    // clear model structure
    _shevaltemp_model_reset(&model);

    // clear coefficient sets
    {
        size_t i;
        for (i = 0; i < SHEVALTEMP_MAX_COEF_LIST_SIZE; ++i)
            _coefset_reset(&coefsets[i]);
    }

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOiO|iiiOdO:shevaltemp", keywords,
        &obj_t, &obj_x, &degree,
        &obj_coef_list,
        &ct_in, &ct_out, &mode,
        &obj_is_internal, &scale_potential, &obj_scale
    ))
        goto exit;

    is_internal = (obj_is_internal == NULL) || PyObject_IsTrue(obj_is_internal);

    // check the type of the coordinate transformation
    if (CT_INVALID == _check_coord_type(ct_in, keywords[7]))
        goto exit;

    if (CT_INVALID == _check_coord_type(ct_out, keywords[8]))
        goto exit;

    // check the operation mode
    if (SM_INVALID == _check_sheval_mode(mode, keywords[9]))
        goto exit;

    // check the coefficients list
    if (!PyList_Check(obj_coef_list))
    {
        PyErr_Format(PyExc_ValueError, "The '%s' parameter is expected to be a list!", keywords[3]);
        goto exit;
    }

    if (PyList_Size(obj_coef_list) < 1)
    {
        PyErr_Format(PyExc_ValueError, "The '%s' list is expected to have at least one item!", keywords[3]);
        goto exit;
    }

    if (PyList_Size(obj_coef_list) > SHEVALTEMP_MAX_COEF_LIST_SIZE)
    {
        PyErr_Format(PyExc_ValueError, "The '%s' list cannot have more that %d iitems!", keywords[3], SHEVALTEMP_MAX_COEF_LIST_SIZE);
        goto exit;
    }

    // parse coefficient sets
    {
        size_t i, n = PyList_Size(obj_coef_list);

        for (i = 0; i < n; ++i)
        {
            if(parse_coefset(&coefsets[i], i, PyList_GetItem(obj_coef_list, i)))
                goto exit;

            ncoefset += 1;
        }

    }

    // cast the objects to arrays
    if (NULL == (arr_t = _get_as_double_array(obj_t, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_x = _get_as_double_array(obj_x, 1, 0, NPY_ARRAY_ALIGNED, keywords[1])))
        goto exit;

    if (degree < 0)
    {
        PyErr_Format(PyExc_ValueError, "Invalid degree %d!", degree);
        goto exit;
    }

    // check the last dimension of the position array
    if (_check_array_dim_eq(arr_x, -1, 3, keywords[0]))
        goto exit;

    // check that times match the spatial coordinates
    {
        npy_intp i, ndim = MIN(PyArray_NDIM(arr_t), PyArray_NDIM(arr_x) - 1);
        for (i = 0; i < ndim; ++i) {
            if (PyArray_DIMS(arr_t)[i] != PyArray_DIMS(arr_x)[i]) {
                PyErr_Format(PyExc_ValueError, "Shape of %s array "\
                    " does not match shape of %s!", keywords[0], keywords[1]);
                goto exit;
            }
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

    // create new output arrays
    {
        npy_intp i, ndim = MAX(PyArray_NDIM(arr_t), PyArray_NDIM(arr_x) - 1);
        npy_intp *dims_src = PyArray_NDIM(arr_t) >= PyArray_NDIM(arr_x) - 1 ? PyArray_DIMS(arr_t) : PyArray_DIMS(arr_x);
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

    // evaluate model

    if(_shevaltemp_model_init(
        &model, is_internal, degree, ct_in, ct_out,
        coefsets, ncoefset, scale_potential, scale_gradient
    ))
        goto exit;

    switch(mode)
    {
        case SM_POTENTIAL:
            _shevaltemp1(
                _array_to_arrd(arr_t),
                _array_to_arrd(arr_x),
                _array_to_arrd(arr_pot),
                 &model
            );
            retval = (PyObject*) arr_pot;
            break;

        case SM_GRADIENT:
            _shevaltemp2(
                _array_to_arrd(arr_t),
                _array_to_arrd(arr_x),
                _array_to_arrd(arr_grd),
                 &model
            );
            retval = (PyObject*) arr_grd;
            break;

        case SM_POTENTIAL_AND_GRADIENT:
            _shevaltemp3(
                _array_to_arrd(arr_t),
                _array_to_arrd(arr_x),
                _array_to_arrd(arr_pot),
                _array_to_arrd(arr_grd),
                 &model
            );
            if (NULL == (retval = Py_BuildValue("NN", (PyObject*) arr_pot, (PyObject*) arr_grd)))
                goto exit;
            break;
    }

  exit:

    _shevaltemp_model_destroy(&model);

    {
        size_t i;
        for (i = 0; i < ncoefset; ++i)
            _coefset_destroy(&coefsets[i]);
    }

    // decrease reference counters to the arrays
    if (arr_t) Py_DECREF(arr_t);
    if (arr_x) Py_DECREF(arr_x);
    if (arr_scale) Py_DECREF(arr_scale);
    if (!retval && arr_grd) Py_DECREF(arr_grd);
    if (!retval && arr_pot) Py_DECREF(arr_pot);

    return retval;
}


static int parse_coefset(COEF_SET *coefset, size_t idx, PyObject *obj_coef_list_item) {

    static char *keywords[] = {
        "coef_t", "coef_gh", "coef_nm", "spline_order", NULL
    };

    int spline_order = SHEVALTEMP_DEFAULT_SPLINE_ORDER;
    PyObject *obj_ct = NULL; // coef_t object
    PyObject *obj_cv = NULL; // coef_gh object
    PyObject *obj_cm = NULL; // coef_nm object
    PyArrayObject *arr_ct = NULL; // coef_t array
    PyArrayObject *arr_cv = NULL; // coef_gh array
    PyArrayObject *arr_cm = NULL; // coef_nm array

    static char label[128];

    // check tuple item
    if (!PyTuple_Check(obj_coef_list_item))
    {
        PyErr_Format(PyExc_ValueError, "The 'coef_set_list' item #%ld is expected to be a tuple!", idx);
        goto error;
    }

    if (PyTuple_Size(obj_coef_list_item) != 4)
    {
        PyErr_Format(PyExc_ValueError, "The 'coef_set_list' item #%ld tuple is expected to have 4 items!", idx);
        goto error;
    }

    // extract array objects
    if (NULL == (obj_ct = PyTuple_GetItem(obj_coef_list_item, 0)))
        goto error;

    if (NULL == (obj_cv = PyTuple_GetItem(obj_coef_list_item, 1)))
        goto error;

    if (NULL == (obj_cm = PyTuple_GetItem(obj_coef_list_item, 2)))
        goto error;

    // extract and parse spline order integer flag
    {
        long value = -1;
        PyObject *obj = NULL;

        if (NULL == (obj = PyTuple_GetItem(obj_coef_list_item, 3)))
            goto error;

        if (!PyLong_Check(obj))
        {
            PyErr_Format(PyExc_ValueError, "The spline order is expected to be an integer value! (#%ld)", idx);
            goto error;
        }

        value = PyLong_AsLong(obj);
        if (NULL != PyErr_Occurred())
            goto error;

        if ((value < SHEVALTEMP_MIN_SPLINE_ORDER) || (value > SHEVALTEMP_MAX_SPLINE_ORDER))
        {
            PyErr_Format(PyExc_ValueError, "Unsupported spline order %d! (#%ld)", spline_order, idx);
            goto error;
        }

        spline_order = value;
    }

    // parse coefficients array
    snprintf(label, 128, "%s#%ld", keywords[0], idx);
    if (NULL == (arr_ct = _get_as_double_array(obj_ct, 1, 1, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_IN_ARRAY, label)))
        goto error;

    snprintf(label, 128, "%s#%ld", keywords[1], idx);
    if (NULL == (arr_cv = _get_as_double_array(obj_cv, 2, 2, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_IN_ARRAY, label)))
        goto error;

    snprintf(label, 128, "%s#%ld", keywords[2], idx);
    if (NULL == (arr_cm = _get_as_int_array(obj_cm, 2, 2, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_IN_ARRAY, label)))
        goto error;

    // check dimensions of the coefficient arrays

    snprintf(label, 128, "%s#%ld", keywords[1], idx);

    if (_check_array_dim_eq(arr_cv, 0, PyArray_DIMS(arr_cm)[0], label))
        goto error;

    if (_check_array_dim_eq(arr_cv, 1, PyArray_DIMS(arr_ct)[0], label))
        goto error;

    snprintf(label, 128, "%s#%ld", keywords[2], idx);

    if (_check_array_dim_eq(arr_cm, 1, 2, label))
        goto error;

    _coefset_init(coefset, arr_ct, arr_cv, arr_cm, spline_order);

    return 0;

  error:
    if (arr_ct) Py_DECREF(arr_ct);
    if (arr_cv) Py_DECREF(arr_cv);
    if (arr_cm) Py_DECREF(arr_cm);

    return 1;
}


#endif  /* PYMM_SHEVAL_H */
