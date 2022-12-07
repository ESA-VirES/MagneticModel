/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings
 * - coefficients evaluation from the 2D Fourier series (MIO models)
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

#ifndef PYMM_FOURIER2D_H
#define PYMM_FOURIER2D_H

#include "math.h"
#include "fourier_series.h"
#include "pymm_aux.h"

#ifndef MIN
#define MIN(a,b) ((b)<(a)?(b):(a))
#endif

#ifndef MAX
#define MAX(a,b) ((b)>(a)?(b):(a))
#endif


/* 2D Fourier series - auxiliary structure */
typedef struct Fourier2D {
    ptrdiff_t min_degree1;
    ptrdiff_t max_degree1;
    ptrdiff_t min_degree2;
    ptrdiff_t max_degree2;
    ptrdiff_t max_abs_degree;
    double scale1;
    double scale2;
    double (*cos_sin_tmp)[2];
    double (*cos_sin_1)[2];
    double (*cos_sin_2)[2];
    double (*cos_sin_12)[2];
} FOURIER_2D;

static void _fourier2d_reset(FOURIER_2D *f2d);
static void _fourier2d_destroy(FOURIER_2D *f2d);
static int _fourier2d_init(
    FOURIER_2D *f2d,
    const ptrdiff_t min_degree1,
    const ptrdiff_t size1,
    const ptrdiff_t min_degree2,
    const ptrdiff_t size2,
    const double scale1,
    const double scale2
);


static void _fourier2d_eval_coeff(
    ARRAY_DATA *arrd_t1,
    ARRAY_DATA *arrd_t2,
    ARRAY_DATA *arrd_c,
    ARRAY_DATA *arrd_c0,
    FOURIER_2D *f2d
);

/* Python function definition */

#define DOC_FOURIER2D "\n"\
"   coeff = fourier2d(time1, time2, coeff0, min_degree1=0, min_degree2=0, scale1=1.0, scale2=1.0)\n"\
"\n"\
"     Interpolate time series of SH coefficients.\n"\
"\n"\
"     The input parameters are:\n"\
"       time1 - seasonal time variable\n"\
"       time2 - diurnal time variable\n"\
"       coeff0[...,n1,n2,2] - series of the 2D Fourier series coefficients\n"\
"       min_degree1 - seasonal time minimal Fourier series min. degree\n"\
"               (max. degree is determined from the coeff0 dimension n1)\n"\
"       min_degree2 - diurnal time minimal Fourier series min. degree\n"\
"               (max. degree is determined from the coeff0 dimension n2)\n"\
"       scale1 - seasonal time scaling factor\n"\
"       scale2 - diurnal time scaling factor\n"\
"\n"\
"     Output:\n"\
"       coeff[...] - coefficients at the given diurnal and seasonal times\n"\
"\n"

static PyObject* fourier2d(PyObject *self, PyObject *args, PyObject *kwdict)
{
    static char *keywords[] = {
        "time1", "time2", "coeff", "min_degree1", "min_degree2",
        "scale1", "scale2", NULL
    };

    int min_degree1 = 0; // dimension #1 - min. degree
    int min_degree2 = 0; // dimension #2 - min. degree
    double scale1 = 1.0; // time #1 scaling factor
    double scale2 = 1.0; // time #2 scaling factor

    PyArrayObject *arr_t1 = NULL; // time #1 (dimension #1 variable) (input)
    PyArrayObject *arr_t2 = NULL; // time #2 (dimension #2 variable) (input)
    PyArrayObject *arr_c0 = NULL; // 2D Fourier series coefficients (input)
    PyArrayObject *arr_c = NULL; // calculated coefficients (output)

    PyObject *obj_t1 = NULL; // input object
    PyObject *obj_t2 = NULL; // input object
    PyObject *obj_c0 = NULL; // input object
    PyObject *retval = NULL; // output value

    FOURIER_2D f2d = {0};

    // parse input arguments
    if (!PyArg_ParseTupleAndKeywords(
        args, kwdict, "OOO|iidd:fourier2d", keywords,
        &obj_t1, &obj_t2, &obj_c0, &min_degree1, &min_degree2, &scale1, &scale2
    ))
        goto exit;

    // cast the input objects to arrays
    if (NULL == (arr_t1 = _get_as_double_array(obj_t1, 0, 0, NPY_ARRAY_ALIGNED, keywords[0])))
        goto exit;

    if (NULL == (arr_t2 = _get_as_double_array(obj_t2, 0, 0, NPY_ARRAY_ALIGNED, keywords[1])))
        goto exit;

    if (NULL == (arr_c0 = _get_as_double_array(obj_c0, 3, 0, NPY_ARRAY_C_CONTIGUOUS|NPY_ARRAY_ALIGNED, keywords[2])))
        goto exit;

    // check the array dimensions
    if (PyArray_NDIM(arr_t1) != PyArray_NDIM(arr_t2))
    {
        PyErr_Format(PyExc_ValueError, "Mismatch of %s and %s dimensions!", keywords[0], keywords[1]);
        goto exit;
    }

    {
        npy_intp i;

        for (i = 0; i < PyArray_NDIM(arr_t1); ++i)
        {
            if (PyArray_DIMS(arr_t1)[i] != PyArray_DIMS(arr_t2)[i])
            {
                PyErr_Format(PyExc_ValueError, "Mismatch of %s and %s shapes!", keywords[0], keywords[1]);
                goto exit;
            }
        }
    }

    if (PyArray_DIMS(arr_c0)[PyArray_NDIM(arr_c0)-1] != 2) {
        PyErr_Format(PyExc_ValueError, "Incorrect shape of %s!", keywords[2]);
        goto exit;
    }

    // create the output array
    {
        npy_intp ndim = PyArray_NDIM(arr_t1) + PyArray_NDIM(arr_c0) - 3;
        npy_intp dims[ndim];
        npy_intp i;

        for (i = 0; i < PyArray_NDIM(arr_t1); ++i)
        {
            dims[i] = PyArray_DIMS(arr_t1)[i];
        }

        for (i = 0; i < PyArray_NDIM(arr_c0) - 3; ++i)
        {
            dims[i + PyArray_NDIM(arr_t1)] = PyArray_DIMS(arr_c0)[i];
        }

        if (NULL == (arr_c = _get_new_array(ndim, dims, NPY_DOUBLE))) {
            goto exit;
        }
    }

    // evaluate coefficients
    {
        npy_intp size1 = PyArray_DIMS(arr_c0)[PyArray_NDIM(arr_c0)-3];
        npy_intp size2 = PyArray_DIMS(arr_c0)[PyArray_NDIM(arr_c0)-2];

        if (_fourier2d_init(&f2d, min_degree1, size1, min_degree2, size2, scale1, scale2))
            goto exit;

        ARRAY_DATA arrd_t1 = _array_to_arrd(arr_t1);
        ARRAY_DATA arrd_t2 = _array_to_arrd(arr_t2);
        ARRAY_DATA arrd_c = _array_to_arrd(arr_c);
        ARRAY_DATA arrd_c0 = _array_to_arrd(arr_c0);

        _fourier2d_eval_coeff(&arrd_t1, &arrd_t2, &arrd_c, &arrd_c0, &f2d);
    }

    // get return value
    retval = (PyObject*) arr_c;

  exit:

    _fourier2d_destroy(&f2d);

    // decrease reference counters to the arrays
    if (arr_t1) Py_DECREF(arr_t1);
    if (arr_t2) Py_DECREF(arr_t2);
    if (arr_c0) Py_DECREF(arr_c0);
    if (!retval && arr_c) Py_DECREF(arr_c);

    return retval;
}


/*
 * high level nD-array recursive coefficient 2D Fourier expansion
 */
static void _fourier2d_cossin(FOURIER_2D *f2d, const double t1, const double t2);
static double _fourier2d_eval(FOURIER_2D *f2d, const double (*ab)[2]);
static void _fourier2d_eval_coeff_single_time(
    ARRAY_DATA *arrd_c,
    ARRAY_DATA *arrd_c0,
    FOURIER_2D *f2d
);


static void _fourier2d_eval_coeff(
    ARRAY_DATA *arrd_t1,
    ARRAY_DATA *arrd_t2,
    ARRAY_DATA *arrd_c,
    ARRAY_DATA *arrd_c0,
    FOURIER_2D *f2d
)
{

    if (arrd_t1->ndim > 0)
    {
        npy_intp i, n = arrd_t1->dim[0];

        for(i = 0; i < n; ++i)
        {
            ARRAY_DATA arrd_t1_item = _get_arrd_item_nocheck(arrd_t1, i);
            ARRAY_DATA arrd_t2_item = _get_arrd_item_nocheck(arrd_t2, i);
            ARRAY_DATA arrd_c_item = _get_arrd_item_nocheck(arrd_c, i);

            _fourier2d_eval_coeff(
                &arrd_t1_item,
                &arrd_t2_item,
                &arrd_c_item,
                arrd_c0,
                f2d
            );
        }
    }
    else
    {
        const double t1 = *((double*)arrd_t1->data);
        const double t2 = *((double*)arrd_t2->data);

        _fourier2d_cossin(f2d, t1, t2);

        _fourier2d_eval_coeff_single_time(arrd_c, arrd_c0, f2d);
    }
}


static void _fourier2d_eval_coeff_single_time(
    ARRAY_DATA *arrd_c,
    ARRAY_DATA *arrd_c0,
    FOURIER_2D *f2d
) {
    if (arrd_c->ndim > 0)
    {
        npy_intp i, n = arrd_c->dim[0];

        for(i = 0; i < n; ++i)
        {
            ARRAY_DATA arrd_c_item = _get_arrd_item_nocheck(arrd_c, i);
            ARRAY_DATA arrd_c0_item = _get_arrd_item_nocheck(arrd_c0, i);

            _fourier2d_eval_coeff_single_time(&arrd_c_item, &arrd_c0_item, f2d);
        }
    }
    else
    {
        double *c = ((double*)arrd_c->data);
        const double (*ab)[2] = ((double(*)[2])arrd_c0->data);

        *c = _fourier2d_eval(f2d, ab);
    }
}


/* 2D Fourier series - reset */

static void _fourier2d_reset(FOURIER_2D *f2d)
{
    memset(f2d, 0, sizeof(FOURIER_2D));
}

/* 2D Fourier series - destruction */

static void _fourier2d_destroy(FOURIER_2D *f2d)
{
    if(NULL != f2d->cos_sin_tmp) free(f2d->cos_sin_tmp);
    if(NULL != f2d->cos_sin_1) free(f2d->cos_sin_1);
    if(NULL != f2d->cos_sin_2) free(f2d->cos_sin_2);
    if(NULL != f2d->cos_sin_12) free(f2d->cos_sin_12);

    _fourier2d_reset(f2d);
}

/* 2D Fourier series - initialization */

static int _fourier2d_init(
    FOURIER_2D *f2d,
    const ptrdiff_t min_degree1,
    const ptrdiff_t size1,
    const ptrdiff_t min_degree2,
    const ptrdiff_t size2,
    const double scale1,
    const double scale2
)
{
    ptrdiff_t max_degree1 = min_degree1 + size1 - 1;
    ptrdiff_t max_degree2 = min_degree2 + size2 - 1;
    ptrdiff_t max_abs_degree;

    _fourier2d_reset(f2d);

    f2d->min_degree1 = min_degree1;
    f2d->max_degree1 = max_degree1;
    f2d->min_degree2 = min_degree2;
    f2d->max_degree2 = max_degree2;
    f2d->scale1 = scale1;
    f2d->scale2 = scale2;

    max_abs_degree = 1;
    max_abs_degree = MAX(max_abs_degree, abs(min_degree1));
    max_abs_degree = MAX(max_abs_degree, abs(max_degree1));
    max_abs_degree = MAX(max_abs_degree, abs(min_degree1));
    max_abs_degree = MAX(max_abs_degree, abs(max_degree2));
    f2d->max_abs_degree = max_abs_degree;

    if (NULL == (f2d->cos_sin_tmp = (double(*)[2])malloc((max_abs_degree+1)*sizeof(double[2]))))
    {
        PyErr_Format(PyExc_MemoryError, "_fourier2d_reset: cos_sin_tmp");
        goto error;
    }

    if (NULL == (f2d->cos_sin_1 = (double(*)[2])malloc(size1*sizeof(double[2]))))
    {
        PyErr_Format(PyExc_MemoryError, "_fourier2d_reset: cos_sin_1");
        goto error;
    }

    if (NULL == (f2d->cos_sin_2 = (double(*)[2])malloc(size2*sizeof(double[2]))))
    {
        PyErr_Format(PyExc_MemoryError, "_fourier2d_reset: cos_sin_2");
        goto error;
    }

    if (NULL == (f2d->cos_sin_12 = (double(*)[2])malloc(size1*size2*sizeof(double[2]))))
    {
        PyErr_Format(PyExc_MemoryError, "_fourier2d_reset: cos_sin_12");
        goto error;
    }

    return 0;

  error:
    _fourier2d_destroy(f2d);
    return 1;
}

/* 2D Fourier series - evaluate series from the coefficients and sin/cos matrix */

static double _fourier2d_eval(FOURIER_2D *f2d, const double (*ab)[2])
{
    return fs_eval_2d(
        ab, f2d->cos_sin_12,
        f2d->max_degree1 - f2d->min_degree1 + 1,
        f2d->max_degree2 - f2d->min_degree2 + 1
    );
}

/* 2D Fourier series - evaluate sin/cos matrix for the given coordinates */

static void _fourier2d_cossin(FOURIER_2D *f2d, double t1, double t2)
{
        fs_cos_sin(
            f2d->cos_sin_tmp,
            MAX(abs(f2d->min_degree1), abs(f2d->max_degree1)),
            t1 * f2d->scale1
        );

        fs_cos_sin_neg(
            f2d->cos_sin_1, f2d->cos_sin_tmp,
            f2d->min_degree1, f2d->max_degree1
        );

        fs_cos_sin(
            f2d->cos_sin_tmp,
            MAX(abs(f2d->min_degree2), abs(f2d->max_degree2)),
            t2 * f2d->scale2
        );

        fs_cos_sin_neg(
            f2d->cos_sin_2, f2d->cos_sin_tmp,
            f2d->min_degree2, f2d->max_degree2
        );

        fs_cos_sin_2d(
            f2d->cos_sin_12, f2d->cos_sin_1, f2d->cos_sin_2,
            f2d->max_degree1 - f2d->min_degree1 + 1,
            f2d->max_degree2 - f2d->min_degree2 + 1
        );
}

#endif /*PYMM_FOURIER2D_H*/
