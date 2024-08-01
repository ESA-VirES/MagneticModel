/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - auxiliary functions
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

#ifndef PYMM_AUX_H
#define PYMM_AUX_H

/*
 * parse signed long integer value
 */
static int _parse_long_value(long *target, PyObject *obj, const char *label)
{
    if (obj == NULL)
        return 1;

    if (!PyLong_Check(obj))
    {
        PyErr_Format(PyExc_ValueError, "The %s parameter is expected to be an integer value!", label);
        return 1;
    }

    *target = PyLong_AsLong(obj);
    if (NULL != PyErr_Occurred())
        return 1;

    return 0;
}

/*
 * parse signed integer value
 */

static int _parse_int_value(int *target, PyObject *obj, const char *label)
{
    long value;

    if (_parse_long_value(&value, obj, label))
        return 1;

    if ((long)((int)value) != value) {
        PyErr_Format(PyExc_ValueError, "The %s parameter cannot be safely cast to the requested integer type!", label);
        return 1;
    }

    *target = (int)value;

    return 0;
}


/*
 * parse double precision floating point value
 */
static int _parse_double_value(double *target, PyObject *obj, const char *label)
{
    if (obj == NULL)
        return 1;

    if (!PyFloat_Check(obj))
    {
        PyErr_Format(PyExc_ValueError, "The %s parameter is expected to be an float value!", label);
        return 1;
    }

    *target = PyFloat_AsDouble(obj);
    if (NULL != PyErr_Occurred())
        return 1;

    return 0;
}

/*
 * Check the input python object and convert it to a NumPy array of a speciefied
 * type, ensuring the native byte-order.Returns NULL if the conversion failed.
 */

static PyArrayObject* _get_as_array(PyObject *data, int typenum,
        int dmin, int dmax, int reqs, const char *label)
{
    PyArray_Descr *dtype = PyArray_DescrFromType(typenum);
    return (PyArrayObject*) PyArray_FromAny(data, dtype, dmin, dmax, reqs, NULL);
}

/*
 * Check the input python object and convert it to an integer NumPy
 * array ensuring the native byte-order.
 * Returns NULL if the conversion failed.
 */

static PyArrayObject* _get_as_int_array(PyObject *data, int dmin, int dmax,
                int reqs, const char *label)
{
    return _get_as_array(data, NPY_INT32, dmin, dmax, reqs, label);
}


/*
 * Check the input python object and convert it to a double precision NumPy
 * array ensuring the native byte-order.
 * Returns NULL if the conversion failed.
 */

static PyArrayObject* _get_as_double_array(PyObject *data, int dmin, int dmax,
                int reqs, const char *label)
{
    return _get_as_array(data, NPY_FLOAT64, dmin, dmax, reqs, label);
}


/*
 * Get new allocated NumPy array.
 */
static PyArrayObject* _get_new_array(npy_intp ndim, const npy_intp *dims, int typenum) {
    return (PyArrayObject*) PyArray_EMPTY(ndim, dims, typenum, 0);
}


/*
 * Get new allocated NumPy double precision float array. The first (N-1)
 * dimensions as read from the array of dimensions (allowing easily set the
 * same shape as the input matrix). The last Nth dimension is overridden by
 * the 'dim_last' value.
 */

static PyArrayObject* _get_new_double_array(npy_intp ndim, const npy_intp *dims, npy_intp dim_last)
{
    npy_intp i;
    npy_intp dims_new[MAX_OUT_ARRAY_NDIM];

    if (ndim > MAX_OUT_ARRAY_NDIM)
    {
        PyErr_Format(PyExc_ValueError, "The output array dimension "\
            " %d exceeds the allowed maximum value %d!", (int)ndim, MAX_OUT_ARRAY_NDIM);
        return NULL;
    }

    for (i = 0; i < (ndim-1); ++i)
        dims_new[i] = dims[i];
    if (ndim >= 1)
        dims_new[ndim-1] = dim_last;

    return _get_new_array(ndim, dims_new, NPY_DOUBLE);
}


/*
 * Check that the dimension of the array has the required value.
 */

static int _check_array_dim_eq(PyArrayObject *arr, int dim, size_t size, const char *label)
{
    if (dim < 0)
        dim += PyArray_NDIM(arr);
    int rv = PyArray_DIM(arr, dim) != size;
    if (rv)
        PyErr_Format(PyExc_ValueError, "The dimension #%d of '%s'"\
            " %ld is not equal to the required value %ld!", dim, label,
            (size_t)PyArray_DIM(arr, dim), size);
    return rv;
}


/*
 * Check that the dimension of the array is greater then or equal to the required value.
 */

static int _check_array_dim_le(PyArrayObject *arr, int dim, size_t size, const char *label)
{
    if (dim < 0)
        dim += PyArray_NDIM(arr);
    int rv = PyArray_DIM(arr, dim) < size;
    if (rv)
        PyErr_Format(PyExc_ValueError, "The dimension #%d of '%s'"\
            " %ld is lower than the minimum allowed value %ld!", dim, label,
            (size_t)PyArray_DIM(arr, dim), size);
    return rv;
}


/*
 * Check that the array dimensions match the required values
 */

static int _check_arr_dims_all_eq(PyArrayObject *arr, npy_intp ndim, const npy_intp *dims, const char *label)
{
    npy_intp dim;

    if (PyArray_NDIM(arr) != ndim)
    {
        PyErr_Format(PyExc_ValueError, "The number of dimensions of '%s'"\
            " %ld does not match the required value %ld!", label,
            (size_t)(PyArray_NDIM(arr)), (size_t)ndim);
        return 1;
    }

    for (dim = 0; dim < ndim; ++dim)
    {
        if (PyArray_DIM(arr, dim) != dims[dim])
        {
            PyErr_Format(PyExc_ValueError, "The dimensions #%ld of '%s'"\
            " %ld does not match the required value %ld!", (size_t)dim, label,
            (size_t)PyArray_DIM(arr, dim), (size_t)dims[dim]);
            return 1;
        }
    }

    return 0;
}


/*
 * Check that two arrays have the same shape
 */

static int _check_equal_shape(PyArrayObject *arr0, PyArrayObject *arr1, const char *label0, const char *label1)
{
    npy_intp dim;

    if (PyArray_NDIM(arr0) != PyArray_NDIM(arr1))
    {
        PyErr_Format(
            PyExc_ValueError, "Array dimension mismatch between '%s' and '%s'!",
            label0, label1
        );
        return 1;
    }

    for (dim = 0; dim < PyArray_NDIM(arr0); ++dim)
    {
        if (PyArray_DIM(arr0, dim) != PyArray_DIM(arr1, dim))
        {
            PyErr_Format(
                PyExc_ValueError, "Array shape mismatch between '%s' and '%s'!",
                label0, label1
            );
            return 1;
        }
    }

    return 0;
}


/*
 *  Extract 1D Numpy array or broadcast scalar to a 1D C array.
 */

static int _extract_1d_double_array(double *out, size_t size, PyArrayObject *arr_in, const char *label)
{
    npy_intp ndim = PyArray_NDIM(arr_in);
    void *data = PyArray_DATA(arr_in);

    if (ndim == 0)
    {
        size_t i;
        for (i = 0 ; i < size; ++i)
            out[i] = *((double*)data);
    }
    else if ((ndim == 1) && (PyArray_DIM(arr_in, 0) == size))
    {
        npy_intp stride = PyArray_STRIDE(arr_in, 0);
        size_t i;
        for (i = 0 ; i < size; ++i)
            out[i] = *((double*)(data + i*stride));
    }
    else
    {
        PyErr_Format(PyExc_ValueError, "Invalid %s dimension!", label);
        return 1;
    }

    return 0;
}


/*
 * Extraction of the lower dimensional parts of the arrays.
 */

typedef struct {
    void *data;
    npy_intp ndim;
    const npy_intp *dim;
    const npy_intp *stride;
} ARRAY_DATA;

static ARRAY_DATA _array_to_arrd(PyArrayObject *arr)
{
    ARRAY_DATA arrd = {
        PyArray_DATA(arr),
        PyArray_NDIM(arr),
        PyArray_DIMS(arr),
        PyArray_STRIDES(arr)
    };
    return arrd;
}

static ARRAY_DATA _get_arrd_item_nocheck(const ARRAY_DATA *arrd, npy_intp idx) {
    // extract sub-dimension
    ARRAY_DATA arrd_sub = {
        arrd->data + idx*arrd->stride[0],
        arrd->ndim - 1,
        arrd->dim + 1,
        arrd->stride + 1
    };
    return arrd_sub;
}

static ARRAY_DATA _get_arrd_item_with_guard(const ARRAY_DATA *arrd, npy_intp idx, npy_intp guard)
{
    if (arrd->ndim <= guard)
        return *arrd;

    return _get_arrd_item_nocheck(arrd, idx);
}


static ARRAY_DATA _get_arrd_item(const ARRAY_DATA *arrd, npy_intp idx)
{
    return _get_arrd_item_with_guard(arrd, idx, 0);
}

static ARRAY_DATA _get_arrd_vector_item(const ARRAY_DATA *arrd, npy_intp idx)
{
    return _get_arrd_item_with_guard(arrd, idx, 1);
}

#endif  /* PYMM_AUX_H */
