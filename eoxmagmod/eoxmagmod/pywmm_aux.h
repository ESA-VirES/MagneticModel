/*-----------------------------------------------------------------------------
 *
 * World Magnetic Model - C python bindings - auxiliary functions
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

#ifndef PYWMM_AUX_H
#define PYWMM_AUX_H

/*
 * Check the input python object and convert it to a double precision NumPy
 * array ensuring the native byte-order.
 */
static PyObject* _get_as_double_array(PyObject *data, int dmin, int dmax,
                int reqs, const char *label)
{
    PyArray_Descr *dtype = PyArray_DescrFromType(NPY_FLOAT64);
    PyObject *arr = PyArray_FromAny(data, dtype, dmin, dmax, reqs, NULL);
    /*
    if (NULL == arr)
        PyErr_Format(PyExc_ValueError, "Failed to cast %s to an array!", label);
    */
    return arr;
}

/*
 * Get new allocated NumPy array. The first (N-1) dimenstions as read from
 * the array of dimenstions (allowing easily set the same shape as the input
 * matrix). The last Nth dimension is overriden by the 'dim_last' value.
 */
static PyObject* _get_new_double_array(npy_intp ndim, const npy_intp *dims, npy_intp dim_last)
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
    dims_new[ndim-1] = dim_last;

    return PyArray_EMPTY(ndim, dims_new, NPY_DOUBLE, 0);
}

/*
 * Check that the dimension of the array has the required value.
 */
static int _check_array_dim_eq(PyObject *arr, int dim, size_t size, const char *label)
{
    if (dim < 0)
        dim += PyArray_NDIM(arr);
    int rv = PyArray_DIM(arr, dim) != size;
    if (rv)
        PyErr_Format(PyExc_ValueError, "The dimension #%d of '%s'"\
            " %ld is not equal the allowed value %ld!", dim, label,
            (size_t)PyArray_DIM(arr, dim), size);
    return rv;
}

static int _check_array_dim_le(PyObject *arr, int dim, size_t size, const char *label)
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
 * Extraction of the lower dimensional parts of the arrays.
 */

typedef struct {
    void *data;
    npy_intp ndim;
    const npy_intp *dim;
    const npy_intp *stride;
} ARRAY_DATA;

static ARRAY_DATA _array_to_arrd(PyObject *arr)
{
    ARRAY_DATA arrd = {
        PyArray_DATA(arr),
        PyArray_NDIM(arr),
        PyArray_DIMS(arr),
        PyArray_STRIDES(arr)
    };
    return arrd;
}

static ARRAY_DATA _get_arrd_item(const ARRAY_DATA *arrd, npy_intp idx)
{
    ARRAY_DATA arrd_sub = {
        arrd->data + idx*arrd->stride[0],
        arrd->ndim - 1,
        arrd->dim + 1,
        arrd->stride + 1
    };
    return arrd_sub;
}

#endif  /* PYWMM_AUX_H */
