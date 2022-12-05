/*-----------------------------------------------------------------------------
 *
 * Geomagnetic Model - C python bindings - spherical-harmonics evaluation
 * - common subroutines and definitions
 *
 * Author: Martin Paces <martin.paces@eox.at>
 *
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

#ifndef PYMM_SPHAR_COMMON_H
#define PYMM_SPHAR_COMMON_H 1

#include "pymm_aux.h"

#ifndef MAX_DEGREE
#define MAX_DEGREE 2147483647
#endif


// compare dimension from two arrays - return 1 if not equal
int _compare_dimensions(npy_intp *dims1, npy_intp *dims2, npy_intp ndim)
{
    npy_intp i;

    for (i = 0; i < ndim; ++i)
    {
        if (dims1[i] != dims2[i])
            return 1;
    }

    return 0;
}


// extract largest common array shape
void _extract_common_shape(
    int *arg_idx_out, npy_intp *ndim_out, npy_intp **dims_out,
    PyArrayObject **arrays, npy_intp *arrays_ndim, int narrays
)
{

    int arg_idx = -1;
    npy_intp ndim = 0;
    npy_intp *dims = NULL;

    if (narrays > 0)
    {
        int i;

        arg_idx = 0;
        ndim = PyArray_NDIM(arrays[0]) - arrays_ndim[0];
        dims = PyArray_DIMS(arrays[0]);

        for (i = 1; i < narrays; ++i)
        {
            npy_intp ndim_tmp = PyArray_NDIM(arrays[i]) - arrays_ndim[i];

            if (ndim < ndim_tmp)
            {
                arg_idx = i;
                ndim = ndim_tmp;
                dims = PyArray_DIMS(arrays[i]);
            }
        }
    }

    *arg_idx_out = arg_idx;
    *ndim_out = ndim;
    *dims_out = dims;
}


// extract maximum degree from an array of possible values
void _get_max_degree(
    int *arg_idx_out, npy_intp *max_degree_out, npy_intp degrees[], int ndegrees
)
{
    npy_intp max_degree = -1;
    int arg_idx = -1;

    if (ndegrees > 0)
    { 
        int i;
        arg_idx = 0;
        max_degree = degrees[0];

        for (i = 1; i < ndegrees; ++i)
        {
            if (degrees[i] < max_degree)
            {
                arg_idx = i;
                max_degree = degrees[i];
            }
        }
    }

    *arg_idx_out = arg_idx;
    *max_degree_out = max_degree;
}
    

npy_intp _size_to_degree(npy_intp size)
{
    npy_intp degree = (npy_intp)(0.5 *(sqrt(8*size + 1) - 3));
    if (((degree+3)*(degree+2))/2 == size)
        return degree + 1;
    return degree;
}

#endif /*PYMM_SPHAR_COMMON_H*/
