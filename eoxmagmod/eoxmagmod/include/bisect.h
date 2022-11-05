/**
 * @file math_aux.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Bisection interval search.
 *
 * This file contains definitions of auxiliary mathematical subroutines
 * to be used all over the whole program.
 *
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
 */

#ifndef BISECT_H
#define BISECT_H 1

#include <stddef.h>
#include <math.h>

/**
 * Common function type of the bisect functions.
 */
typedef ptrdiff_t (*f_bisect)(const double x, const double* v, const size_t n);


/**
 * @brief Find interval the given value fits in using the bisection method.
 *
 * Perform bisect search of an interval from the given sorted array of
 * nodes v of size n, defining n-1 intervals.
 *
 * Function returns the index i of the interval the value x fits in so that
 * v[i] <= x < v[i+1].
 *
 * The index is set to -1 if x < v[-1] or n is 0.
 * The index is set to n-1 if x >= v[n-1] or x is Nan.
 */

static ptrdiff_t bisect_right(const double x, const double *v, const size_t n)
{
    ptrdiff_t i_low, i_high;

    if ((n < 1)||!(x >= v[0])) {
        // no intervals (n == 0), x is below the lower bound or NaN
        return -1;
    }
    if (x >= v[n-1]) {
        // x is above the upper bound
        return n - 1;
    }

    // proceed with the halving of the intervals

    i_low = 0;
    i_high = n - 2;

    while (i_low < i_high) {
        ptrdiff_t i = i_low + 1 + (i_high - i_low) / 2;
        if (x >= v[i]) {
            i_low = i;
        } else {
            i_high = i - 1;
        }
    }

    return i_low;
}


/**
 * @brief Find interval the given value fits in using the bisection method.
 *
 * Perform bisect search of an interval from the given sorted array of
 * nodes v of size n, defining n-1 intervals.
 *
 * Function returns the index i of the interval the value x fits in so that
 * v[i] < x <= v[i+1].
 *
 * The index is set to -1 if x <= 0 v[-1] or n is 0.
 * The index is set to n-1 if x > v[n-1] or x is NaN.
 */

static ptrdiff_t bisect_left(const double x, const double *v, const size_t n)
{
    ptrdiff_t i_low, i_high;

    if ((n < 1)||!(x > v[0])) {
        // no intervals (n == 0), x is below the lower bound or NaN
        return -1;
    }
    if (x > v[n-1]) {
        // x is above the upper bound
        return n - 1;
    }

    // proceed with the halving of the intervals

    i_low = 0;
    i_high = n - 2;

    while (i_low < i_high) {
        ptrdiff_t i = i_low + 1 + (i_high - i_low) / 2;
        if (x <= v[i]) {
            i_high = i - 1;
        } else {
            i_low = i;
        }
    }

    return i_high;
}

#endif /*BISECT_H*/
