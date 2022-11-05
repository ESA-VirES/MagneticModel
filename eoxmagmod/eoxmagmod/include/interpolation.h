/**
 * @file math_aux.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Data interpolation
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

#ifndef INTERPOLATION_H
#define INTERPOLATION_H 1

#include <stddef.h>
#include <bisect.h>

#define MAX_SPLINE_ORDER 6

#define CLIP(v, v_min, v_max) ((v)<(v_min)?(v_min):((v)>(v_max)?(v_max):(v)))

/**
 * Common interpolation basis structures (interval index and basis functions).
 */

typedef struct {
    ptrdiff_t i0;
    ptrdiff_t i;
    size_t order;
    double b[MAX_SPLINE_ORDER];
} INTERP_BASIS;


/**
 * @brief common interpolation basis function type
 *
 */

typedef INTERP_BASIS (*f_get_interp_basis) (const double x, const double *v, const size_t n);

/**
 * @brief common interpolation evaluation function type
 *
 */

typedef double (*f_interp_eval)(INTERP_BASIS *basis, const double *v);

/**
 * @brief Get constant interpolation basis
 *
 */

static INTERP_BASIS get_interp0_basis (const double x, const double *v, const size_t n)
{
    ptrdiff_t idx0 = bisect_right(x, v, n);
    ptrdiff_t idx = CLIP(idx0, 0, n-1);

    INTERP_BASIS basis = {idx0, idx, 1, {1.0}};

    return basis;
}

/**
 * @brief Get linear interpolation basis (interval index and basis functions)
 *
 * First, perform bisect search of an interval from the given sorted array of
 * nodes v of size n, defining n-1 intervals.
 *
 * Once found, calculate the linear interpolation basis functions.
 * For x outside the lower and upper bounds the basis functions are
 * extrapolated.
 */

static INTERP_BASIS get_interp1_basis (const double x, const double *v, const size_t n)
{
    ptrdiff_t idx0 = bisect_right(x, v, n);
    ptrdiff_t idx = CLIP(idx0, 0, n-2);
    double alpha = (x - v[idx]) / (v[idx+1] - v[idx]);

    INTERP_BASIS basis = {idx0, idx, 2, {1.0 - alpha, alpha}};

    return basis;
}

/**
 * @brief Get first derivative of linear interpolation basis (interval index
 * and basis functions)
 *
 * First, perform bisect search of an interval from the given sorted array of
 * nodes v of size n, defining n-1 intervals.
 *
 * Once found, calculate the first derivative linear interpolation basis
 * functions. For x outside the lower and upper bounds the basis functions are
 * extrapolated.
 */

static INTERP_BASIS get_interp1d1_basis (const double x, const double *v, const size_t n)
{
    ptrdiff_t idx0 = bisect_right(x, v, n);
    ptrdiff_t idx = CLIP(idx0, 0, n-2);
    double alpha = 1.0 / (v[idx+1] - v[idx]);

    INTERP_BASIS basis = {idx0, idx, 2, {-alpha, alpha}};

    return basis;
}

/**
 * @brief Evaluate constant interpolation.
 */

static double interp0_eval(INTERP_BASIS *basis, const double *v) {
    return v[basis->i];
}

/**
 * @brief Evaluate linear interpolation.
 */

static double interp1_eval(INTERP_BASIS *basis, const double *v) {
    return (basis->b[0]*v[basis->i] + basis->b[1]*v[basis->i+1]);
}

#undef CLIP

#endif /*INTERPOLATION_H*/
