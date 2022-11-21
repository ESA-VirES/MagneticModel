/**
 * @file fourier_series.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Fourier series subroutines
 *
 * Various utilities needed by the spherical harmonic model evaluation.
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

#ifndef FOURIER_SERIES_H
#define FOURIER_SERIES_H 1

/**
 * @brief Evaluate series of sines and cosines.
 *
 * Evaluate the series of sines and cosines for the given parameter.
 *
 *   sin(i*t) for i = 0...i_size
 *   cos(i*t) for i = 0...i_size
 *
 *   where 1 <= i_size
 *
 * This subroutine uses a faster evaluation based on pure recurrent
 * addition/substration and multiplication:
 *
 *   sin(i*t) = cos((i-1)*t)*sin(t) + sin((i-1)*t)*cos(t)
 *   cos(i*t) = cos((i-1)*t)*cos(t) - sin((i-1)*t)*sin(t)
 *
 * @param[out] sin_it Array of the sin(i*t) values [i_size]
 * @param[out] cos_it Array of the cos(i*t) values [i_size]
 * @param      i_size Arrays size (>= 2)
 * @param      t      Parameter value
 */

static void fs_sincos(
    double *sin_it, double *cos_it,
    const ptrdiff_t i_size, const double t
)
{
    ptrdiff_t i;
    const double sin_t = sin(t);
    const double cos_t = cos(t);
    double _sin_it = sin_t;
    double _cos_it = cos_t;

    sin_it[0] = 0.0;
    cos_it[0] = 1.0;

    sin_it[1] = _sin_it;
    cos_it[1] = _cos_it;

    for (i = 2; i <= i_size; ++i)
    {
        double _sin_it_new = _cos_it*sin_t + _sin_it*cos_t;
        double _cos_it_new = _cos_it*cos_t - _sin_it*sin_t;

        sin_it[i] = _sin_it = _sin_it_new;
        cos_it[i] = _cos_it = _cos_it_new;
    }
}


/**
 * @brief Evaluate series of sines and cosines.
 *
 * Evaluate the series of sines and cosines for the given parameter.
 *
 *   cos(i*t) for i = 0...i_size
 *   sin(i*t) for i = 0...i_size
 *
 * This subroutine contains the reference (slow) implementation evaluation
 * sine and cosine functions for each term of the series.
 *
 * @param[out] sin_it Array of the sin(i*t) values [i_size]
 * @param[out] cos_it Array of the cos(i*t) values [i_size]
 * @param      i_size Arrays dimensions (>= 1)
 * @param      t      Parameter value
 */

static void fs_sincos_ref(
    double *sin_it, double *cos_it,
    const ptrdiff_t i_size, const double t
)
{
    ptrdiff_t i;

    for (i = 0; i <= i_size; ++i)
    {
        sin_it[i] = sin(i*t);
        cos_it[i] = cos(i*t);
    }
}


/**
 * @brief Evaluate series of sines and cosines.
 *
 * Evaluate the series of sines and cosines for the given parameter.
 *
 *   sin(i*t) for i = i_min...i_max
 *   cos(i*t) for i = i_min...i_max
 *
 *   where i_min <= i_max and both can be negative.
 *
 *   The values are retrieved from a precalculated series
 *
 *   sin(i*t) for i = 0...order
 *   cos(i*t) for i = 0...order
 *
 *   where order = max(abs(i_min), abs(i_max))
 *
 *   (cos(-i*t) = cos(i*t) and sin(-i*t) = -sin(i*t))
 *
 * @param[out] sin_it     Array of the sin(i*t) values [i_max-i_min+1]
 * @param[out] cos_it     Array of the cos(i*t) values [i_max-i_min+1]
 * @param[in]  sin_it_pos Array of the sin(i*t) values for positive i [max(abs(i_min), abs(i_max))+1]
 * @param[in]  cos_it_pos Array of the cos(i*t) values for positive i [max(abs(i_min), abs(i_max))+1]
 * @param      i_min      Min. degree value
 * @param      i_max      Max. degree value
 */

static void fs_sincos_neg(
    double *sin_it, double *cos_it,
    const double *sin_it_pos, const double *cos_it_pos,
    const ptrdiff_t i_min, const ptrdiff_t i_max
)
{
    ptrdiff_t i;

    for (i = i_min; i <= i_max; ++i)
    {
        const ptrdiff_t idx = i - i_min;
        const double tmp = sin_it_pos[abs(i)];
        sin_it[idx] = (i < 0 ? -tmp : tmp);
        cos_it[idx] = cos_it_pos[abs(i)];
    }
}


/**
 * @brief Evaluate matrix of sines and cosines for the 2D Fourier series
 *
 * Evaluate matrix of sines and cosines for the 2D Fourier series
 *
 *   sin(i*t + j*s) = cos(i*t)*sin(j*s) + sin(i*t)*cos(j*s)
 *   cos(i*t + j*s) = cos(i*t)*cos(j*s) - sin(i*t)*sin(j*s)
 *
 *   The values are retrieved from a precalculated series
 *
 *   sin(i*t) for i = i_min...i_max
 *   cos(i*t) for i = i_min...i_max
 *
 *   sin(j*s) for j = j_min...j_max
 *   cos(j*s) for j = j_min...j_max
 *
 * @param[out] sin_it_js Array of the calculated sin(i*t+j*s) values [j_size, i_size]
 * @param[out] cos_it_js Array of the calculated cos(i*t+j*s) values [j_size, i_size]
 * @param[in]  sin_it    Array of the sin(i*t) values [i_size]
 * @param[in]  cos_it    Array of the cos(i*t) values [i_size]
 * @param[in]  sin_js    Array of the sin(j*s) values [j_size]
 * @param[in]  cos_js    Array of the cos(j*s) values [j_size]
 * @param      i_size    Array size
 * @param      j_size    Array size
 */

static void fs_sincos_2d(
    double *sin_it_js, double *cos_it_js,
    const double *sin_it, const double *cos_it,
    const double *sin_js, const double *cos_js,
    const ptrdiff_t i_size, const ptrdiff_t j_size
)
{
    ptrdiff_t i, j;

    for (i = 0; i < i_size; ++i)
    {
        const ptrdiff_t i_stride = i * j_size;

        for (j = 0; j < j_size; ++j)
        {
            sin_it_js[i_stride + j] = cos_it[i]*sin_js[j] + sin_it[i]*cos_js[j];
            cos_it_js[i_stride + j] = cos_it[i]*cos_js[j] - sin_it[i]*sin_js[j];
        }
    }
}


/**
 * @brief Evaluate 2D Fourier series
 *
 * Evaluate 2D Fourier series in the form of
 *
 *   f(t,s) = sum(A_ij*cos(i*t + j*s) + B_ij*sin(i*t + j*s))
 *
 *   for i = i_min...i_max and j = j_min...j_max
 *
 * @param[in] ab        Array of the 2D Fourier series coefficients [i_size, j_size, 2]
 * @param[in] sin_it_js Array of the sin(i*t+j*s) values [j_size, i_size]
 * @param[in] cos_it_js Array of the cos(i*t+j*s) values [j_size, i_size]
 * @param     i_size    Array size
 * @param     j_size    Array size
 * @return    Evaluated function f(t,s)
 */

static double fs_eval_2d(
    const double *ab, const double *sin_it_js, const double *cos_it_js,
    ptrdiff_t i_size, ptrdiff_t j_size
)
{
    double f = 0.0;
    ptrdiff_t i;
    ptrdiff_t size = i_size * j_size;

    for (i = 0; i < size; ++i)
    {
        f += ab[2*i]*cos_it_js[i] + ab[2*i+1]*sin_it_js[i];
    }

    return f;
}

#endif  /*FOURIER_SERIES_H*/
