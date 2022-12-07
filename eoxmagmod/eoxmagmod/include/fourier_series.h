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

#include <stddef.h>

/**
 * @brief Evaluate series of sines and cosines (fast implementation).
 *
 * Evaluate the series of cosines and sines for the given parameter.
 *
 *   (cos(i*x), sin(i*x)) for i = 0...n
 *
 *   where 0 <= n
 *
 * This subroutine uses a faster evaluation based on pure recurrent
 * addition/substration and multiplication:
 *
 *   sin(i*x) = cos((i-1)*x)*sin(x) + sin((i-1)*x)*cos(x)
 *   cos(i*x) = cos((i-1)*x)*cos(x) - sin((i-1)*x)*sin(x)
 *
 * @param[out] cos_sin_ix Array of the sin(i*x) values [(n+1), 2]
 * @param      n number of terms
 * @param      t parameter value
 */

static void fs_cos_sin(
    double (*restrict cos_sin_ix)[2],
    const size_t n,
    const double x
)
{
    if (n >= 0)
    {
        cos_sin_ix[0][0] = 1.0;
        cos_sin_ix[0][1] = 0.0;

        if (n >= 1)
        {
            size_t i;
            const double cos_x = cos(x);
            const double sin_x = sin(x);
            double cos_ix = cos_x;
            double sin_ix = sin_x;

            cos_sin_ix[1][0] = cos_ix;
            cos_sin_ix[1][1] = sin_ix;

            for (i = 2; i <= n; i++)
            {
                const double cos_ix_new = cos_ix*cos_x - sin_ix*sin_x;
                const double sin_ix_new = cos_ix*sin_x + sin_ix*cos_x;

                cos_sin_ix[i][0] = cos_ix = cos_ix_new;
                cos_sin_ix[i][1] = sin_ix = sin_ix_new;
            }
        }
    }
}


/**
 * @brief Evaluate series of cosines and sines (reference implementation).
 *
 * Evaluate the series of cosines and sines for the given parameter.
 *
 *   (cos(i*x), sin(i*x)) for i = 0...degree
 *
 *   where 0 <= degree
 *
 * This subroutine contains the reference (slow) implementation evaluation
 * cosine and sine functions for each term of the series.
 *
 * @param[out] cos_sin_it Array of the sin(i*t) values [(degree+1), 2]
 * @param      n Number of terms
 * @param      t Parameter value
 */

static void fs_cos_sin_ref(
    double (*restrict cos_sin_it)[2],
    const size_t n,
    const double t
)
{
    size_t i;

    for (i = 0; i <= n; ++i)
    {
        cos_sin_it[i][0] = cos(i*t);
        cos_sin_it[i][1] = sin(i*t);
    }
}


/**
 * @brief Evaluate series of cosined and sines.
 *
 * Evaluate the series of sines and cosines for the given parameter.
 *
 *   (cos(i*x), sin(i*x)) for i = n_min...n_max
 *
 *   where n_min <= n_max and both can be negative.
 *
 *   The values are retrieved from a precalculated series
 *
 *   (cos(i*x), sin(i*x)) for i = 0...n
 *
 *   where n = max(abs(n_min), abs(n_max))
 *
 *   cos(-i*x) =  cos(i*x)
 *   sin(-i*x) = -sin(i*x)
 *
 * @param[out] cos_sin_ix  Array of the (cos(i*x), sin(i*x)) pairs [n_max-n_min+1, 2]
 * @param[in]  cos_sin_ix_positive  Array of the (cos(i*x), sin(i*x)) pairs for positive i [max(abs(n_min), abs(n_max))+1, 2]
 * @param      n_min  Min. term of the series.
 * @param      n_max  Max. term of the series.
 */

static void fs_cos_sin_neg(
    double (*restrict cos_sin_ix)[2],
    const double (*restrict cos_sin_ix_positive)[2],
    const ptrdiff_t n_min, const ptrdiff_t n_max
)
{
    ptrdiff_t i;

    for (i = n_min; i <= n_max; ++i)
    {
        const ptrdiff_t idx_positive = abs(i);
        const ptrdiff_t idx = (i - n_min);

        const double cos_ix_positive = cos_sin_ix_positive[idx_positive][0];
        const double sin_ix_positive = cos_sin_ix_positive[idx_positive][1];

        cos_sin_ix[idx][0] = cos_ix_positive;
        cos_sin_ix[idx][1] = (i < 0 ? -sin_ix_positive : sin_ix_positive);
    }
}


/**
 * @brief Evaluate matrix of cosines and sines for the 2D Fourier series
 *
 * Evaluate matrix of sines and cosines for the 2D Fourier series
 *
 *   (cos(i*x + j*y), sin(i*x + j*y)) for i = n_min...n_max and j = m_min...m_max
 *
 *   cos(i*x + j*y) = cos(i*x)*cos(j*y) - sin(i*x)*sin(j*y)
 *   sin(i*x + j*y) = cos(i*x)*sin(j*y) + sin(i*x)*cos(j*y)
 *
 *   The values are retrieved from a precalculated series
 *
 *   (cos(i*x), sin(i*x)) for i = i_min...i_max
 *   (cos(j*y), sin(j*y)) for j = j_min...j_max
 *
 * @param[out] cos_sin_ix_jy Array of the calculated cos/sin values [n*m, 2]
 * @param[in]  cos_sin_ix Array of the sin(i*x) values [n, 2]
 * @param[in]  cos_sin_jy Array of the sin(j*y) values [m, 2]
 * @param      n_size Array size
 * @param      m_size Array size
 */

static void fs_cos_sin_2d(
    double (*restrict cos_sin_ix_jy)[2],
    const double (*restrict cos_sin_ix)[2],
    const double (*restrict cos_sin_jy)[2],
    const size_t n,
    const size_t m
)
{
    size_t i, j;

    for (i = 0; i < n; i++)
    {
        const size_t i_stride = i * m;

        for (j = 0; j < m; j++)
        {
            const size_t idx = i_stride + j;

            cos_sin_ix_jy[idx][0] = cos_sin_ix[i][0]*cos_sin_jy[j][0] - cos_sin_ix[i][1]*cos_sin_jy[j][1];
            cos_sin_ix_jy[idx][1] = cos_sin_ix[i][0]*cos_sin_jy[j][1] + cos_sin_ix[i][1]*cos_sin_jy[j][0];
        }
    }
}

/**
 * @brief Evaluate 2D Fourier series
 *
 * Evaluate 2D Fourier series in the form of
 *
 *   f(x,y) = sum(A_ij*cos(i*x + j*y) + B_ij*sin(i*x + j*y))
 *
 *   for i = n_min...n_max and j = m_min...m_max
 *
 * @param[in] ab Array of the 2D Fourier series coefficients [n, m, 2]
 * @param[in] cos_sin_ix_jy Array of the (cos(i*x+j*y), sin(i*x+j*y)) values [m, n, 2]
 * @param     n Array size n = n_max - n_min + 1
 * @param     m Array size m = n_max - n_min + 1
 * @return    Evaluated function f(x,y)
 */

static double fs_eval_2d(
    const double (*restrict ab)[2], const double (*restrict cos_sin_ix_jy)[2],
    size_t n, size_t m
)
{
    double f = 0.0;
    size_t i;
    const size_t size = n*m;

    for (i = 0; i < size; ++i)
    {
        f += ab[i][0]*cos_sin_ix_jy[i][0] + ab[i][1]*cos_sin_ix_jy[i][1];
    }

    return f;
}

#endif  /*FOURIER_SERIES_H*/
