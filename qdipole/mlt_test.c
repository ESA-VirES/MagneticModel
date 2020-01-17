/*-----------------------------------------------------------------------------
 *
 * Magnetic Local Time - Fortran wrapper test
 *
 * Project: EOX Magnetic Model
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2016 EOX IT Services GmbH
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

#include <stdio.h>
#include "cqdipole.h"

#define N 9

int main(int argc, char* argv[])
{
    int i;
    int n_data = N;
    double time[N] = {
        0.0, 0.0, 0.0, 0.0, // MJD200 origin 2000-01-01T00:00:00
        5923.00, 5923.25, 5923.50, 5923.75, 5924.00, // spring equinox 2016
    };
    double qdlon[N] = {
        0.0, 90.0, 180.0, -90.0,
        0.0, 0.0, 0.0, 0.0, 0.0
    };
    double mlt[N];
    const char *fname = "apexsh_1980-2025.txt";

    if (argc > 1)
    {
        fname = argv[1];
    }
    printf("Using: %s\n", fname);

    for (i = 0; i < N; ++i)
    {
        mlt[i] = 0.0;
    }

    c_eval_mlt(mlt, qdlon, time, n_data, fname);

    for (i = 0; i < N; ++i)
    {
        printf("#%i\n", i+1);
        printf("    time:  %g\n", time[i]);
        printf("    qdlon: %g\n", qdlon[i]);
        printf("    mlt:   %g\n", mlt[i]);
    }

    return 0;
}
