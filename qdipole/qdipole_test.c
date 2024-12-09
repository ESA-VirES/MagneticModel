/*-----------------------------------------------------------------------------
 *
 * Quasi-Dipole Magnetic Coordinates - Fortran wrapper test
 *
 * Project: EOX Magnetic Model
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2015 EOX IT Services GmbH
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

#define N 6

int main(int argc, char* argv[])
{
    int i;
    int n_data = N;
    double time[N] = { 2012.5, 2013.5, 2014.5, 2019.5, 2024.5, 2029.5 };
    double gclat[N] = { 0.0, 30.0, 75.0, -70.0, -30.0, 80.0 };
    double gclon[N] = { 0.0, 45.0, -90.0, 90.0, -45.0, 120.0 };
    double gcrad[N] = { 7000.0, 6371.0, 6371.0, 6371.0, 6371.0, 7371.0 };
    double qdlon[N], qdlat[N], xmlt[N];
    double f11[N], f12[N], f21[N], f22[N];
    const char *fname = "apexsh_1980-2030.txt";

    if (argc > 1)
    {
        fname = argv[1];
    }
    printf("Using: %s\n", fname);

    for (i = 0; i < N; ++i)
    {
        qdlon[i] = 0.0;
        qdlat[i] = 0.0;
        xmlt[i] = 0.0;
        f11[i] = 0.0;
        f12[i] = 0.0;
        f21[i] = 0.0;
        f22[i] = 0.0;
    }

    c_make_apex(qdlat, qdlon, xmlt, f11, f12, f21, f22,
               time, gcrad, gclat, gclon, n_data, fname);

    for (i = 0; i < N; ++i)
    {
        printf("#%i\n", i+1);
        printf("    time:  %g\n", time[i]);
        printf("    gclat: %g\n", gclat[i]);
        printf("    gclon: %g\n", gclon[i]);
        printf("    gcrad: %g\n", gcrad[i]);
        printf("    qdlon: %g\n", qdlon[i]);
        printf("    qdlat: %g\n", qdlat[i]);
        printf("    xmlt:  %g\n", xmlt[i]);
        printf("    f11: %g", f11[i]);
        printf("    f12: %g\n", f12[i]);
        printf("    f21: %g", f21[i]);
        printf("    f22: %g\n", f22[i]);
    }

    return 0;
}
