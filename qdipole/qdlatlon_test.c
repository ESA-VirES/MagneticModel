/*-----------------------------------------------------------------------------
 *
 * Quasi-Dipole Magnetic Coordinates - Fortran wrapper test
 *
 * Project: EOX Magnetic Model
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2016-2024 EOX IT Services GmbH
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
#include <string.h>
#include "cqdipole.h"

#define N 6

const char *basename(char const *path)
{
    const char *position = strrchr(path, '/');
    if (position != NULL)
    {
        return position + 1;
    }
    return path;
}


int main(int argc, char* argv[])
{
    int status;
    int i;
    int n_data = N;
    double time[N] = { 2012.5, 2013.5, 2014.5, 2019.5, 2024.5, 2029.5 };
    double gclat[N] = { 0.0, 30.0, 75.0, -70.0, -30.0, 80.0 };
    double gclon[N] = { 0.0, 45.0, -90.0, 90.0, -45.0, 120.0 };
    double gcrad[N] = { 7000.0, 6371.0, 6371.0, 6371.0, 6371.0, 7371.0 };
    double qdlon1[N], qdlat1[N];
    double qdlon2[N], qdlat2[N];
    double f11[N], f12[N], f21[N], f22[N], f[N];
    const char *fname = "apexsh_1980-2030.txt";

    if (argc > 1)
    {
        fname = argv[1];
    }
    printf("Using: %s\n", basename(fname));

    if (strlen(fname) > get_qdipole_max_fname_lenght())
    {
        fprintf(stderr, "ERROR: Filename is too long and exceeds the maximum allowed %d bytes! filename = %s\n", get_qdipole_max_fname_lenght(), fname);
        return 1;
    }

    for (i = 0; i < N; ++i)
    {
        qdlon1[i] = 0.0;
        qdlat1[i] = 0.0;
        qdlon2[i] = 0.0;
        qdlat2[i] = 0.0;
        f11[i] = 0.0;
        f12[i] = 0.0;
        f21[i] = 0.0;
        f22[i] = 0.0;
        f[i] = 0.0;
    }

    status = c_eval_qdlatlon(qdlat1, qdlon1, time, gcrad, gclat, gclon, n_data, fname);

    if (status) {
        fprintf(stderr, "ERROR: Call to c_eval_qdlatlon() failed with an error! error_code = %d\n", status);
        return 1;
    }

    status = c_eval_qdlatlonvb(qdlat2, qdlon2, f11, f12, f21, f22, f,
               time, gcrad, gclat, gclon, n_data, fname);

    if (status) {
        fprintf(stderr, "ERROR: Call to c_eval_qdlatlonvb() failed with an error! error_code = %d\n", status);
        return 1;
    }

    for (i = 0; i < N; ++i)
    {
        printf("#%i\n", i+1);
        printf("    time:  %g\n", time[i]);
        printf("    gclat: %g\n", gclat[i]);
        printf("    gclon: %g\n", gclon[i]);
        printf("    gcrad: %g\n", gcrad[i]);
        printf("    qdlon: %g %g\n", qdlon1[i], qdlon2[i]);
        printf("    qdlat: %g %g\n", qdlat1[i], qdlat2[i]);
        printf("    f11: %g", f11[i]);
        printf("    f12: %g\n", f12[i]);
        printf("    f21: %g", f21[i]);
        printf("    f22: %g\n", f22[i]);
        printf("    f: %g\n", f[i]);
    }

    return 0;
}

