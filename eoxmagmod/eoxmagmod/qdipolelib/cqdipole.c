/*-----------------------------------------------------------------------------
 *
 * Qausi-Dipole Magnetic Coordinates - C interface (wrapper arround the F code)
 *
 * Project: World Magnetic Model - python interface
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

#include <string.h>
#include "cqdipole.h"

/* NOTE: the C/Fortran interface conventions may differ between various compilers */
void make_apex_(double*, double*, double*, double *, double *, double *, double*,
        const double*, const double*, const double*, const double*, 
        const int*, const char*);

void c_make_apex(
    double* qdlat, double* qdlon, double* xmlt,
    double* f11, double* f12, double* f21, double* f22,
    const double* time, const double* gcrad, const double* gclat,
    const double* gclon, int n_data, const char *fname)
{
    /* NOTE: Fotran code expects the file name as a 128-character string.*/
    char fname128[129];
    strncpy(fname128, fname, 128);
    fname128[128] = '\0' ; 

    /* call Fortran subroutine */
    make_apex_(qdlat, qdlon, xmlt, f11, f12, f21, f22,
               time, gcrad, gclat, gclon, &n_data, fname128);
}
