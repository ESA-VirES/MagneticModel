/*-----------------------------------------------------------------------------
 *
 * Quasi-Dipole Magnetic Coordinates - C wrapper around the Fortran code
 *
 * Project: EOX Magnetic Model
 * Author: Martin Paces <martin.paces@eox.at>
 *
 *-----------------------------------------------------------------------------
 * Copyright (C) 2015-2024 EOX IT Services GmbH
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

/*
 * NOTE: The C/Fortran interface conventions may differ between various
 *       compilers. This code uses the GCC conventions.
*/

#include <string.h>
#include "cqdipole.h"
#include "qdipole_conf.h"

#define MAX_PATH_LENGTH 1023


static int check_string_lenght(const char *srt, const size_t size)
{
    return strlen(srt) > size;
}


static void copy_string(char *dst, const char *src, const size_t size)
{
    // copy string truncated to the given destination size
    strncpy(dst, src, size - 1);
    dst[size - 1] = '\0';
}


size_t get_qdipole_max_fname_lenght()
{
    return MAX_PATH_LENGTH;
}


const char* get_qdipole_version()
{
    return PACKAGE_VERSION;
}


void eval_mlt_(double*, const double*, const double*, const int*);

int c_eval_mlt(double *t_mlt, const double *qdlon, const double *t_mjd2k,
                const int n_data)
{
    /* call the Fortran subroutine */
    eval_mlt_(t_mlt, qdlon, t_mjd2k, &n_data);

    return 0;
}


void eval_subsol_(double*, double*, const double*, const int*);

int c_eval_subsol(double *sbsllat, double *sbsllon, const double *time_mjd2k,
                   const int n_data)
{
    /* call the Fortran subroutine */
    eval_subsol_(sbsllat, sbsllon, time_mjd2k, &n_data);

    return 0;
}


void eval_qdlatlon_(double*, double*, const double*, const double*,
                    const double*, const double*, const int*, const char*);

int c_eval_qdlatlon(
     double *qdlat, double *qdlon, const double *time_dy, const double *gcrad,
     const double *gclat, const double *gclon, const int n_data,
     const char *coeff_file)
{
    /* NOTE: Fortran code expects the file name as a sized size string. */
    char filename[MAX_PATH_LENGTH + 1];
    if (check_string_lenght(coeff_file, MAX_PATH_LENGTH)) return 1;
    copy_string(filename, coeff_file, MAX_PATH_LENGTH + 1);

    /* call the Fortran subroutine */
    eval_qdlatlon_(qdlat, qdlon, time_dy, gcrad, gclat, gclon, &n_data, filename);

    return 0;
}


void eval_qdlatlonvb_(double*, double*, double*, double*, double*, double*,
                      double*, const double*, const double*, const double*,
                      const double*, const int*, const char*);

int c_eval_qdlatlonvb(
     double *qdlat, double *qdlon, double *f11, double *f12, double *f21,
     double *f22, double *f, const double *time_dy, const double *gcrad,
     const double *gclat, const double *gclon,
     const int n_data, const char *coeff_file)
{
    /* NOTE: Fortran code expects the file name as a sized size string. */
    char filename[MAX_PATH_LENGTH + 1];
    if (check_string_lenght(coeff_file, MAX_PATH_LENGTH)) return 1;
    copy_string(filename, coeff_file, MAX_PATH_LENGTH + 1);

    /* call the Fortran subroutine */
    eval_qdlatlonvb_(qdlat, qdlon, f11, f12, f21, f22, f,
                     time_dy, gcrad, gclat, gclon, &n_data, filename);

    return 0;
}
