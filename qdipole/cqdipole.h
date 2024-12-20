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

#ifndef CQDIPOLE_H
#define CQDIPOLE_H

/**
 * @brief Get filename size limit.
 * @returns maximum length of the filename path
 */

size_t get_qdipole_max_fname_lenght();


/**
 * @brief Get library version.
 * @returns String containing version of the library
 */

const char* get_qdipole_version();


/**
 * @brief Batch evaluation of the Quasi-Dipole apex coordinates with the vector base
 * @returns 0 on success, non-zero on failure
 */

int c_eval_qdlatlonvb(
    double* qdlat, double* qdlon, double* f11, double* f12, double* f21,
    double* f22, double* f, const double* time_dy, const double* gcrad,
    const double* gclat, const double* gclon, const int n_data,
    const char* coeff_file
);


/**
 * @brief Batch evaluation of the Quasi-Dipole apex coordinates
 * @returns 0 on success, non-zero on failure
 */

int c_eval_qdlatlon(
    double* qdlat, double* qdlon, const double* time_dy, const double* gcrad,
    const double* gclat, const double* gclon, const int n_data,
    const char* coeff_file
);


/**
 * @brief Batch evaluation of the Magnetic Local Time
 * @returns 0 on success, non-zero on failure
 */

int c_eval_mlt(
    double *t_mlt, const double *qdlon, const double *t_mjd2k, const int n_data
);


/**
 * @brief Batch evaluation of the sub-solar point latitude and longitude
 * @returns 0 on success, non-zero on failure
 */

int c_eval_subsol(
    double* sbsllat, double* sbsllon, const double* time_mjd2k, const int n_data
);

#endif /*CQDIPOLE_H*/
