/*-----------------------------------------------------------------------------
 *
 * Quasi-Dipole Magnetic Coordinates - C wrapper around the Fortran code
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

#ifndef CQDIPOLE_H
#define CQDIPOLE_H

/**
 * @brief Batch evaluation of the Quasi-Dipole apex coordinates
 *
 */

void c_make_apex(
    double* qdlat, double* qdlon, double* xmlt,
    double* f11, double* f12, double* f21, double* f22,
    const double* time, const double* gcrad, const double* gclat,
    const double* gclon, const int n_data, const char *fname);


/**
 * @brief Batch evaluation of the sub-solar point latitude and longitude
 *
 */

void c_eval_subsol(double* sbsllat, double* sbsllon, const double* time_mjd2k,
                   const int n_data);

#endif /*CQDIPOLE_H*/
