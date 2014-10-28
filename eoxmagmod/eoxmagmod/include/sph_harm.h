/**
 * @file sph_harm.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Spherical Harmonics
 *
 * Various utilities needed by the spherical harmonic model evaluation.
 *
 * Copyright (C) 2014 EOX IT Services GmbH
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

#ifndef SPH_HARM_H
#define SPH_HARM_H 1

#include <math.h>

/**
 * @brief Evaluate the series of relative radius powers.
 *
 * Evaluate the series of relative radius powers for given relative radius ratio
 * 'relrad' (rad/rad0)
 *
 * (rad0/rad)**(i+2) for i = 0..degree
 *
 */

static void rel_rad_pow(double *rrp, int degree, double relrad)
{
    int i;
    double rr = 1.0/relrad;
    double rrp_last = rr;

    for (i = 0; i <= degree; ++i)
    {
        rrp_last *= rr;
        rrp[i] = rrp_last;
    }
}

/**
 * @brief Evaluate the series of azimut angle (lonigitude) sines and cosines.
 *
 * Evaluate the series of azimut angle (lonigitude - 'lon') sines and cosines.
 *
 *   cos(i*lon) for i = 0...degree
 *   sin(i*lon) for i = 0...degree
 *
 * This subroutine uses a faster evaluation based on pure recurrent
 * addition/substration and multiplication:
 *
 *   sin(i*lon) = cos((i-1)*lon)*sin(lon) + sin((i-1)*lon)*cos(lon)
 *   cos(i*lon) = cos((i-1)*lon)*cos(lon) - sin((i-1)*lon)*sin(lon)
 *
 * The input angle must be in radians.
 */

static void azmsincos(double *lonsin, double *loncos, int degree, double lon)
{
    int i;
    const double sin_lon = sin(lon);
    const double cos_lon = cos(lon);
    double sl, sl_new, cl, cl_new;

    lonsin[0] = 0.0;
    loncos[0] = 1.0;
    lonsin[1] = sl = sin_lon;
    loncos[1] = cl = cos_lon;

    for (i = 2; i <= degree; ++i)
    {
        sl_new = cl*sin_lon + sl*cos_lon;
        cl_new = cl*cos_lon - sl*sin_lon;
        lonsin[i] = sl = sl_new;
        loncos[i] = cl = cl_new;
    }
}

/**
 * @brief Evaluate the series of azimut angle (lonigitude) sines and cosines.
 *
 * Evaluate the series of azimut angle (lonigitude - 'lon') sines and cosines.
 *
 *   cos(i*lon) for i = 0...degree
 *   sin(i*lon) for i = 0...degree
 *
 * This subroutine contains the reference (slow) implementation evaluation
 * sine and cosine functions for each term of the series.
 *
 * The input angle must be in radians.
 */

static void azmsincos_ref(double *lonsin, double *loncos, int degree, double lon)
{
    int i;
    for (i = 0; i <= degree; ++i)
    {
        const double ilon = i*lon;
        lonsin[i] = sin(ilon);
        loncos[i] = cos(ilon);
    }
}

/**
 * @brief Evaluate the scalar potential and its (sherica) gradient/
 *
 * Sperical harmonic evaluation of the scalar potential and
 * the gradient in the spherical coordinates.
 *
 *  outputs:
 *    vpot - value of the scalar potential
 *    dvel - elevation component of the gradient
 *    dvaz - azimut component of the gradient
 *    dvrd - radial component of the gradient
 *
 *  inputs:
 *    degree - degree of the model
 *    mode - 3 - evaluate both gradient and potential
 *           1 - evaluate both potential only
 *           2 - evaluate both gradient only
 *    elv - elevation angle in radians (coordinate - needed in mode 3 and 2)
 *    rad - radius (coordinate needed in mode 3 and 1)
 *    cg, ch - spherical harmonic coeficients [(degree+1)*(degree+2)/2]
 *    lp, ldp - Legendre associative function and their derivtives (with
 *              respect to the elevation coordinate) [(degree+1)*(degree+2)/2]
 *    rrp - relative radius power series [degree+1]
 *    lsin, lcos - series of azimut angle sines and cosines [degree+1]
 *
 */

static void sph_harm_eval(double *vpot, double *dvel, double *dvaz, double *dvrd,
    int degree, int mode, double elv, double rad,
    const double *cg, const double *ch, const double *lp, const double *ldp,
    const double *rrp, const double *lsin, const double *lcos)
// the evaluation
{
    int i, j;
    const double sin_elv = sin(elv);
    const double cos_elv = cos(elv);
    double _vpot = 0.0, _dvel = 0.0, _dvaz = 0.0, _dvrd = 0.0;

    for (i = 1; i <= degree; ++i)
    {
        const int i_off = (i*(i+1))/2;

        for (j = 0; j <= i; ++j)
        {
            const int idx = i_off + j;
            const double tmp0 = (cg[idx]*lcos[j] + ch[idx]*lsin[j])*rrp[i];
            const double tmp1 = (cg[idx]*lsin[j] - ch[idx]*lcos[j])*rrp[i];

            _vpot += tmp0 * lp[idx];
            _dvel -= tmp0 * ldp[idx];
            _dvaz += tmp1 * lp[idx] * j;
            _dvrd -= tmp0 * lp[idx] * (i+1);

        }
    }

    if (mode&0x1)
    {
        *vpot = _vpot * rad;
    }

    if (mode&0x2)
    {
        *dvel = _dvel;
        *dvaz = _dvaz / cos_elv;
        *dvrd = _dvrd;

        // handling of the poles
        if (fabs(cos_elv) < 1e-10)
        {
            const double lsin1 = lsin[1], lcos1 = lcos[1];
            double sqn3, sqn1 = 1.0;
            double ps2, ps1 = 1.0, ps0 = 1.0,

            // i = 1
            _dvaz = (cg[2]*lsin1 - ch[2]*lcos1) * rrp[1];

            for (i = 2; i <= degree; ++i)
            {
                const int idx = 1 + (i*(i+1))/2;
                #define FDIV(a,b) ((double)(a)/(double)(b))
                const double tmp = FDIV((i-1)*(i-1)-1, (2*i-1)*(2*i-3));

                // evaluate ratio between the Gauss-normalised and Smidth
                // quasi-normalised associated Legendre functions.
                //  Equivalent to: sqrt((j==0?1:2)*(i-j)!/(i+j!))*(2i-1)!!/(i-j)!
                sqn1 = sqn1 * FDIV(2*i-1, i);
                sqn3 = sqn1 * sqrt(FDIV(i*2, i+1));
                #undef FDIV
                ps2 = ps1;
                ps1 = ps0;
                ps0 = sin_elv*ps1 - tmp*ps2;

                _dvaz += (cg[idx]*lsin1 - ch[idx]*lcos1) * rrp[i] * ps0 * sqn3;
            }

            *dvaz = _dvaz;
        }
    }
}

#endif  /*SPH_HARM_H*/
