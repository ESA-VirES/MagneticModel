/**
 * @file spherical_harmonics.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Spherical Harmonics
 *
 * Various utilities needed by the spherical harmonic model evaluation.
 *
 * Copyright (C) 2014-2022 EOX IT Services GmbH
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

#ifndef SPHERICAL_HARMONICS_H
#define SPHERICAL_HARMONICS_H 1

#define SHC_LEG_HIGH_SCALE 1.0e+280
#define FDIV(a,b) ((double)(a)/(double)(b))

#include <math.h>

/**
 * @brief Allocate and evaluate pre-calculated square-root series
 *
 *   sqrt(i) for i = 0..(2*degree+1)
 *
 *   NOTE: The returned pointer *MUST* be passed to free().
 */

static double * shc_presqrt(int degree)
{
    const int nval = 2*degree + 2;
    int i;
    double * psqrt = malloc(nval*sizeof(double));

    if (NULL == psqrt)
        return NULL;

    for (i = 0; i < nval ; ++i)
    {
        psqrt[i] = sqrt((double) i);
    }

    return psqrt;
}


/**
 * @brief Allocate and evaluate pre-calculated reciprocal square-root series
 *
 *   sqrt(i) for i = 0..(2*degree+1)
 *
 *   NOTE: The returned pointer *MUST* be passed to free().
 */

static double * shc_prersqrt(int degree)
{
    const int nval = 2*degree + 2;
    int i;
    double * prsqrt = malloc(nval*sizeof(double));

    if (NULL == prsqrt)
        return NULL;

    for (i = 0; i < nval ; ++i)
    {
        prsqrt[i] = 1.0 / sqrt((double) i);
    }

    return prsqrt;
}


/**
 * @brief Evaluate the associated Legendre functions
 *
 * Evaluate the Schmidt semi-normalised associated Legendre functions.
 *
 * This version is susceptible to overflow for higher orders (degree > 20)
 * especially close to the poles.
 */

static void shc_legendre_low(double *lp, double *ldp, int degree, double sin_elv, double cos_elv, const double *psqrt, const double *prsqrt)
{
    int i, j;

    // Gauss normalised associated Legendre functions (aLf)
    // Note the ldp sign change because of the derivative being calculated
    // with respect to the elevation (latitude) instead of the
    // inclination (co-latitude).
    lp[0] = 1.0;
    lp[1] = sin_elv ;
    lp[2] = cos_elv;

    ldp[0] = 0.0 ;
    ldp[1] = cos_elv;
    ldp[2] = -sin_elv;

    for (i = 2; i <= degree; ++i)
    {
        const int i_off = (i*(i+1))/2;
        const int i1sq = (i-1)*(i-1);
        const double w_dnm = 1.0/((2*i-1)*(2*i-3));

        for (j = 0; j < (i-1); ++j)
        {
            const int idx0 = i_off + j ;
            const int idx1 = i_off + j - i;
            const int idx2 = i_off + j - i - i + 1;
            const double w = w_dnm*(i1sq - j*j);

            lp [idx0] = sin_elv*lp [idx1] - w*lp [idx2];
            ldp[idx0] = sin_elv*ldp[idx1] - w*ldp[idx2] + cos_elv*lp[idx1];
        }

        {
            const int idx0 = i_off + i;
            const int idx1 = i_off - 1;

            lp [idx0-1] = sin_elv*lp [idx1];
            ldp[idx0-1] = sin_elv*ldp[idx1] + cos_elv*lp[idx1];

            lp [idx0] = cos_elv*lp [idx1];
            ldp[idx0] = cos_elv*ldp[idx1] - sin_elv*lp[idx1];
        }
    }

    // Convert Gauss normalised aLf to Schmidt quasi-normalised aLf
    {
        double sqn1, sqn0 = 1.0 ;

        for (i = 2; i <= degree; ++i)
        {
            const int i_off = (i*(i+1))/2;

            sqn0 *= FDIV(2*i-1,i);
            lp[i_off] *= sqn0;
            ldp[i_off] *= sqn0;

            sqn1 = sqn0*psqrt[2*i]*prsqrt[i+1];

            lp[i_off+1] *= sqn1;
            ldp[i_off+1] *= sqn1;

            for (j = 2; j <= i; ++j)
            {
                sqn1 *= psqrt[i-j+1]*prsqrt[i+j];
                lp[i_off+j] *= sqn1;
                ldp[i_off+j] *= sqn1;
            }
        }
    }
}


/**
 * @brief Evaluate the associated Legendre functions -
 *
 * Evaluate the Schmidt semi-normalised associated Legendre functions.
 *
 * This subroutine tries to compensate the effect of underflow for high degree
 * values near the poles by scaling of the of the powers of the .
 * [Holmes and Featherstone 2002, J. Geodesy, 76, 279-299]
 */


static void shc_legendre_high(double *lp, double *ldp, int degree, double sin_elv, double cos_elv, const double *psqrt, const double *prsqrt)
{
    int i, j, idx, idx0;
    const double scale = SHC_LEG_HIGH_SCALE;
    const double rscale = 1.0/SHC_LEG_HIGH_SCALE;
    const double rcos_elv = 1.0/cos_elv ;
    double pm2, pm1, pm0, pmm;
    double pow_cos_elv;

    // Gauss normalised associated Legendre functions (aLf)
    lp[0] = pm2 = 1.0;
    lp[1] = pm1 = sin_elv;

    ldp[0] = 0.0 ;
    ldp[1] = cos_elv;

    for (i = 2, idx = 3; i <= degree; idx += ++i)
    {
        const double f1 = FDIV(2*i-1,i);
        const double f2 = FDIV(i-1,i);

        lp[idx] = pm0 = f1*sin_elv*pm1 - f2*pm2;
        ldp[idx] = i*rcos_elv*(pm1 - sin_elv*pm0);

        pm2 = pm1;
        pm1 = pm0;
    }

    pmm = rscale*psqrt[2];
    pow_cos_elv = scale;

    for (j = 1, idx0 = 2; j < degree; idx0 += ++j + 1)
    {
        const double pow_cos_elv_last = pow_cos_elv;
        pow_cos_elv *= cos_elv;

        // i = j
        pmm *= psqrt[2*j+1]*prsqrt[2*j];
        pm2 = pmm*prsqrt[2*j+1];
        pm1 = sin_elv*pmm;

        lp[idx0] = pow_cos_elv*pm2;
        ldp[idx0] = -pow_cos_elv_last*(j*sin_elv*pm2);

        // i = j + 1
        idx = idx0 + j + 1;
        lp[idx] = pow_cos_elv*pm1;
        ldp[idx] = -pow_cos_elv_last*((j+1)*sin_elv*pm1 - pmm);

        for(i = j + 2, idx += j + 2; i <= degree; idx += ++i)
        {
            const double f0 = psqrt[i+j]*psqrt[i-j];
            const double rf0 = prsqrt[i+j]*prsqrt[i-j];
            const double f1 = (2*i-1);
            const double f2 = psqrt[i-j-1]*psqrt[i+j-1];
            const double plm = (f1*sin_elv*pm1 - f2*pm2)*rf0;

            lp[idx] = pow_cos_elv*plm;
            ldp[idx] = -pow_cos_elv_last*(i*sin_elv*plm - f0*pm1);

            pm2 = pm1;
            pm1 = plm;
        }
    }

    // i = j = degree
    pm2 = pmm*prsqrt[2*j];

    lp[idx0] = pow_cos_elv*cos_elv*pm2;
    ldp[idx0] = -pow_cos_elv*(j*sin_elv*pm2);
}


/**
 * @brief Evaluate the associated Legendre functions
 *
 * Evaluate the Schmidt semi-normalised associated Legendre functions and
 * their derivatives with respect to the elevation (latitude) coordinate.
 *
 * parameters:
 *
 *      lp ... output array of the associated Legendre functions
 *      ldp ... output array of the associated Legendre functions' derivatives
 *      degree ... model degree
 *      elv ... elevation (latitude) angle in radians
 *      psqrt ... array of the pre-calculated square roots
 */

static void shc_legendre(double *lp, double *ldp, int degree, double elv, const double *psqrt, const double *prsqrt)
{
    const double sin_elv = sin(elv);
    //const double cos_elv = cos(elv);
    const double cos_elv = sqrt((1.0-sin_elv)*(1.0+sin_elv));

    if ((degree <= 16)||(fabs(cos_elv) < 1e-10))
        // low degree model + poles
        shc_legendre_low(lp, ldp, degree, sin_elv, cos_elv, psqrt, prsqrt);
    else
        // high degree model
        shc_legendre_high(lp, ldp, degree, sin_elv, cos_elv, psqrt, prsqrt);
}


/**
 * @brief Evaluate the "internal" series of relative radius powers.
 *
 * Evaluate the series of relative radius powers for given relative radius ratio
 * 'relrad' (rad/rad0)
 *
 * (rad0/rad)**(i+2) for i = 0..degree
 *
 */

static void shc_relradpow_internal(double *rrp, int degree, double relrad)
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
 * @brief Evaluate the "external" series of relative radius powers.
 *
 * Evaluate the series of relative radius powers for given relative radius ratio
 * 'relrad' (rad/rad0)
 *
 * (rad/rad0)**(i-1) for i = 0..degree
 *
 */

static void shc_relradpow_external(double *rrp, int degree, double relrad)
{
    int i;
    double rrp_val = 1.0;

    if (0 <= degree)
        rrp[0] = rrp_val/relrad;

    if (1 <= degree)
        rrp[1] = rrp_val;

    for (i = 2; i <= degree; ++i)
    {
        rrp_val *= relrad;
        rrp[i] = rrp_val;
    }
}


/**
 * @brief Evaluate the series of azimuth angle (longitude) sines and cosines.
 *
 * Evaluate the series of azimuth angle (longitude - 'lon') sines and cosines.
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

static void shc_azmsincos(double *lonsin, double *loncos, int degree, double lon)
{
    int i;
    const double sin_lon = sin(lon);
    const double cos_lon = cos(lon);
    double sl, sl_new, cl, cl_new;

    lonsin[0] = 0.0;
    loncos[0] = 1.0;

    if (degree > 0)
    {
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
}


/**
 * @brief Evaluate the series of azimuth angle (longitude) sines and cosines.
 *
 * Evaluate the series of azimuth angle (longitude - 'lon') sines and cosines.
 *
 *   cos(i*lon) for i = 0...degree
 *   sin(i*lon) for i = 0...degree
 *
 * This subroutine contains the reference (slow) implementation evaluation
 * sine and cosine functions for each term of the series.
 *
 * The input angle must be in radians.
 */

static void shc_azmsincos_ref(double *lonsin, double *loncos, int degree, double lon)
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
 * @brief Evaluate the scalar potential
 *
 * Spherical harmonic evaluation of the scalar potential
 *
 *  outputs:
 *    vpot - value of the scalar potential
 *
 *  inputs:
 *    degree      - degree of the model
 *    rad         - radius
 *    cg, ch      - spherical harmonic coefficients [(degree+1)*(degree+2)/2]
 *    lp          - Legendre associative function [(degree+1)*(degree+2)/2]
 *    rrp         - relative radius power series [degree+1]
 *    lsin, lcos  - series of azimuth angle sines and cosines [degree+1]
 */

static void shc_eval_v(
    double *vpot, int degree, double rad,
    const double *cg, const double *ch, const double *lp,
    const double *rrp, const double *lsin, const double *lcos
)
{
    int i, j;
    double _vpot;

    // i = 0
    _vpot = cg[0]*rrp[0];

    for (i = 1; i <= degree; ++i)
    {
        double _vpot_i;
        const int i_off = (i*(i+1))/2;

        _vpot_i = cg[i_off] * lp[i_off];

        for (j = 1; j <= i; ++j)
        {
            const int idx = i_off + j;

            _vpot_i += (cg[idx]*lcos[j] + ch[idx]*lsin[j]) * lp[i_off+j];
        }

        _vpot += rrp[i] * _vpot_i;
    }

    *vpot = _vpot * rad;
}


/**
 * @brief Evaluate the azimuth component of the gradient of potential near poles
 *
 * Spherical harmonic evaluation of the azimuth component of the gradient
 * in the vicinity of the poles.
 *
 *  outputs:
 *    dvaz - azimuth component of the gradient
 *
 *  inputs:
 *    degree - degree of the model
 *    elv - elevation angle in radians
 *    cg, ch - spherical harmonic coefficients [(degree+1)*(degree+2)/2]
 *    rrp - relative radius power series [degree+1]
 *    lsin, lcos - series of azimuth angle sines and cosines [degree+1]
 */

static void shc_eval_dvaz_pole(
    double *dvaz,
    int degree, double elv,
    const double *cg, const double *ch,
    const double *rrp, const double *lsin, const double *lcos
)
{
    int i;
    const double sin_elv = sin(elv);

    const double lsin1 = lsin[1], lcos1 = lcos[1];
    double sqn3, sqn1 = 1.0;
    double ps2, ps1 = 1.0, ps0 = 1.0,

    // i = 1
    _dvaz = (cg[2]*lsin1 - ch[2]*lcos1) * rrp[1];

    for (i = 2; i <= degree; ++i)
    {
        const int idx = 1 + (i*(i+1))/2;
        const double tmp = FDIV((i-1)*(i-1)-1, (2*i-1)*(2*i-3));

        // evaluate ratio between the Gauss-normalised and Schmidt
        // quasi-normalised associated Legendre functions.
        //  Equivalent to: sqrt((j==0?1:2)*(i-j)!/(i+j!))*(2i-1)!!/(i-j)!
        sqn1 = sqn1 * FDIV(2*i-1, i);
        sqn3 = sqn1 * sqrt(FDIV(i*2, i+1));
        ps2 = ps1;
        ps1 = ps0;
        ps0 = sin_elv*ps1 - tmp*ps2;

        _dvaz += (cg[idx]*lsin1 - ch[idx]*lcos1) * rrp[i] * ps0 * sqn3;
    }

    *dvaz = -_dvaz;
}


/**
 * @brief Evaluate the gradient of potential
 *
 * Spherical harmonic evaluation of the gradient in the spherical coordinates.
 *
 *  outputs:
 *    dvel - elevation component of the gradient
 *    dvaz - azimuth component of the gradient
 *    dvrd - radial component of the gradient
 *
 *  inputs:
 *    degree - degree of the model
 *    elv - elevation angle in radians
 *    cg, ch - spherical harmonic coefficients [(degree+1)*(degree+2)/2]
 *    lp, ldp - Legendre associative function and their derivatives (with
 *              respect to the elevation coordinate) [(degree+1)*(degree+2)/2]
 *    rrp - relative radius power series [degree+1]
 *    lsin, lcos - series of azimuth angle sines and cosines [degree+1]
 *    is_internal - boolean flag indicating type of the evaluated field
                    set true for the internal or false for the external field.
 */

static void shc_eval_dv(
    double *dvel, double *dvaz, double *dvrd,
    int degree, double elv,
    const double *cg, const double *ch, const double *lp, const double *ldp,
    const double *rrp, const double *lsin, const double *lcos, int is_internal
)
{
    int i, j;
    const double cos_elv = cos(elv);
    double _dvrd, _dvel, _dvaz;
    const int dr_scale = is_internal ? 1 : -1;
    const int dr_offset = is_internal ? 1 : 0;

    { // i = 0
        const double tmp = cg[0]*rrp[0];
        const int i_dr = dr_scale*dr_offset;

        _dvel = 0.0;          // north-ward
        _dvaz = 0.0;          // east-ward
        _dvrd = tmp * -i_dr ; // up-ward
    }

    for (i = 1; i <= degree; ++i)
    {
        double _dvel_i, _dvaz_i, _dvrd_i;
        const int i_off = (i*(i+1))/2;
        const int i_dr = dr_scale*(i + dr_offset);

        _dvel_i = cg[i_off] * ldp[i_off];        // north-ward
        _dvaz_i = 0.0;                           // east-ward
        _dvrd_i = cg[i_off] * lp[i_off] * i_dr ; // up-ward

        for (j = 1; j <= i; ++j)
        {
            const int idx = i_off + j;
            const double tmp0 = cg[idx]*lcos[j] + ch[idx]*lsin[j];
            const double tmp1 = cg[idx]*lsin[j] - ch[idx]*lcos[j];

            _dvel_i += tmp0 * ldp[idx];         // north-ward
            _dvaz_i += tmp1 * lp[idx] * j;      // east-ward
            _dvrd_i += tmp0 * lp[idx] * i_dr ;  // up-ward
        }

        _dvel += rrp[i] * _dvel_i;
        _dvaz -= rrp[i] * _dvaz_i;
        _dvrd -= rrp[i] * _dvrd_i;
    }

    *dvel = _dvel;
    *dvaz = _dvaz / cos_elv;
    *dvrd = _dvrd;

    // special handling of the azimuth component in the vicinity of the poles
    if ((degree > 0) && (fabs(cos_elv) < 1e-10))
    {
        shc_eval_dvaz_pole(dvaz, degree, elv, cg, ch, rrp, lsin, lcos);
    }
}


/**
 * @brief Evaluate the scalar potential and its (spherical) gradient/
 *
 * Spherical harmonic evaluation of the scalar potential and
 * the gradient in the spherical coordinates.
 *
 *  outputs:
 *    vpot - value of the scalar potential
 *    dvel - elevation component of the gradient
 *    dvaz - azimuth component of the gradient
 *    dvrd - radial component of the gradient
 *
 *  inputs:
 *    degree - degree of the model
 *    elv - elevation angle in radians
 *    rad - radius
 *    cg, ch - spherical harmonic coefficients [(degree+1)*(degree+2)/2]
 *    lp, ldp - Legendre associative function and their derivatives (with
 *              respect to the elevation coordinate) [(degree+1)*(degree+2)/2]
 *    rrp - relative radius power series [degree+1]
 *    lsin, lcos - series of azimuth angle sines and cosines [degree+1]
 *    is_internal - boolean flag indicating type of the evaluated field
                    set true for the internal or false for the external field.
 */

static void shc_eval_v_dv(
    double *vpot, double *dvel, double *dvaz, double *dvrd,
    int degree, double elv, double rad,
    const double *cg, const double *ch, const double *lp, const double *ldp,
    const double *rrp, const double *lsin, const double *lcos, int is_internal
)
{
    int i, j;
    const double cos_elv = cos(elv);
    double _vpot, _dvrd, _dvel, _dvaz;
    const int dr_scale = is_internal ? 1 : -1;
    const int dr_offset = is_internal ? 1 : 0;

    { // i = 0
        const double tmp = cg[0]*rrp[0];
        const int i_dr = dr_scale*dr_offset;

        _vpot = tmp;
        _dvel = 0.0;          // north-ward
        _dvaz = 0.0;          // east-ward
        _dvrd = tmp * -i_dr ; // up-ward
    }

    for (i = 1; i <= degree; ++i)
    {
        double _vpot_i, _dvel_i, _dvaz_i, _dvrd_i;
        const int i_off = (i*(i+1))/2;
        const int i_dr = dr_scale*(i + dr_offset);

        _vpot_i = cg[i_off] * lp[i_off];
        _dvel_i = cg[i_off] * ldp[i_off];        // north-ward
        _dvaz_i = 0.0;                           // east-ward
        _dvrd_i = cg[i_off] * lp[i_off] * i_dr ; // up-ward

        for (j = 1; j <= i; ++j)
        {
            const int idx = i_off + j;
            const double tmp0 = cg[idx]*lcos[j] + ch[idx]*lsin[j];
            const double tmp1 = cg[idx]*lsin[j] - ch[idx]*lcos[j];

            _vpot_i += tmp0 * lp[i_off+j];
            _dvel_i += tmp0 * ldp[idx];         // north-ward
            _dvaz_i += tmp1 * lp[idx] * j;      // east-ward
            _dvrd_i += tmp0 * lp[idx] * i_dr ;  // up-ward
        }

        _vpot += rrp[i] * _vpot_i;
        _dvel += rrp[i] * _dvel_i;
        _dvaz -= rrp[i] * _dvaz_i;
        _dvrd -= rrp[i] * _dvrd_i;
    }

    *vpot = _vpot * rad;
    *dvel = _dvel;
    *dvaz = _dvaz / cos_elv;
    *dvrd = _dvrd;

    // special handling of the azimuth component in the vicinity of the poles
    if ((degree > 0) && (fabs(cos_elv) < 1e-10))
    {
        shc_eval_dvaz_pole(dvaz, degree, elv, cg, ch, rrp, lsin, lcos);
    }
}


/**
 * @brief Evaluate the scalar potential and its (spherical) gradient/
 *
 * Spherical harmonic evaluation of the scalar potential and
 * the gradient in the spherical coordinates.
 *
 *  outputs:
 *    vpot - value of the scalar potential
 *    dvel - elevation component of the gradient
 *    dvaz - azimuth component of the gradient
 *    dvrd - radial component of the gradient
 *
 *  inputs:
 *    degree - degree of the model
 *    mode - 3 - evaluate both gradient and potential
 *           1 - evaluate potential only
 *           2 - evaluate gradient only
 *    elv - elevation angle in radians (coordinate - needed in mode 3 and 2)
 *    rad - radius (coordinate needed in mode 3 and 1)
 *    cg, ch - spherical harmonic coefficients [(degree+1)*(degree+2)/2]
 *    lp, ldp - Legendre associative function and their derivatives (with
 *              respect to the elevation coordinate) [(degree+1)*(degree+2)/2]
 *    rrp - relative radius power series [degree+1]
 *    lsin, lcos - series of azimuth angle sines and cosines [degree+1]
 *    is_internal - boolean flag indicating type of the evaluated field
                    set true for the internal or false for the external field.
 */

static void shc_eval(
    double *vpot, double *dvel, double *dvaz, double *dvrd,
    int degree, int mode, double elv, double rad,
    const double *cg, const double *ch, const double *lp, const double *ldp,
    const double *rrp, const double *lsin, const double *lcos, int is_internal
)
// the evaluation
{
    switch (mode & 0x3)
    {
        case 0x1:
            shc_eval_v(
                vpot, degree, rad,
                cg, ch, lp, rrp, lsin, lcos
            );
            break;
        case 0x2:
            shc_eval_dv(
                dvel, dvaz, dvrd, degree, elv,
                cg, ch, lp, ldp, rrp, lsin, lcos, is_internal
            );
            break;
        case 0x3:
            shc_eval_v_dv(
                vpot, dvel, dvaz, dvrd, degree, elv, rad,
                cg, ch, lp, ldp, rrp, lsin, lcos, is_internal
            );
            break;
    }
}

#undef FDIV
#endif  /*SPHERICAL_HARMONICS_H*/
