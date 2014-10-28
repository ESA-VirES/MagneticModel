/**
 * @file geo_conv.h
 * @author Martin Paces <martin.paces@eox.at>
 * @brief Geo-coordinates conversions.
 *
 * This file contains various geo-coordinates conversion utilities.
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

#ifndef GEO_CONV_H
#define GEO_CONV_H 1

#include <math.h>
#include "math_aux.h"

#define DG2RAD (M_PI/180.0)
#define RAD2DG (180.0*M_1_PI)

#define WGS84_A 6378.137
//#define WGS84_B 6356.752314245
#define WGS84_B 6356.7523142
#define WGS84_INVFLAT 298.257223563
#define WGS84_FLAT (1.0/WGS84_INVFLAT)
#define WGS84_EPS2 (1.0 - (WGS84_B*WGS84_B)/(WGS84_A*WGS84_A))
#define WGS84_EPS sqrt(WGS84_EPS2)
#define WGS84_RADIUS 6371.2

/**
 * @brief Convert geodetic to geocentric cartesian coordinates.
 *
 * Covert geodetic coordinates (latitude, longitude, elevation(above ellipsoid)
 * to geocentric cartesian (ECEF) coordinates.
 *
 * The input geodetic coordinates shall be in degrees. The units of the output
 * cartesian coordinates are the same as the one of the ellipsid semi-major axis.
 *
 * Required ellipsoid parameters:
 *  elp_a - semi-major axis
 *  elp_e2 - squared first excentricity
 */

static void geodetic2geocentric_cart(
    double *x, double *y, double *z,
    double lat, double lon, double elv,
    double elp_a, double elp_e2)
{
    double lat_rad = DG2RAD*lat;
    double lon_rad = DG2RAD*lon;
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);

    /* radius of curvature of the ellipsoid */
    double rc = elp_a / sqrt(1.0 - elp_e2*sin_lat*sin_lat);

    /* radius on the equator plane */
    double re = (rc + elv)*cos_lat;

    *x = re*cos(lon_rad);
    *y = re*sin(lon_rad);
    *z = (rc*(1.0 - elp_e2) + elv)*sin_lat;
}


/**
 * @brief Convert geocentric coordinates to geocentric sherical coordinates
 *
 * Covert geodetic coordinates (latitude, longitude, elevation(above ellipsoid)
 * to geocentric sperical coordinates (radius, elevation/theta, azimut/phi).
 *
 * The input geodetic coordinates shall be in degrees. The units of the output
 * cartesian coordinates as well as the unit of the radial coordinate are
 * the same as the one of the ellipsid semi-major axis. The output sperical
 * coordinates are in radians.
 *
 * Required ellipsoid parameters:
 *  elp_a - semi-major axis
 *  elp_e2 - squared first excentricity
 */

static void geodetic2geocentric_sph(
    double *r, double *th, double *ph,
    double lat, double lon, double elv,
    double elp_a, double elp_e2)
{
    double lat_rad = DG2RAD*lat;
    double lon_rad = DG2RAD*lon;
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);

    /* radius of curvature of the ellipsoid */
    double rc = elp_a / sqrt(1.0 - elp_e2*sin_lat*sin_lat);

    /* radius on the equator plane and the cartesian height*/
    double re = (rc + elv)*cos_lat;
    double  z = (rc*(1.0 - elp_e2) + elv)*sin_lat;

    *r = norm2d(re, z);
    *th = asin(z/(*r));
    *ph = lon_rad;
}


/**
 * @brief Convert geodetic coordinates to geocentric cartesian and shperical ones.
 *
 * Covert geodetic coordinates (latitude, longitude, elevation(above ellipsoid))
 * to geocentric coordinates both cartesian (x, y, z) and sperical (radius,
 * elevation/theta, azimut/phi). The input geodetic coordinates shall be in
 * degrees. The output sperical coordinates are in radians.
 * The input geodetic coordinates shall be in degrees. The unit of the output
 * the radial coordinate is  the same as the one of the ellipsid semi-major
 * axis. The output sperical coordinates are in radians.
 *
 * Use when both, shperical and cartesian, coordinates are needed.
 * Note that this function is more efficent that two calls to
 *  geodetic2geocentric_cart
 *  geodetic2geocentric_sph
 *
 * Required ellipsoid parameters:
 *  elp_a - semi-major axis
 *  elp_e2 - squared first excentricity
 */

static void geodetic2geocentric(
    double *x, double *y, double *z,
    double *r, double *th, double *ph,
    double lat, double lon, double elv,
    double elp_a, double elp_e2)
{
    double lat_rad = DG2RAD*lat;
    double lon_rad = DG2RAD*lon;
    double sin_lat = sin(lat_rad);
    double cos_lat = cos(lat_rad);

    /* radius of curvature of the ellipsoid */
    double rc = elp_a / sqrt(1.0 - elp_e2*sin_lat*sin_lat);

    /* radius on the equator plane */
    double re = (rc + elv)*cos_lat;

    *x = re*cos(lon_rad);
    *y = re*sin(lon_rad);
    *z = (rc*(1.0 - elp_e2) + elv)*sin_lat;

    *r = norm2d(re, *z);
    *th = asin(*z/(*r));
    *ph = lon_rad;
}

/**
 * @brief Convert spherical to cartesian coordinates
 *
 * Covert spherical coordinates(radius, elevation/theta, azimut/phi)
 * to cartesian coordinates (x, y, z).
 *
 * The spherical coordinates shall be in radians.
 */

static void sph2cart(
    double *x, double *y, double *z,
    double r, double th, double ph)
{
    double sin_th = sin(th);
    double cos_th = cos(th);
    double sin_ph = sin(ph);
    double cos_ph = cos(ph);

    /* radius on the azimut plane (z=0)*/
    double ra = r*cos_th;

    *x = ra*cos_ph;
    *y = ra*sin_ph;
    *z = r*sin_th;
}

/**
 * @brief Convert cartesian to spherical coordinates
 *
 * Covert cartesian coordinates (x, y, z) to spherical
 * coordinates (radius, elevation/theta, azimut/phi)
 *
 * The spherical coordinates are produced in radians.
 */
static void cart2sph(
    double *r, double *th, double *ph,
    double x, double y, double z)
{
    *r = norm3d(x, y, z);
    *th = asin(z/(*r));
    *ph = (y >= 0 ? M_PI_2 : -M_PI_2) - atan(x/y);
}


/**
 * @brief Convert cartesian geocentric to geodetic coordinates
 *
 * Covert cartesian (ECEF) coordinates (x, y, z) to geodetic coordinates
 * (latitude, longitude, elevation(above ellipsoid).
 *
 * The geodetic coordinates are produced in degrees. The unit of the elevation
 * is the same as the unit of the cartesian coordinates and of the ellipsoid
 * semi-major axis.
 */

static void geocentric_cart2geodetic(
    double *lat, double *lon, double *elv,
    double x, double y, double z,
    double elp_a, double elp_e2)
{
    double p = norm2d(x, y);
    /* Ferraris's solution */
    double pa = p/elp_a;
    double za = z/elp_a;
    double pa2 = pa*pa;
    double za2 = za*za;
    double ee4 = elp_e2*elp_e2;
    double rkp0 = (1.0 - elp_e2);
    double zt = rkp0*za2;
    double rh = (pa2 + zt - ee4)/6.0;
    double ss = 0.25*zt*ee4*pa2;
    double rh3 = rh*rh*rh;
    double tmp = rh3 + ss + sqrt(ss*(ss+2.0*rh3));
    double tt = copysign(pow(fabs(tmp), 1.0/3.0), tmp);
    double uu = rh + tt + rh*rh/tt;
    double vv = sqrt(uu*uu + ee4*zt);
    double ww = 0.5*elp_e2*(uu + vv - zt)/vv;
    double kp = 1.0 + elp_e2*(sqrt(uu + vv + ww*ww) + ww)/(uu + vv);

    *lat = RAD2DG*atan(kp*z/p);
    *lon = (y >= 0 ? +90.0 : -90.0) - RAD2DG*atan(x/y);
    *elv = norm2d(p, z*kp)*(1.0/kp - rkp0)/elp_e2;
}

/**
 * @brief Convert spherical geocentric to geodetic coordinates
 *
 * Convert geocentric spherical coordinates (radius, elevation/theta, azimut/phi)
 * to geodetic coordinates (latitude, longitude, elevation(above ellipsoid).
 *
 * The spherical coordinates are expect to be in radians.
 * The geodetic coordinates are produced in degrees. The unit of the elevation
 * is the same as the unit of the radial coordinate and of the ellipsoid
 * semi-major axis.
 */

static void geocentric_sph2geodetic(
    double *lat, double *lon, double *elv,
    double r, double th, double ph,
    double elp_a, double elp_e2)
{
    double z = r*sin(th);
    double p = r*cos(th);
    /* Ferrari's solution */
    double pa = p/elp_a;
    double za = z/elp_a;
    double pa2 = pa*pa;
    double za2 = za*za;
    double ee4 = elp_e2*elp_e2;
    double rkp0 = (1.0 - elp_e2);
    double zt = rkp0*za2;
    double rh = (pa2 + zt - ee4)/6.0;
    double ss = 0.25*zt*ee4*pa2;
    double rh3 = rh*rh*rh;
    double tmp = rh3 + ss + sqrt(ss*(ss+2.0*rh3));
    double tt = copysign(pow(fabs(tmp), 1.0/3.0), tmp);
    double uu = rh + tt + rh*rh/tt;
    double vv = sqrt(uu*uu + ee4*zt);
    double ww = 0.5*elp_e2*(uu + vv - zt)/vv;
    double kp = 1.0 + elp_e2*(sqrt(uu + vv + ww*ww) + ww)/(uu + vv);

    *lat = RAD2DG*atan(kp*z/p);
    *lon = RAD2DG*ph;
    *elv = norm2d(p, z*kp)*(1.0/kp - rkp0)/elp_e2;
}

#endif  /*GEO_CONV_H*/

