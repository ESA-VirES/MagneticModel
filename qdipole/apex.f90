!  FILE NAME: apex.f

      SUBROUTINE APEX (DATE,DLAT,DLON,ALT, &
                       A,ALAT,ALON,BMAG,XMAG,YMAG,ZMAG,V)
!          Calculate apex radius, latitude, longitude; and magnetic field and
!          scalar magnetic potential.
!
!          INPUTS:
!            DATE = Year and fraction (1990.0 = 1990 January 1, 0 UT)
!            DLAT = Geodetic latitude in degrees
!            DLON = Geodetic longitude in degrees
!            ALT = Altitude in km
!
!          RETURNS:
!            A    = (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
!                   A is analogous to the L value in invariant coordinates.
!            ALAT = Apex latitude in degrees (negative in S. magnetic hemisphere)
!            ALON = Apex longitude (geomagnetic longitude of apex) in degrees
!            BMAG = geomagnetic field magnitude (nT)
!            XMAG = geomagnetic field component (nT): north
!            YMAG = geomagnetic field component (nT): east
!            ZMAG = geomagnetic field component (nT): downward
!            V    = geomagnetic potential (T-m)
!
!          COMMON BLOCKS:
!            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
!
!          DIPOLE has IGRF variables obtained from routines in magfld.f:
!            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
!            ELON  = East longitude of geomagnetic dipole north pole (deg)
!            VP    = Magnitude (T-m) of dipole component of magnetic potential at
!                    geomagnetic pole and geocentric radius of 6371.0088 km
!            CTP   = cosine of COLAT
!            STP   = sine   of COLAT
!
!------------------------------------------------------------------------------
!          HISTORY:
!          Aug 1994: First version completed on the 22nd by A.D. Richmond.
!          May 1999: Revise DS calculation in LINAPX to avoid divide by zero.
!          Apr 2004: - Change definition of earth's equatorial radius (REQ)
!                      from the IAU-1966 spheroid (6378.160 km) to the WGS-1984
!                      spheroid (6378.137 km); see description below.
!                    - Revise comments toward a consistent format so they are
!                      easy to read.
!                    - Replace computed GO TO in ITRACE with IF blocks.
!                    - Refine FNDAPX to insure |Bdown/Btot| < 1.E-6 at apex
!          Nov 2009: Change definition of earth's mean radius (RE) from 6371.2
!                    to the WGS84 value (6371.0088), by J.T. Emmert, NRL
!
!------------------------------------------------------------------------------
!                Reference Spheroid Change                March 2004
!
!   Apex geomagnetic coordinates are based on the International Reference
!   Geomagnetic Field (IGRF) which involves the earth's shape when converting
!   geographic latitude and altitude to geocentric coordinates.  For this
!   purpose, the earth is assumed to be an ellipsoid which is fatter at the
!   equator than the poles.  The World Geodetic System 1984 spheroid
!   (WGS-1984) is recommended in the recent release of IGRF-9 because it is
!   used to position current satellite magnetic data (EOS Volume 84 Number 46
!   November 18 2003).  This differs from previous IGRF releases which favored
!   the International Astronomical Union 1966 spheroid (IAU-1966) so the Apex
!   program conversion from geographic to geocentric coordinates in subroutine
!   CONVRT of file magfld.f has been revised from the IAU-1966 spheroid to the
!   WGS-1984 spheroid.
!
!   The spheroid used to prepare earlier IGRF releases is not always known but
!   changing spheroids now produces differences at the earth's surface less
!   than 1 nT, significantly less than other error sources: viz., 9 nT RMS
!   error for older measurements reporting 1 nT accuracy, 200 nT in the
!   vicinity of magnetized rocks, or as much as 1000 nT during and after a
!   geomagnetic storm (www.ngdc.noaa.gov/IAGA/vmod/index.html).
!
!   The elliptical shape is characterized in subroutine CONVRT by eccentricity
!   (e) which is related to the the earth's equatorial radius (a) and polar
!   radius (b) by
!
!        e**2  = 1 - (b/a)**2     (1)
!
!   This term is part of an eighth order Lagrange expansion formula (Astron.
!   J.  Vol. 66, p. 15-16, 1961) designed to give eight digit conversion
!   accuracy.  The following table summarizes the relevant spheroids:
!
!          a           b           e**2          Source
!        ----------- -----------  -----------  --------------
!          -           -        0.006722670  Astron J. 1961
!        6378.160 km 6356.775 km  0.006701642  IAU-1966
!        6378.137 km 6356.752 km  0.006694478  WGS-1984
!
!   The previous formulation in CONVRT used the oblateness factor (epsilon),
!   a surrogate for eccentricity, where
!
!        e**2 = 2*epsilon - epsilon**2
!
!   with oblateness revised from the 1961 paper's suggested 1/297 to 1/298.25
!   in accordance with the IAU-1966 spheroid.  Now CONVRT is reformulated to
!   use equation 1 with the WGS-1984 spheroid's major axis (a) and minor axis
!   (b) for which epsilon is approximately 1/298.2528.
!
!   In addition to earth's equatorial radius and polar radius, the reference
!   radius (Re) of 6371.2 km is explicit in the IGRF formula for magnetic
!   potential and implicit in the derived magnetic field coefficients.  The
!   reference radius has not changed in IGRF releases.
!
!------------------------------------------------------------------------------

      PARAMETER (RE = 6371.0088, DTOR = .01745329251994330)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP

      CALL COFRM (DATE)
      CALL DYPOL (CLATP,POLON,VPOL)
      COLAT = CLATP
      CTP   = COS(CLATP*DTOR)
      STP   = SQRT(1. - CTP*CTP)
      ELON  = POLON
      VP    = VPOL
      CALL LINAPX (DLAT,DLON,ALT, A,ALAT,ALON,XMAG,YMAG,ZMAG,BMAG)
      XMAG = XMAG*1.E5                 ! convert from gauss to nT
      YMAG = YMAG*1.E5
      ZMAG = ZMAG*1.E5
      BMAG = BMAG*1.E5
      CALL GD2CART (DLAT,DLON,ALT,X,Y,Z)
      CALL FELDG (3, X/RE,Y/RE,Z/RE, BX,BY,BZ,V)
      RETURN
      END

      SUBROUTINE LINAPX (GDLAT,GLON,ALT, A,ALAT,ALON,XMAG,YMAG,ZMAG,F)

!          Transform geographic coordinates to Apex coordinates.
!
!          INPUTS:
!            GDLAT = Latitude  (degrees, positive northward)
!            GLON  = Longitude (degrees, positive eastward)
!            ALT   = Height of starting point (km above mean sea level)
!
!          OUTPUTS:
!            A     = (Apex height + REQ)/REQ, where REQ = equatorial Earth radius.
!                    A is analogous to the L value in invariant coordinates.
!            ALAT  = Apex Lat. (deg)
!            ALON  = Apex Lon. (deg)
!            XMAG  = Geomagnetic field component (gauss): north
!            YMAG  = Geomagnetic field component (gauss): east
!            ZMAG  = Geomagnetic field component (gauss): down
!            F     = Geomagnetic field magnitude (gauss)
!
!          Trace the geomagnetic field line from the given location to find the
!          apex of the field line.  Before starting iterations to trace along
!          the field line: (1) Establish a step size (DS, arc length in km)
!          based on the geomagnetic dipole latitude; (2) determine the step
!          direction from the sign of the vertical component of the geomagnetic
!          field; and (3) convert to geocentric cartesian coordinates.  Each
!          iteration increments a step count (NSTP) and calls ITRACE to move
!          along the the field line until reaching the iteration count limit
!          (MAXS) or passing the apex (IAPX=2) and then calling FNDAPX to
!          determine the apex location from the last three step locations
!          (YAPX); however, if reaching the iteration limit, apex coordinates
!          are calculated by DIPAPX which assumes a simplified dipole field.
!
!          COMMON BLOCKS:
!            COMMON /APXIN/   YAPX(3,3)
!            COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
!            COMMON /FLDCOMD/ BX, BY, BZ, BB
!            COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
!
!          APXIN has step locations determined in ITRACE:
!            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
!                    three points about the apex.  Set in subroutine ITRACE.
!
!          DIPOLE has IGRF variables obtained from routines in magfld.f:
!            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
!            ELON  = East longitude of geomagnetic dipole north pole (deg)
!            VP    = Magnitude (T-m) of dipole component of magnetic potential at
!                    geomagnetic pole and geocentric radius of 6371.0088 km
!            CTP   = cosine of COLAT
!            STP   = sine   of COLAT
!
!          FLDCOMD has geomagnetic field at current trace point:
!            BX    = X component (Gauss)
!            BY    = Y component (Gauss)
!            BZ    = Z component (Gauss)
!            BB    = Magnitude   (Gauss)
!
!          ITRA has field line tracing variables determined in LINAPX:
!            NSTP  = Step count.
!            Y     = Array containing current tracing point cartesian coordinates.
!            YOLD  = Array containing previous tracing point cartesian coordinates.
!            SGN   = Determines direction of trace.
!            DS    = Step size (arc length in km).
!
!          REFERENCES:
!            Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
!            Greenbelt, Maryland
!
!          EXTERNALS:
!            GD2CART = Convert geodetic to geocentric cartesian coordinates (in magfld.f)
!            CONVRT  = Convert geodetic to geocentric cylindrical or geocentric spherical
!                      and back (in magfld.f).
!            FELDG   = Obtain IGRF magnetic field components (in magfld.f).
!            ITRACE  = Follow a geomagnetic field line
!            DIPAPX  = Compute apex coordinates assuming a geomagnetic dipole field
!            FNDAPX  = Compute apex coordinates from the last three traced field line points
!
!------------------------------------------------------------------------------
!          HISTORY:
!          Oct 1973: Initial version completed on the 29th by Wally Clark, NOAA
!                    ERL Lab.
!          Feb 1988: Revised on the 1st by Harsh Anand Passi, NCAR.
!          Aug 1994: Revision by A. D. Richmond, NCAR.
!          Nov 2009: Change definition of earth's mean radius (RE) from 6371.2
!                    to the WGS84 value (6371.0088), by J.T. Emmert, NRL.

      PARAMETER (MAXS = 200, RTOD = 57.2957795130823,   RE =6371.0088, &
                             DTOR = .01745329251994330, REQ=6378.137)
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /APXIN/   YAPX(3,3)
      COMMON /DIPOLE/  COLAT,ELON,VP,CTP,STP
      COMMON /ITRA/    NSTP, Y(3), YP(3), SGN, DS

!          Set step size based on the geomagnetic dipole latitude of the starting point
      CALL CONVRT (2,GDLAT,ALT,GCLAT,R)
      SINGML = CTP*SIN(GCLAT*DTOR) + STP*COS(GCLAT*DTOR)* &
                                                   COS((GLON-ELON)*DTOR)
!          May 1999: avoid possible divide by zero (when SINGML = 1.): the old version
!          limited DS to its value at 60 deg GM latitude with: DS = .06*R/(1.-SINGML*SINGML) - 370.
!                                                              IF (DS .GT. 1186.) DS = 1186.
      CGML2 = AMAX1 (0.25,1.-SINGML*SINGML)
      DS = .06*R/CGML2 - 370.

!          Initialize YAPX array
      DO 4 J=1,3
      DO 4 I=1,3
    4 YAPX(I,J) = 0.

!          Convert from geodetic to earth centered cartesian coordinates
      CALL GD2CART (GDLAT,GLON,ALT,Y(1),Y(2),Y(3))
      NSTP = 0

!          Get magnetic field components to determine the direction for
!          tracing the field line
      CALL FELDG (1,GDLAT,GLON,ALT,XMAG,YMAG,ZMAG,F)
      SGN = SIGN (1.,-ZMAG)

!          Use cartesian coordinates to get magnetic field components
!          (from which gradients steer the tracing)
   10 CALL FELDG (2, Y(1)/RE,Y(2)/RE,Y(3)/RE, BX,BY,BZ,BB)
      NSTP = NSTP + 1

      IF (NSTP .LT. MAXS) THEN
      CALL ITRACE (IAPX)                               ! trace along fie
      IF (IAPX .EQ. 1) GO TO 10
      CALL FNDAPX (ALT,ZMAG,A,ALAT,ALON)               ! (IAPX=2) => pas
      ELSE
      RHO = SQRT (Y(1)*Y(1) + Y(2)*Y(2))               ! too many steps;
      CALL CONVRT (3,XLAT,HT,RHO,Y(3))
      XLON = RTOD*ATAN2 (Y(2),Y(1))
      CALL FELDG (1,XLAT,XLON,HT,BNRTH,BEAST,BDOWN,BABS)
      CALL DIPAPX  (XLAT,XLON,HT,BNRTH,BEAST,BDOWN,A,ALON)
      ALAT = -SGN*RTOD*ACOS (SQRT(1./A))
      ENDIF

      RETURN
      END

      SUBROUTINE ITRACE (IAPX)

!          Follow a geomagnetic field line until passing its apex
!
!          INPUTS:
!            (all are in common blocks)
!          OUTPUTS:
!            IAPX = 2 (when apex passed) or 1 (not)
!
!          This uses the 4-point Adams formula after initialization.
!          First 7 iterations advance point by 3 steps.
!
!          COMMON BLOCKS:
!            COMMON /APXIN/   YAPX(3,3)
!            COMMON /FLDCOMD/ BX, BY, BZ, BB
!            COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
!
!          APXIN has step locations determined in ITRACE:
!            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
!                    three points about the apex.  Set in subroutine ITRACE.
!
!          FLDCOMD has geomagnetic field at current trace point:
!            BX    = X component (Gauss)
!            BY    = Y component (Gauss)
!            BZ    = Z component (Gauss)
!            BB    = Magnitude   (Gauss)
!
!          ITRA has field line tracing variables determined in LINAPX:
!            NSTP  = Step count.
!            Y     = Array containing current tracing point cartesian coordinates.
!            YOLD  = Array containing previous tracing point cartesian coordinates.
!            SGN   = Determines direction of trace.
!            DS    = Step size (arc length in km).
!
!          REFERENCES:
!            Stassinopoulos E. G. , Mead Gilbert D., X-841-72-17 (1971) GSFC,
!            Greenbelt, Maryland
!------------------------------------------------------------------------------
!          HISTORY:
!          Oct 1973: Initial version completed on the 29th by W. Clark, NOAA ERL
!                    Laboratory.
!          Feb 1988: Revised by H. Passi, NCAR.
!          Apr 2004: Replace computed GO TO with IF blocks because some compilers
!                    are threatening to remove this old feature
!

      COMMON /APXIN/   YAPX(3,3)
      COMMON /FLDCOMD/ BX, BY, BZ, BB
      COMMON /ITRA/    NSTP, Y(3), YOLD(3), SGN, DS
      DIMENSION YP(3,4)
      SAVE

!          Statement function
      RDUS(D,E,F) = SQRT (D**2 + E**2 + F**2)

      IAPX = 1

!          Cartesian component magnetic field (partial) derivitives steer the trace
      YP(1,4) = SGN*BX/BB
      YP(2,4) = SGN*BY/BB
      YP(3,4) = SGN*BZ/BB

      IF (NSTP .LE. 7) THEN
      DO 10 I=1,3
      IF (NSTP .EQ. 1) THEN
        D2        = DS/2.
        D6        = DS/6.
        D12       = DS/12.
        D24       = DS/24.
        YP(I,1)   = YP(I,4)
        YOLD(I)   = Y(I)
        YAPX(I,1) = Y(I)
        Y(I)      = YOLD(I) + DS*YP(I,1)

      ELSE IF (NSTP .EQ. 2) THEN
        YP(I,2) = YP(I,4)
        Y(I)    = YOLD(I) + D2*(YP(I,2)+YP(I,1))

      ELSE IF (NSTP .EQ. 3) THEN
        Y(I) = YOLD(I) + D6*(2.*YP(I,4)+YP(I,2)+3.*YP(I,1))

      ELSE IF (NSTP .EQ. 4) THEN
        YP(I,2)   = YP(I,4)
        YAPX(I,2) = Y(I)
        YOLD(I)   = Y(I)
        Y(I)      = YOLD(I) + D2*(3.*YP(I,2)-YP(I,1))

      ELSE IF (NSTP .EQ. 5) THEN
        Y(I) = YOLD(I) + D12*(5.*YP(I,4)+8.*YP(I,2)-YP(I,1))

      ELSE IF (NSTP .EQ. 6) THEN
        YP(I,3)   = YP(I,4)
        YOLD(I)   = Y(I)
        YAPX(I,3) = Y(I)
        Y(I)      = YOLD(I) + D12*(23.*YP(I,3)-16.*YP(I,2)+5.*YP(I,1))

      ELSE IF (NSTP .EQ. 7) THEN
        YAPX(I,1) = YAPX(I, 2)
        YAPX(I,2) = YAPX(I, 3)
        Y(I)      = YOLD(I) + D24*(9.*YP(I,4) + 19.*YP(I,3) - &
                                     5.*YP(I,2) +     YP(I,1))
        YAPX(I,3) = Y(I)
      ENDIF
   10   CONTINUE
      IF (NSTP .EQ. 6 .OR. NSTP .EQ. 7) THEN        ! signal if apex pas
        RC = RDUS (YAPX(1,3), YAPX(2,3), YAPX(3,3))
        RP = RDUS (YAPX(1,2), YAPX(2,2), YAPX(3,2))
        IF (RC .LT. RP) IAPX = 2
      ENDIF

      ELSE                 ! NSTP > 7

      DO 30 I=1,3
      YAPX(I,1) = YAPX(I,2)
      YAPX(I,2) = Y(I)
      YOLD(I)   = Y(I)
      Y(I)      = YOLD(I) + D24*(55.*YP(I,4) - 59.*YP(I,3) + &
                                   37.*YP(I,2) -  9.*YP(I,1))
      YAPX(I,3) = Y(I)

      DO 20 J=1,3
   20   YP(I,J) = YP(I,J+1)
   30   CONTINUE
      RC = RDUS (   Y(1),    Y(2),    Y(3))
      RP = RDUS (YOLD(1), YOLD(2), YOLD(3))
      IF (RC .LT. RP) IAPX = 2
      ENDIF

      RETURN
      END

       SUBROUTINE FNDAPX (ALT,ZMAG,A,ALAT,ALON)

!          Find apex coordinates once tracing (in subroutine ITRACE) has
!          signalled that the apex has been passed.
!          INPUTS:
!            ALT  = Altitude of starting point
!            ZMAG = Downward component of geomagnetic field at starting point
!          OUTPUT
!            A    = Apex radius, defined as (Apex height + Req)/Req, where
!                   Req = equatorial Earth radius.
!                   A is analogous to the L value in invariant coordinates.
!            ALAT = Apex Lat. (deg)
!            ALON = Apex Lon. (deg)
!
!          COMMON BLOCKS:
!            COMMON /APXIN/  YAPX(3,3)
!            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
!
!          APXIN has step locations determined in ITRACE:
!            YAPX  = Matrix of cartesian coordinates (loaded columnwise) of the
!                    three points about the apex.  Set in subroutine ITRACE.
!
!          DIPOLE has IGRF variables obtained from routines in magfld.f:
!            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
!            ELON  = East longitude of geomagnetic dipole north pole (deg)
!            VP    = Magnitude (T-m) of dipole component of magnetic potential at
!                    geomagnetic pole and geocentric radius of 6371.0088 km
!            CTP   = cosine of COLAT
!            STP   = sine   of COLAT
!
!          EXTERNALS:
!            FINT = Second degree interpolation routine
!------------------------------------------------------------------------------
!          HISTORY:
!          Oct 1973: Initial version completed on the 23rd by Clark, W., NOAA
!                    Boulder.
!          Aug 1994: Revision on the 3rd by A.D. Richmond, NCAR
!          Apr 2004: Repair problem noted by Dan Weimer where the apex location
!                    produced by FINT may still have a non-zero vertical magnetic
!                    field component.

      PARAMETER (RTOD = 57.2957795130823, &
                 DTOR = .01745329251994330, REQ=6378.137)
      COMMON /APXIN/  YAPX(3,3)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
      DIMENSION BD(3), Y(3)

!          Get geodetic height and vertical (downward) component of the magnetic
!          field at last three points found by ITRACE
      DO 10 I=1,3
      RHO  = SQRT (YAPX(1,I)**2 + YAPX(2,I)**2)
      CALL CONVRT (3,GDLT,HT, RHO,YAPX(3,I))
      GDLN = RTOD*ATAN2 (YAPX(2,I),YAPX(1,I))
   10 CALL FELDG (1,GDLT,GDLN,HT, BN,BE,BD(I),BMAG)

!          Interpolate to where Bdown=0 to find cartesian coordinates at dip equator
      NITR = 0
   20 Y(1) = FINT (BD(1),BD(2),BD(3),YAPX(1,1),YAPX(1,2),YAPX(1,3), 0.)
      Y(2) = FINT (BD(1),BD(2),BD(3),YAPX(2,1),YAPX(2,2),YAPX(2,3), 0.)
      Y(3) = FINT (BD(1),BD(2),BD(3),YAPX(3,1),YAPX(3,2),YAPX(3,3), 0.)

!          Insure negligible Bdown or
!
!            |Bdown/Btot| < 2.E-6
!
!          For instance, Bdown must be less than 0.1 nT at low altitudes where
!          Btot ~ 50000 nT.  This ratio can be exceeded when interpolation is
!          not accurate; i.e., when the middle of the three points interpolated
!          is too far from the dip equator.  The three points were initially
!          defined with equal spacing by ITRACE, so replacing point 2 with the
!          most recently fit location will reduce the interpolation span.
      RHO  = SQRT (Y(1)**2 + Y(2)**2)
      GDLN = RTOD*ATAN2 (Y(2),Y(1))
      CALL CONVRT (3,GDLT,HTA, RHO,Y(3))
      CALL FELDG (1,GDLT,GDLN,HTA, BNA,BEA,BDA,BA)
      ABDOB = ABS(BDA/BA)

      IF (ABDOB .GT. 2.E-6) THEN
      IF (NITR .LT. 4) THEN        ! 4 was chosen because tests rarely r
        NITR      = NITR + 1
        YAPX(1,2) = Y(1)
        YAPX(2,2) = Y(2)
        YAPX(3,2) = Y(3)
        BD(2)     = BDA
        GO TO 20
      ELSE
        WRITE (0,'(''APEX: Imprecise fit of apex: |Bdown/B| ='',1PE7.1)') ABDOB
      ENDIF
      ENDIF

!          Ensure altitude of the Apex is at least the initial altitude when
!          defining the Apex radius then use it to define the Apex latitude whose
!          hemisphere (sign) is inferred from the sign of the dip angle at the
!          starting point
      A = (REQ + AMAX1(ALT,HTA)) / REQ
      IF (A .LT. 1.) THEN
      WRITE (0,'(''APEX: A can not be less than 1; A, REQ, HTA: '',1P3E15.7)') A,REQ,HTA
      CALL EXIT (1)
      ENDIF
      RASQ = ACOS (SQRT(1./A))*RTOD
      ALAT = SIGN (RASQ,ZMAG)

!          ALON is the dipole longitude of the apex and is defined using
!          spherical coordinates where
!            GP   = geographic pole.
!            GM   = geomagnetic pole (colatitude COLAT, east longitude ELON).
!            XLON = longitude of apex.
!            TE   = colatitude of apex.
!            ANG  = longitude angle from GM to apex.
!            TP   = colatitude of GM.
!            TF   = arc length between GM and apex.
!            PA   = ALON be geomagnetic longitude, i.e., Pi minus angle measured
!                   counterclockwise from arc GM-apex to arc GM-GP.
!          then, spherical-trigonometry formulas for the functions of the angles
!          are as shown below.  Notation uses C=cos, S=sin and STFCPA = sin(TF) * cos(PA),
!                                                              STFSPA = sin(TF) * sin(PA)
      XLON = ATAN2 (Y(2),Y(1))
      ANG  = XLON-ELON*DTOR
      CANG = COS (ANG)
      SANG = SIN (ANG)
      R    = SQRT (Y(1)**2 + Y(2)**2 + Y(3)**2)
      CTE  = Y(3)/R
      STE  = SQRT (1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      ALON = ATAN2 (STFSPA,STFCPA)*RTOD
      RETURN
      END

      SUBROUTINE DIPAPX (GDLAT,GDLON,ALT,BNORTH,BEAST,BDOWN, A,ALON)

!          Compute A, ALON from local magnetic field using dipole and spherical
!          approximation.
!
!          INPUTS:
!            GDLAT  = geodetic latitude, degrees
!            GDLON  = geodetic longitude, degrees
!            ALT    = altitude, km
!            BNORTH = geodetic northward magnetic field component (any units)
!            BEAST  = eastward magnetic field component
!            BDOWN  = geodetic downward magnetic field component
!          OUTPUTS:
!            A      = apex radius, 1 + h_A/R_eq
!            ALON   = apex longitude, degrees
!
!          Use spherical coordinates and define:
!            GP    = geographic pole.
!            GM    = geomagnetic pole (colatitude COLAT, east longitude ELON).
!            G     = point at GDLAT,GDLON.
!            E     = point on sphere below apex of dipolar field line passing
!                    through G.
!            TD    = dipole colatitude of point G, found by applying dipole
!                    formula for dip angle to actual dip angle.
!            B     = Pi plus local declination angle.  B is in the direction
!                    from G to E.
!            TG    = colatitude of G.
!            ANG   = longitude angle from GM to G.
!            TE    = colatitude of E.
!            TP    = colatitude of GM.
!            A     = longitude angle from G to E.
!            APANG = A + ANG
!            PA    = geomagnetic longitude, i.e., Pi minus angle measured
!                    counterclockwise from arc GM-E to arc GM-GP.
!            TF    = arc length between GM and E.
!          Then, using notation C=cos, S=sin, COT=cot, spherical-trigonometry
!          formulas for the functions of the angles are as shown below.  Note:
!            STFCPA = sin(TF) * cos(PA)
!            STFSPA = sin(TF) * sin(PA)
!
!          COMMON BLOCKS:
!            COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP
!
!          DIPOLE has IGRF variables obtained from routines in magfld.f:
!            COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
!            ELON  = East longitude of geomagnetic dipole north pole (deg)
!            VP    = Magnitude (T-m) of dipole component of magnetic potential at
!                    geomagnetic pole and geocentric radius of 6371.0088 km
!            CTP   = cosine of COLAT
!            STP   = sine   of COLAT
!------------------------------------------------------------------------------
!          HISTORY:
!          May 1994:  Completed on the 1st by A. D. Richmond
!          Nov 2009: Change definition of earth's mean radius (RE) from 6371.2
!                    to the WGS84 value (6371.0088), by J.T. Emmert, NRL.

      PARAMETER (RTOD = 57.2957795130823,   RE =6371.0088, &
                 DTOR = .01745329251994330, REQ=6378.137)
      COMMON /DIPOLE/ COLAT,ELON,VP,CTP,STP

      BHOR = SQRT(BNORTH*BNORTH + BEAST*BEAST)
      IF (BHOR .EQ. 0.) THEN
      ALON = 0.
      A    = 1.E34
      RETURN
      ENDIF
      COTTD  = BDOWN*.5/BHOR
      STD    = 1./SQRT(1.+COTTD*COTTD)
      CTD    = COTTD*STD
      SB     = -BEAST /BHOR
      CB     = -BNORTH/BHOR
      CTG    = SIN (GDLAT*DTOR)
      STG    = COS (GDLAT*DTOR)
      ANG    = (GDLON-ELON)*DTOR
      SANG   = SIN(ANG)
      CANG   = COS(ANG)
      CTE    = CTG*STD + STG*CTD*CB
      STE    = SQRT(1. - CTE*CTE)
      SA     = SB*CTD/STE
      CA     = (STD*STG - CTD*CTG*CB)/STE
      CAPANG = CA*CANG - SA*SANG
      SAPANG = CA*SANG + SA*CANG
      STFCPA = STE*CTP*CAPANG - CTE*STP
      STFSPA = SAPANG*STE
      ALON = ATAN2 (STFSPA,STFCPA)*RTOD
      R    = ALT + RE
      HA   = ALT + R*COTTD*COTTD
      A    = 1. + HA/REQ
      RETURN
      END

      FUNCTION FINT (X1,X2,X3,Y1,Y2,Y3, XFIT)
!          Second degree interpolation used by FNDAPX
!          INPUTS:
!            X1   = point 1 ordinate value
!            X2   = point 2 ordinate value
!            X3   = point 3 ordinate value
!            Y1   = point 1 abscissa value
!            Y2   = point 2 abscissa value
!            Y3   = point 3 abscissa value
!            XFIT = ordinate value to fit
!          RETURNS:
!            YFIT = abscissa value corresponding to XFIT
!
!          MODIFICATIONS:
!          Apr 2004: Change from subroutine to function, rename variables and
!                    add intermediates which are otherwise calculated twice
      X12 = X1-X2
      X13 = X1-X3
      X23 = X2-X3
      XF1 = XFIT-X1
      XF2 = XFIT-X2
      XF3 = XFIT-X3

      FINT = (Y1*X23*XF2*XF3 - Y2*X13*XF1*XF3 + Y3*X12*XF1*XF2) / &
                                                           (X12*X13*X23)
      RETURN
      END

      SUBROUTINE GD2CART (GDLAT,GLON,ALT,X,Y,Z)
!          Convert geodetic to cartesian coordinates by calling CONVRT
!          940503 A. D. Richmond
      PARAMETER (DTOR = 0.01745329251994330)
      CALL CONVRT (1,GDLAT,ALT,RHO,Z)
      ANG = GLON*DTOR
      X = RHO*COS(ANG)
      Y = RHO*SIN(ANG)
      RETURN
      END

      SUBROUTINE CONVRT (I,GDLAT,ALT,X1,X2)
!          Convert space point from geodetic to geocentric or vice versa.
!
!          I is an input flag controlling the meaning and direction of the
!            remaining formal arguments:
!
!          I = 1  (convert from geodetic to cylindrical geocentric)
!            INPUTS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!            RETURNS:
!              X1    = Distance from Earth's rotation axis (km)
!              X2    = Distance above (north of) Earth's equatorial plane (km)
!
!          I = 2  (convert from geodetic to spherical geocentric)
!            INPUTS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!            RETURNS:
!              X1    = Geocentric latitude (deg)
!              X2    = Geocentric distance (km)
!
!          I = 3  (convert from cylindrical geocentric to geodetic)
!            INPUTS:
!              X1    = Distance from Earth's rotation axis (km)
!              X2    = Distance from Earth's equatorial plane (km)
!            RETURNS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!
!          I = 4  (convert from spherical geocentric to geodetic)
!            INPUTS:
!              X1    = Geocentric latitude (deg)
!              X2    = Geocentric distance (km)
!            RETURNS:
!              GDLAT = Geodetic latitude (deg)
!              ALT   = Altitude above reference ellipsoid (km)
!
!
!          HISTORY:
!          940503 (A. D. Richmond):  Based on a routine originally written
!          by V. B. Wickwar.
!
!          Mar 2004: (Barnes) Revise spheroid definition to WGS-1984 to conform
!          with IGRF-9 release (EOS Volume 84 Number 46 November 18 2003).
!
!          REFERENCE: ASTRON. J. VOL. 66, p. 15-16, 1961

!          E2  = square of eccentricity of ellipse
!          REP = earth's polar      radius (km)
!          REQ = earth's equatorial radius (km)
      PARAMETER (RTOD = 57.2957795130823, DTOR = 0.01745329251994330, &
                 REP  = 6356.752, REQ = 6378.137, E2 = 1.-(REP/REQ)**2, &
           E4 = E2*E2, E6 = E4*E2, E8 = E4*E4, OME2REQ = (1.-E2)*REQ, &
           A21 =     (512.*E2 + 128.*E4 + 60.*E6 + 35.*E8)/1024. , &
           A22 =     (                        E6 +     E8)/  32. , &
           A23 = -3.*(                     4.*E6 +  3.*E8)/ 256. , &
           A41 =    -(           64.*E4 + 48.*E6 + 35.*E8)/1024. , &
           A42 =     (            4.*E4 +  2.*E6 +     E8)/  16. , &
           A43 =                                   15.*E8 / 256. , &
           A44 =                                      -E8 /  16. , &
           A61 =  3.*(                     4.*E6 +  5.*E8)/1024. , &
           A62 = -3.*(                        E6 +     E8)/  32. , &
           A63 = 35.*(                     4.*E6 +  3.*E8)/ 768. , &
           A81 =                                   -5.*E8 /2048. , &
           A82 =                                   64.*E8 /2048. , &
           A83 =                                 -252.*E8 /2048. , &
           A84 =                                  320.*E8 /2048. )

      IF (I .GE. 3) GO TO 300

!          Geodetic to geocentric

!          Compute RHO,Z
      SINLAT = SIN(GDLAT*DTOR)
      COSLAT = SQRT(1.-SINLAT*SINLAT)
      D      = SQRT(1.-E2*SINLAT*SINLAT)
      Z      = (ALT+OME2REQ/D)*SINLAT
      RHO    = (ALT+REQ/D)*COSLAT
      X1 = RHO
      X2 = Z
      IF (I .EQ. 1) RETURN

!          Compute GCLAT,RKM
      RKM   = SQRT(Z*Z + RHO*RHO)
      GCLAT = RTOD*ATAN2(Z,RHO)
      X1 = GCLAT
      X2 = RKM
      RETURN

!          Geocentric to geodetic
  300 IF (I .EQ. 3) THEN
         RHO = X1
         Z = X2
         RKM = SQRT(Z*Z+RHO*RHO)
         SCL = Z/RKM
         GCLAT = ASIN(SCL)*RTOD
      ELSEIF (I .EQ. 4) THEN
         GCLAT = X1
         RKM = X2
         SCL = SIN(GCLAT*DTOR)
      ELSE
         RETURN
      ENDIF

      RI = REQ/RKM
      A2 = RI*(A21+RI*(A22+RI* A23))
      A4 = RI*(A41+RI*(A42+RI*(A43+RI*A44)))
      A6 = RI*(A61+RI*(A62+RI* A63))
      A8 = RI*(A81+RI*(A82+RI*(A83+RI*A84)))
      CCL = SQRT(1.-SCL*SCL)
      S2CL = 2.*SCL*CCL
      C2CL = 2.*CCL*CCL-1.
      S4CL = 2.*S2CL*C2CL
      C4CL = 2.*C2CL*C2CL-1.
      S8CL = 2.*S4CL*C4CL
      S6CL = S2CL*C4CL+C2CL*S4CL
      DLTCL = S2CL*A2+S4CL*A4+S6CL*A6+S8CL*A8
      GDLAT = DLTCL*RTOD+GCLAT
      SGL = SIN(GDLAT*DTOR)
      ALT = RKM*COS(DLTCL)-REQ*SQRT(1.-E2*SGL*SGL)
      RETURN
      END


      SUBROUTINE COFRM (DATE)
!          Define the International Geomagnetic Reference Field (IGRF) as a
!          scalar potential field using a truncated series expansion with
!          Schmidt semi-normalized associated Legendre functions of degree n and
!          order m.  The polynomial coefficients are a function of time and are
!          interpolated between five year epochs or extrapolated at a constant
!          rate after the last epoch.
!
!          INPUTS:
!            DATE = yyyy.fraction (UT)
!          OUTPUTS (in comnon block MAGCOF):
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed in COFRM
!
!          It is fatal to supply a DATE before the first epoch.  A warning is
!          issued to Fortran unit 0 (stderr) if DATE is later than the
!          recommended limit, five years after the last epoch.
!
!          HISTORY (blame):
!          Apr 1983:  Written by Vincent B. Wickwar (Utah State Univ.) including
!          secular variation acceleration rate set to zero in case the IGRF
!          definition includes such second time derivitives.  The maximum degree
!          (n) defined was 10.
!
!          Jun 1986:  Updated coefficients adding Definitive Geomagnetic Reference
!          Field (DGRF) for 1980 and IGRF for 1985 (EOS Volume 7 Number 24).  The
!          designation DGRF means coefficients will not change in the future
!          whereas IGRF coefficients are interim pending incorporation of new
!          magnetometer data.  Common block MAG was replaced by MAGCOF, thus
!          removing variables not used in subroutine FELDG.  (Roy Barnes)
!
!          Apr 1992 (Barnes):  Added DGRF 1985 and IGRF 1990 as given in EOS
!          Volume 73 Number 16 April 21 1992.  Other changes were made so future
!          updates should:  (1) Increment NDGY; (2) Append to EPOCH the next IGRF
!          year; (3) Append the next DGRF coefficients to G1DIM and H1DIM; and (4)
!          replace the IGRF initial values (G0, GT) and rates of change indices
!          (H0, HT).
!
!          Apr 1994 (Art Richmond): Computation of GV added, for finding magnetic
!          potential.
!
!          Aug 1995 (Barnes):  Added DGRF for 1990 and IGRF for 1995, which were
!          obtained by anonymous ftp to geomag.gsfc.nasa.gov (cd pub, mget table*)
!          as per instructions from Bob Langel (langel@geomag.gsfc.nasa.gov) with
!          problems reported to baldwin@geomag.gsfc.nasa.gov.
!
!          Oct 1995 (Barnes):  Correct error in IGRF-95 G 7 6 and H 8 7 (see email
!          in folder).  Also found bug whereby coefficients were not being updated
!          in FELDG when IENTY did not change so ICHG was added to flag date
!          changes.  Also, a vestigial switch (IS) was removed from COFRM; it was
!          always zero and involved 3 branch if statements in the main polynomial
!          construction loop now numbered 200.
!
!          Feb 1999 (Barnes):  Explicitly initialize GV(1) in COFRM to avoid the
!          possibility of compiler or loader options initializing memory to
!          something else (e.g., indefinite).  Also simplify the algebra in COFRM
!          with no effect on results.
!
!          Mar 1999 (Barnes):  Removed three branch if's from FELDG and changed
!          statement labels to ascending order.
!
!          Jun 1999 (Barnes):  Corrected RTOD definition in GD2CART.
!
!          May 2000 (Barnes):  Replace IGRF 1995, add IGRF 2000, and extend the
!          earlier DGRF's back to 1900.  The coefficients came from an NGDC web
!          page.  Related documentation is in $APXROOT/docs/igrf.2000.*  where
!          $APXROOT, defined by 'source envapex', is traditionally ~bozo/apex).
!
!          Mar 2004 (Barnes):  Replace 1995 and 2000 coefficients; now both are
!          DGRF.  Coefficients for 2000 are degree 13 with precision increased to
!          tenths nT and accommodating this instigated changes:  (1) degree (NMAX)
!          is now a function of epoch (NMXE) to curtail irrelevant looping over
!          unused high order terms (n > 10 in epochs before 2000) when calculating
!          GB; (2) expand coefficients data statement layout for G1D and H1D,
!          formerly G1DIM and H1DIM; (3) omit secular variation acceleration terms
!          which were always zero; (4) increase array dimensions in common block
!          MAGCOF and associated arrays G and H in FELDG; (5) change earth's shape
!          in CONVRT from the IAU-1966 to the WGS-1984 spheroid; (6) eliminate
!          reference to 'definitive' in variables in COFRM which were not always
!          definitive; (7) change G to GB in COFRM s.t. arrays GB and GV in common
!          block MAGCOF are consistently named in all subroutines; (8) remove
!          unused constants in all five subroutines.  See EOS Volume 84 Number 46
!          November 18 2003, www.ngdc.noaa.gov/IAGA/vmod/igrf.html or local files
!          $APXROOT/docs/igrf.2004.*
!
!          Sept. 2005 (Maute): update with IGRF10 from
!          http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html use script
!          ~maute/apex.d/apex_update/igrf2f Note that the downloaded file the start
!          column of the year in the first line has to be before the start of each
!          number in the same column
!
!          July 2010 (Stolle, cst@space.dtu.dk): update with IGRF11 from
!          http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
!          used script
!          /home/cst/CASTOR/APEX/print_new_coeff.pro
!          increased parameter NEPO = 22 > NEPO = 23
!          and adjusted EPOCH and NMXE
!
!          February 2015 (Rauberg, rauberg@gfz-potsdam.de): update with IGRF12 from
!          http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html
!          used script
!          /home/mag/rauberg/test/prj/me/apex/apex3_2020_igrf12_lib/igrf2apex
!          g(n,m), h(n,m)
!          Updated to DGRF2010
!          Added IGRF2015
!          Updated SV
!          increased parameter NEPO = 23 > NEPO = 24
!          and adjusted EPOCH and NMXE from 23 to 24 elements
!
!          January 2020 (Rauberg, rauberg@gfz-potsdam.de): update with IGRF13 from
!          https://www.ngdc.noaa.gov/IAGA/vmod/coeffs/igrf13coeffs.txt
!          used script
!          /home/rauberg/test/prj/me/apex/apex3_igrf13_nils/igrf_reform.py
!          g(n,m), h(n,m)
!          Updated IGRF2010 to DGRF2010 (was wrong in version 2015)
!          Updated IGRF2015 to DGRF2015
!          Added IGRF2020
!          Updated SV (valid up to 2025)
!          increased parameter NEPO = 24 > NEPO = 25
!          and adjusted EPOCH (2020) and NMXE (13) from 24 to 25 elements
!

      DOUBLE PRECISION F,F0
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
      DATA ICHG /-99999/

!          NEPO = Number of epochs
!          NGH  = Single dimensioned array size of 2D version (GYR or HYR)
!          NGHT = Single dimensioned array size of 2D version (GT  or HT)
      PARAMETER (NEPO = 26, NGH = 225*NEPO, NGHT = 225)
      DIMENSION GYR(15,15,NEPO), HYR(15,15,NEPO), EPOCH(NEPO), &
                GT (15,15),      HT (15,15),       NMXE(NEPO), &
                GY1D(NGH),       HY1D(NGH), &
                GT1D(NGHT),      HT1D(NGHT)
      EQUIVALENCE (GYR(1,1,1),GY1D(1)), (HYR(1,1,1),HY1D(1)), &
                  (GT (1,1),  GT1D(1)), (HT (1,1),  HT1D(1))

      SAVE DATEL, EPOCH, NMXE, GYR, HYR, GT, HT, GY1D, HY1D, GT1D, HT1D
      DATA DATEL /-999./, &
           EPOCH / 1900, 1905, 1910, 1915, 1920, 1925, 1930, 1935, 1940, &
                   1945, 1950, 1955, 1960, 1965, 1970, 1975, 1980, 1985, &
                   1990, 1995, 2000, 2005, 2010, 2015, 2020, 2025/, &
           NMXE  /   10,   10,   10,   10,   10,   10,   10,   10,   10, &
                     10,   10,   10,   10,   10,   10,   10,   10,   10, &
                 10,   13,   13,   13,   13,   13,   13, 13/

!          g(n,m) for 1900
!          Fields across a line are (degree) n=1,13; lines are (order) m=0,13 as indicated
!          in column 6; e.g., for 1965 g(n=3,m=0) = 1297 or g(n=6,m=6) = -111
!
!           1       2       3      4      5      6      7      8      9
!                                        10     11     12     13          (n)
      DATA (GY1D(I),I=1,145) /0, &
        -31543,   -677,  1022,   876,  -184,    63,    70,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2298,   2905, -1469,   628,   328,    61,   -55,     8,    10, &
                                         -4,     0,     0,     0,   3*0, &
                   924,  1256,   660,   264,   -11,     0,    -4,     1, &
                                          2,     0,     0,     0,   4*0, &
                          572,  -361,     5,  -217,    34,    -9,   -11, &
                                         -5,     0,     0,     0,   5*0, &
                                 134,   -86,   -58,   -41,     1,    12, &
                                         -2,     0,     0,     0,   6*0, &
                                        -16,    59,   -21,     2,     1, &
                                          6,     0,     0,     0,   7*0, &
                                               -90,    18,    -9,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        6,     5,     2, &
                                          0,     0,     0,     0,   9*0, &
                                                               8,    -1, &
                                          2,     0,     0,     0,  10*0, &
                                                                     -1/
      DATA (GY1D(I),I=146,225) / &
                                          2,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1905
      DATA (GY1D(I),I=226,370) /0, &
        -31464,   -728,  1037,   880,  -192,    62,    70,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2298,   2928, -1494,   643,   328,    60,   -54,     8,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1041,  1239,   653,   259,   -11,     0,    -4,     1, &
                                          2,     0,     0,     0,   4*0, &
                          635,  -380,    -1,  -221,    33,    -9,   -11, &
                                         -5,     0,     0,     0,   5*0, &
                                 146,   -93,   -57,   -41,     1,    12, &
                                         -2,     0,     0,     0,   6*0, &
                                        -26,    57,   -20,     2,     1, &
                                          6,     0,     0,     0,   7*0, &
                                               -92,    18,    -8,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        6,     5,     2, &
                                          0,     0,     0,     0,   9*0, &
                                                               8,     0, &
                                          2,     0,     0,     0,  10*0, &
                                                                     -1/
      DATA (GY1D(I),I=371,450) / &
                                          2,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1910
      DATA (GY1D(I),I=451,595) /0, &
        -31354,   -769,  1058,   884,  -201,    62,    71,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2297,   2948, -1524,   660,   327,    58,   -54,     8,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1176,  1223,   644,   253,   -11,     1,    -4,     1, &
                                          2,     0,     0,     0,   4*0, &
                          705,  -400,    -9,  -224,    32,    -9,   -11, &
                                         -5,     0,     0,     0,   5*0, &
                                 160,  -102,   -54,   -40,     1,    12, &
                                         -2,     0,     0,     0,   6*0, &
                                        -38,    54,   -19,     2,     1, &
                                          6,     0,     0,     0,   7*0, &
                                               -95,    18,    -8,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        6,     5,     2, &
                                          0,     0,     0,     0,   9*0, &
                                                               8,     0, &
                                          2,     0,     0,     0,  10*0, &
                                                                     -1/
      DATA (GY1D(I),I=596,675) / &
                                          2,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1915
      DATA (GY1D(I),I=676,820) /0, &
        -31212,   -802,  1084,   887,  -211,    61,    72,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2306,   2956, -1559,   678,   327,    57,   -54,     8,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1309,  1212,   631,   245,   -10,     2,    -4,     1, &
                                          2,     0,     0,     0,   4*0, &
                          778,  -416,   -16,  -228,    31,    -9,   -11, &
                                         -5,     0,     0,     0,   5*0, &
                                 178,  -111,   -51,   -38,     2,    12, &
                                         -2,     0,     0,     0,   6*0, &
                                        -51,    49,   -18,     3,     1, &
                                          6,     0,     0,     0,   7*0, &
                                               -98,    19,    -8,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        6,     6,     2, &
                                          0,     0,     0,     0,   9*0, &
                                                               8,     0, &
                                          1,     0,     0,     0,  10*0, &
                                                                     -1/
      DATA (GY1D(I),I=821,900) / &
                                          2,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1920
      DATA (GY1D(I),I=901,1045) /0, &
        -31060,   -839,  1111,   889,  -221,    61,    73,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2317,   2959, -1600,   695,   326,    55,   -54,     7,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1407,  1205,   616,   236,   -10,     2,    -3,     1, &
                                          2,     0,     0,     0,   4*0, &
                          839,  -424,   -23,  -233,    29,    -9,   -11, &
                                         -5,     0,     0,     0,   5*0, &
                                 199,  -119,   -46,   -37,     2,    12, &
                                         -2,     0,     0,     0,   6*0, &
                                        -62,    44,   -16,     4,     1, &
                                          6,     0,     0,     0,   7*0, &
                                              -101,    19,    -7,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        6,     6,     2, &
                                          0,     0,     0,     0,   9*0, &
                                                               8,     0, &
                                          1,     0,     0,     0,  10*0, &
                                                                     -1/
      DATA (GY1D(I),I=1046,1125) / &
                                          3,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1925
      DATA (GY1D(I),I=1126,1270) /0, &
        -30926,   -893,  1140,   891,  -230,    61,    73,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2318,   2969, -1645,   711,   326,    54,   -54,     7,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1471,  1202,   601,   226,    -9,     3,    -3,     1, &
                                          2,     0,     0,     0,   4*0, &
                          881,  -426,   -28,  -238,    27,    -9,   -11, &
                                         -5,     0,     0,     0,   5*0, &
                                 217,  -125,   -40,   -35,     2,    12, &
                                         -2,     0,     0,     0,   6*0, &
                                        -69,    39,   -14,     4,     1, &
                                          6,     0,     0,     0,   7*0, &
                                              -103,    19,    -7,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        6,     7,     2, &
                                          0,     0,     0,     0,   9*0, &
                                                               8,     0, &
                                          1,     0,     0,     0,  10*0, &
                                                                     -1/
      DATA (GY1D(I),I=1271,1350) / &
                                          3,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1930
      DATA (GY1D(I),I=1351,1495) /0, &
        -30805,   -951,  1172,   896,  -237,    60,    74,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2316,   2980, -1692,   727,   327,    53,   -54,     7,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1517,  1205,   584,   218,    -9,     4,    -3,     1, &
                                          2,     0,     0,     0,   4*0, &
                          907,  -422,   -32,  -242,    25,    -9,   -12, &
                                         -5,     0,     0,     0,   5*0, &
                                 234,  -131,   -32,   -34,     2,    12, &
                                         -2,     0,     0,     0,   6*0, &
                                        -74,    32,   -12,     5,     1, &
                                          6,     0,     0,     0,   7*0, &
                                              -104,    18,    -6,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        6,     8,     3, &
                                          0,     0,     0,     0,   9*0, &
                                                               8,     0, &
                                          1,     0,     0,     0,  10*0, &
                                                                     -2/
      DATA (GY1D(I),I=1496,1575) / &
                                          3,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1935
      DATA (GY1D(I),I=1576,1720) /0, &
        -30715,  -1018,  1206,   903,  -241,    59,    74,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2306,   2984, -1740,   744,   329,    53,   -53,     7,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1550,  1215,   565,   211,    -8,     4,    -3,     1, &
                                          2,     0,     0,     0,   4*0, &
                          918,  -415,   -33,  -246,    23,    -9,   -12, &
                                         -5,     0,     0,     0,   5*0, &
                                 249,  -136,   -25,   -33,     1,    11, &
                                         -2,     0,     0,     0,   6*0, &
                                        -76,    25,   -11,     6,     1, &
                                          6,     0,     0,     0,   7*0, &
                                              -106,    18,    -6,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        6,     8,     3, &
                                          0,     0,     0,     0,   9*0, &
                                                               7,     0, &
                                          2,     0,     0,     0,  10*0, &
                                                                     -2/
      DATA (GY1D(I),I=1721,1800) / &
                                          3,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1940
      DATA (GY1D(I),I=1801,1945) /0, &
        -30654,  -1106,  1240,   914,  -241,    57,    74,    11,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2292,   2981, -1790,   762,   334,    54,   -53,     7,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1566,  1232,   550,   208,    -7,     4,    -3,     1, &
                                          2,     0,     0,     0,   4*0, &
                          916,  -405,   -33,  -249,    20,   -10,   -12, &
                                         -5,     0,     0,     0,   5*0, &
                                 265,  -141,   -18,   -31,     1,    11, &
                                         -2,     0,     0,     0,   6*0, &
                                        -76,    18,    -9,     6,     1, &
                                          6,     0,     0,     0,   7*0, &
                                              -107,    17,    -5,    -2, &
                                          4,     0,     0,     0,   8*0, &
                                                        5,     9,     3, &
                                          0,     0,     0,     0,   9*0, &
                                                               7,     1, &
                                          2,     0,     0,     0,  10*0, &
                                                                     -2/
      DATA (GY1D(I),I=1946,2025) / &
                                          3,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1945
      DATA (GY1D(I),I=2026,2170) /0, &
        -30594,  -1244,  1282,   944,  -253,    59,    70,    13,     5, &
                                         -3,     0,     0,     0,   2*0, &
         -2285,   2990, -1834,   776,   346,    57,   -40,     7,   -21, &
                                         11,     0,     0,     0,   3*0, &
                  1578,  1255,   544,   194,     6,     0,    -8,     1, &
                                          1,     0,     0,     0,   4*0, &
                          913,  -421,   -20,  -246,     0,    -5,   -11, &
                                          2,     0,     0,     0,   5*0, &
                                 304,  -142,   -25,   -29,     9,     3, &
                                         -5,     0,     0,     0,   6*0, &
                                        -82,    21,   -10,     7,    16, &
                                         -1,     0,     0,     0,   7*0, &
                                              -104,    15,   -10,    -3, &
                                          8,     0,     0,     0,   8*0, &
                                                       29,     7,    -4, &
                                         -1,     0,     0,     0,   9*0, &
                                                               2,    -3, &
                                         -3,     0,     0,     0,  10*0, &
                                                                     -4/
      DATA (GY1D(I),I=2171,2250) / &
                                          5,     0,     0,     0,  11*0, &
                                         -2,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1950
      DATA (GY1D(I),I=2251,2395) /0, &
        -30554,  -1341,  1297,   954,  -240,    54,    65,    22,     3, &
                                         -8,     0,     0,     0,   2*0, &
         -2250,   2998, -1889,   792,   349,    57,   -55,    15,    -7, &
                                          4,     0,     0,     0,   3*0, &
                  1576,  1274,   528,   211,     4,     2,    -4,    -1, &
                                         -1,     0,     0,     0,   4*0, &
                          896,  -408,   -20,  -247,     1,    -1,   -25, &
                                         13,     0,     0,     0,   5*0, &
                                 303,  -147,   -16,   -40,    11,    10, &
                                         -4,     0,     0,     0,   6*0, &
                                        -76,    12,    -7,    15,     5, &
                                          4,     0,     0,     0,   7*0, &
                                              -105,     5,   -13,    -5, &
                                         12,     0,     0,     0,   8*0, &
                                                       19,     5,    -2, &
                                          3,     0,     0,     0,   9*0, &
                                                              -1,     3, &
                                          2,     0,     0,     0,  10*0, &
                                                                      8/
      DATA (GY1D(I),I=2396,2475) / &
                                         10,     0,     0,     0,  11*0, &
                                          3,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1955
      DATA (GY1D(I),I=2476,2620) /0, &
        -30500,  -1440,  1302,   958,  -229,    47,    65,    11,     4, &
                                         -3,     0,     0,     0,   2*0, &
         -2215,   3003, -1944,   796,   360,    57,   -56,     9,     9, &
                                         -5,     0,     0,     0,   3*0, &
                  1581,  1288,   510,   230,     3,     2,    -6,    -4, &
                                         -1,     0,     0,     0,   4*0, &
                          882,  -397,   -23,  -247,    10,   -14,    -5, &
                                          2,     0,     0,     0,   5*0, &
                                 290,  -152,    -8,   -32,     6,     2, &
                                         -3,     0,     0,     0,   6*0, &
                                        -69,     7,   -11,    10,     4, &
                                          7,     0,     0,     0,   7*0, &
                                              -107,     9,    -7,     1, &
                                          4,     0,     0,     0,   8*0, &
                                                       18,     6,     2, &
                                         -2,     0,     0,     0,   9*0, &
                                                               9,     2, &
                                          6,     0,     0,     0,  10*0, &
                                                                      5/
      DATA (GY1D(I),I=2621,2700) / &
                                         -2,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1960
      DATA (GY1D(I),I=2701,2845) /0, &
        -30421,  -1555,  1302,   957,  -222,    46,    67,    15,     4, &
                                          1,     0,     0,     0,   2*0, &
         -2169,   3002, -1992,   800,   362,    58,   -56,     6,     6, &
                                         -3,     0,     0,     0,   3*0, &
                  1590,  1289,   504,   242,     1,     5,    -4,     0, &
                                          4,     0,     0,     0,   4*0, &
                          878,  -394,   -26,  -237,    15,   -11,    -9, &
                                          0,     0,     0,     0,   5*0, &
                                 269,  -156,    -1,   -32,     2,     1, &
                                         -1,     0,     0,     0,   6*0, &
                                        -63,    -2,    -7,    10,     4, &
                                          4,     0,     0,     0,   7*0, &
                                              -113,    17,    -5,    -1, &
                                          6,     0,     0,     0,   8*0, &
                                                        8,    10,    -2, &
                                          1,     0,     0,     0,   9*0, &
                                                               8,     3, &
                                         -1,     0,     0,     0,  10*0, &
                                                                     -1/
      DATA (GY1D(I),I=2846,2925) / &
                                          2,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1965
      DATA (GY1D(I),I=2926,3070) /0, &
        -30334,  -1662,  1297,   957,  -219,    45,    75,    13,     8, &
                                         -2,     0,     0,     0,   2*0, &
         -2119,   2997, -2038,   804,   358,    61,   -57,     5,    10, &
                                         -3,     0,     0,     0,   3*0, &
                  1594,  1292,   479,   254,     8,     4,    -4,     2, &
                                          2,     0,     0,     0,   4*0, &
                          856,  -390,   -31,  -228,    13,   -14,   -13, &
                                         -5,     0,     0,     0,   5*0, &
                                 252,  -157,     4,   -26,     0,    10, &
                                         -2,     0,     0,     0,   6*0, &
                                        -62,     1,    -6,     8,    -1, &
                                          4,     0,     0,     0,   7*0, &
                                              -111,    13,    -1,    -1, &
                                          4,     0,     0,     0,   8*0, &
                                                        1,    11,     5, &
                                          0,     0,     0,     0,   9*0, &
                                                               4,     1, &
                                          2,     0,     0,     0,  10*0, &
                                                                     -2/
      DATA (GY1D(I),I=3071,3150) / &
                                          2,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1970
      DATA (GY1D(I),I=3151,3295) /0, &
        -30220,  -1781,  1287,   952,  -216,    43,    72,    14,     8, &
                                         -3,     0,     0,     0,   2*0, &
         -2068,   3000, -2091,   800,   359,    64,   -57,     6,    10, &
                                         -3,     0,     0,     0,   3*0, &
                  1611,  1278,   461,   262,    15,     1,    -2,     2, &
                                          2,     0,     0,     0,   4*0, &
                          838,  -395,   -42,  -212,    14,   -13,   -12, &
                                         -5,     0,     0,     0,   5*0, &
                                 234,  -160,     2,   -22,    -3,    10, &
                                         -1,     0,     0,     0,   6*0, &
                                        -56,     3,    -2,     5,    -1, &
                                          6,     0,     0,     0,   7*0, &
                                              -112,    13,     0,     0, &
                                          4,     0,     0,     0,   8*0, &
                                                       -2,    11,     3, &
                                          1,     0,     0,     0,   9*0, &
                                                               3,     1, &
                                          0,     0,     0,     0,  10*0, &
                                                                     -1/
      DATA (GY1D(I),I=3296,3375) / &
                                          3,     0,     0,     0,  11*0, &
                                         -1,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1975
      DATA (GY1D(I),I=3376,3520) /0, &
        -30100,  -1902,  1276,   946,  -218,    45,    71,    14,     7, &
                                         -3,     0,     0,     0,   2*0, &
         -2013,   3010, -2144,   791,   356,    66,   -56,     6,    10, &
                                         -3,     0,     0,     0,   3*0, &
                  1632,  1260,   438,   264,    28,     1,    -1,     2, &
                                          2,     0,     0,     0,   4*0, &
                          830,  -405,   -59,  -198,    16,   -12,   -12, &
                                         -5,     0,     0,     0,   5*0, &
                                 216,  -159,     1,   -14,    -8,    10, &
                                         -2,     0,     0,     0,   6*0, &
                                        -49,     6,     0,     4,    -1, &
                                          5,     0,     0,     0,   7*0, &
                                              -111,    12,     0,    -1, &
                                          4,     0,     0,     0,   8*0, &
                                                       -5,    10,     4, &
                                          1,     0,     0,     0,   9*0, &
                                                               1,     1, &
                                          0,     0,     0,     0,  10*0, &
                                                                     -2/
      DATA (GY1D(I),I=3521,3600) / &
                                          3,     0,     0,     0,  11*0, &
                                         -1,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1980
      DATA (GY1D(I),I=3601,3745) /0, &
        -29992,  -1997,  1281,   938,  -218,    48,    72,    18,     5, &
                                         -4,     0,     0,     0,   2*0, &
         -1956,   3027, -2180,   782,   357,    66,   -59,     6,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1663,  1251,   398,   261,    42,     2,     0,     1, &
                                          2,     0,     0,     0,   4*0, &
                          833,  -419,   -74,  -192,    21,   -11,   -12, &
                                         -5,     0,     0,     0,   5*0, &
                                 199,  -162,     4,   -12,    -7,     9, &
                                         -2,     0,     0,     0,   6*0, &
                                        -48,    14,     1,     4,    -3, &
                                          5,     0,     0,     0,   7*0, &
                                              -108,    11,     3,    -1, &
                                          3,     0,     0,     0,   8*0, &
                                                       -2,     6,     7, &
                                          1,     0,     0,     0,   9*0, &
                                                              -1,     2, &
                                          2,     0,     0,     0,  10*0, &
                                                                     -5/
      DATA (GY1D(I),I=3746,3825) / &
                                          3,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1985
      DATA (GY1D(I),I=3826,3970) /0, &
        -29873,  -2072,  1296,   936,  -214,    53,    74,    21,     5, &
                                         -4,     0,     0,     0,   2*0, &
         -1905,   3044, -2208,   780,   355,    65,   -62,     6,    10, &
                                         -4,     0,     0,     0,   3*0, &
                  1687,  1247,   361,   253,    51,     3,     0,     1, &
                                          3,     0,     0,     0,   4*0, &
                          829,  -424,   -93,  -185,    24,   -11,   -12, &
                                         -5,     0,     0,     0,   5*0, &
                                 170,  -164,     4,    -6,    -9,     9, &
                                         -2,     0,     0,     0,   6*0, &
                                        -46,    16,     4,     4,    -3, &
                                          5,     0,     0,     0,   7*0, &
                                              -102,    10,     4,    -1, &
                                          3,     0,     0,     0,   8*0, &
                                                        0,     4,     7, &
                                          1,     0,     0,     0,   9*0, &
                                                              -4,     1, &
                                          2,     0,     0,     0,  10*0, &
                                                                     -5/
      DATA (GY1D(I),I=3971,4050) / &
                                          3,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1990
      DATA (GY1D(I),I=4051,4195) /0, &
        -29775,  -2131,  1314,   939,  -214,    61,    77,    23,     4, &
                                         -3,     0,     0,     0,   2*0, &
         -1848,   3059, -2239,   780,   353,    65,   -64,     5,     9, &
                                         -4,     0,     0,     0,   3*0, &
                  1686,  1248,   325,   245,    59,     2,    -1,     1, &
                                          2,     0,     0,     0,   4*0, &
                          802,  -423,  -109,  -178,    26,   -10,   -12, &
                                         -5,     0,     0,     0,   5*0, &
                                 141,  -165,     3,    -1,   -12,     9, &
                                         -2,     0,     0,     0,   6*0, &
                                        -36,    18,     5,     3,    -4, &
                                          4,     0,     0,     0,   7*0, &
                                               -96,     9,     4,    -2, &
                                          3,     0,     0,     0,   8*0, &
                                                        0,     2,     7, &
                                          1,     0,     0,     0,   9*0, &
                                                              -6,     1, &
                                          3,     0,     0,     0,  10*0, &
                                                                     -6/
      DATA (GY1D(I),I=4196,4275) / &
                                          3,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 1995
      DATA (GY1D(I),I=4276,4420) /0, &
        -29692,  -2200,  1335,   940,  -214,    68,    77,    25,     4, &
                                         -3,     0,     0,     0,   2*0, &
         -1784,   3070, -2267,   780,   352,    67,   -72,     6,     9, &
                                         -6,     0,     0,     0,   3*0, &
                  1681,  1249,   290,   235,    68,     1,    -6,     3, &
                                          2,     0,     0,     0,   4*0, &
                          759,  -418,  -118,  -170,    28,    -9,   -10, &
                                         -4,     0,     0,     0,   5*0, &
                                 122,  -166,    -1,     5,   -14,     8, &
                                         -1,     0,     0,     0,   6*0, &
                                        -17,    19,     4,     9,    -8, &
                                          4,     0,     0,     0,   7*0, &
                                               -93,     8,     6,    -1, &
                                          2,     0,     0,     0,   8*0, &
                                                       -2,    -5,    10, &
                                          2,     0,     0,     0,   9*0, &
                                                              -7,    -2, &
                                          5,     0,     0,     0,  10*0, &
                                                                     -8/
      DATA (GY1D(I),I=4421,4500) / &
                                          1,     0,     0,     0,  11*0, &
                                          0,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          g(n,m) for 2000
      DATA (GY1D(I),I=4501,4645) /0, &
      -29619.4,-2267.7,1339.6, 932.3,-218.8,  72.3,  79.0,  24.4,   5.0, &
                                       -2.6,   2.7,  -2.2,  -0.2,   2*0, &
       -1728.2, 3068.4,-2288.0, 786.8, 351.4,  68.2, -74.0,   6.6,  9.4, &
                                       -6.0,  -1.7,  -0.3,  -0.9,   3*0, &
                1670.9,1252.1, 250.0, 222.3,  74.2,   0.0,  -9.2,   3.0, &
                                        1.7,  -1.9,   0.2,   0.3,   4*0, &
                        714.5,-403.0,-130.4,-160.9,  33.3,  -7.9,  -8.4, &
                                       -3.1,   1.5,   0.9,   0.1,   5*0, &
                               111.3,-168.6,  -5.9,   9.1, -16.6,   6.3, &
                                       -0.5,  -0.1,  -0.2,  -0.4,   6*0, &
                                      -12.9,  16.9,   6.9,   9.1,  -8.9, &
                                        3.7,   0.1,   0.9,   1.3,   7*0, &
                                             -90.4,   7.3,   7.0,  -1.5, &
                                        1.0,  -0.7,  -0.5,  -0.4,   8*0, &
                                                     -1.2,  -7.9,   9.3, &
                                        2.0,   0.7,   0.3,   0.7,   9*0, &
                                                            -7.0,  -4.3, &
                                        4.2,   1.7,  -0.3,  -0.4,  10*0, &
                                                                  -8.2/
      DATA (GY1D(I),I=4646,4725) / &
                                        0.3,   0.1,  -0.4,   0.3,  11*0, &
                                       -1.1,   1.2,  -0.1,  -0.1,  12*0, &
                                               4.0,  -0.2,   0.4,  13*0, &
                                                     -0.4,   0.0,  14*0, &
                                                             0.1,  16*0/
!          g(n,m) for 2005
      DATA (GY1D(I),I=4726,4870) /0, &
      -29554.63,-2337.24,1336.30,920.55,-227.00,73.60,79.88,24.80, 5.58, &
                                          -2.17,2.95, -2.15,-0.16,  2*0, &
       -1669.05, 3047.69,-2305.83,797.96,354.41,69.56,-74.46,7.62, 9.76, &
                                          -6.12,-1.60, -0.29,-0.88, 3*0, &
                 1657.76, 1246.39,210.65,208.95,76.74,-1.65,-11.73,3.58, &
                                           1.42,-1.88,  0.21, 0.30, 4*0, &
                       672.51,-379.86,-136.54,-151.34,38.73,-6.88,-6.94, &
                                          -2.35, 1.44, 0.89,  0.28, 5*0, &
                               100.00,-168.05, -14.58,12.30,-18.11,5.01, &
                                          -0.15,-0.31,-0.38, -0.43, 6*0, &
                                      -13.55, 14.58, 9.37, 10.17,-10.76, &
                                           3.06, 0.29, 0.96,  1.18, 7*0, &
                                             -86.36, 5.42,  9.36, -1.25, &
                                           0.29,-0.79,-0.30,-0.37,  8*0, &
                                                     1.94,-11.25,  8.76, &
                                           2.06, 0.53, 0.46, 0.75,  9*0, &
                                                           -4.87, -6.66, &
                                           3.77, 1.80, -0.35,-0.26,10*0, &
                                                                  -9.22/
      DATA (GY1D(I),I=4871,4950) / &
                                          -0.21,0.16,  -0.36, 0.35,11*0, &
                                          -2.09,0.96,   0.08,-0.05,12*0, &
                                                3.99,  -0.49, 0.41,13*0, &
                                                       -0.08,-0.10,14*0, &
                                                             -0.18,16*0/
!          g(n,m) for 2010
      DATA (GY1D(I),I=4951,5095) /0, &
      -29496.57,-2396.06,1339.85,912.66,-230.87,72.78,80.44,24.41, 5.50, &
                                          -1.94, 3.05,-2.12,-0.09,  2*0, &
       -1586.42,3026.34,-2326.54,808.97,357.29,68.69,-75.00, 8.21, 9.45, &
                                          -6.24,-1.48,-0.21,-0.89,  3*0, &
                 1668.17,1232.10,166.58,200.26,75.92,-4.55,-14.50, 3.45, &
                                           0.89,-2.03, 0.30, 0.31,  4*0, &
                       633.73,-356.83,-141.05,-141.40,45.24,-5.59,-5.27, &
                                          -1.07, 1.65, 1.04, 0.42,  5*0, &
                                89.40,-163.17,-22.83,14.00,-19.34, 3.13, &
                                          -0.16,-0.51,-0.63,-0.45,  6*0, &
                                         -8.03,13.10,10.46,11.61,-12.38, &
                                           2.45, 0.54, 0.95, 1.08,  7*0, &
                                               -78.09, 1.64,10.85,-0.76, &
                                          -0.33,-0.79,-0.11,-0.31,  8*0, &
                                                      4.92,-14.05, 8.43, &
                                           2.13, 0.37, 0.52, 0.78,  9*0, &
                                                            -3.54,-8.42, &
                                           3.09, 1.79,-0.39,-0.18, 10*0, &
                                                                 -10.08/
      DATA (GY1D(I),I=5096,5175) / &
                                          -1.03, 0.12,-0.37, 0.38, 11*0, &
                                          -2.80, 0.75, 0.21, 0.02, 12*0, &
                                                 3.75,-0.77, 0.42, 13*0, &
                                                       0.04,-0.26, 14*0, &
                                                            -0.26, 16*0/
!          g(n,m) for 2015
  DATA (GY1D(I),I=5176,5320) /0, &
      -29441.46,-2445.88,1350.33,907.42,-232.91,69.55,81.29,23.98, 5.33, &
                                          -2.01, 3.00,-2.09,-0.02,  2*0, &
       -1501.77,3012.20,-2352.26,813.68,360.14,67.57,-75.99, 8.89, 8.83, &
                                          -6.26,-1.40,-0.16,-0.92,  3*0, &
                 1676.35,1225.85,120.49,192.35,72.79,-6.79,-16.78, 3.02, &
                                           0.17,-2.30, 0.46, 0.42,  4*0, &
                       581.69,-334.85,-140.94,-129.85,51.82,-3.16,-3.22, &
                                           0.55, 2.08, 1.23, 0.63,  5*0, &
                                70.38,-157.40,-28.93,15.07,-20.56, 0.67, &
                                          -0.55,-0.79,-0.89,-0.42,  6*0, &
                                          4.30,13.14, 9.32,13.33,-13.20, &
                                           1.70, 0.58, 0.85, 0.96,  7*0, &
                                               -70.85,-2.88,11.76,-0.10, &
                                          -0.67,-0.70, 0.10,-0.19,  8*0, &
                                                      6.61,-15.98, 8.68, &
                                           2.13, 0.14, 0.54, 0.81,  9*0, &
                                                            -2.02,-9.06, &
                                           2.33, 1.70,-0.37,-0.13, 10*0, &
                                                                 -10.54/
  DATA (GY1D(I),I=5321,5400) / &
                                          -1.80,-0.22,-0.43, 0.38, 11*0, &
                                          -3.59, 0.44, 0.22, 0.08, 12*0, &
                                                 3.49,-0.94, 0.46, 13*0, &
                                                      -0.03,-0.35, 14*0, &
                                                            -0.36, 16*0/
!          g(n,m) for 2020
  DATA (GY1D(I),I=5401,5545) /0, &
      -29403.41,-2499.78,1363.00,902.82,-234.42,65.97,80.54,23.66, 5.03, &
                                          -1.84, 2.96,-2.00, 0.08,  2*0, &
       -1451.37,2981.96,-2380.80,809.47,363.26,65.56,-76.63, 9.74, 8.36, &
                                          -6.25,-1.36,-0.13,-0.93,  3*0, &
                  1676.85,1236.06,86.18,187.86,72.96,-8.23,-17.49, 2.84, &
                                          -0.11,-2.51, 0.43, 0.53,  4*0, &
                       525.60,-309.47,-140.73,-121.57,56.45,-0.49,-1.48, &
                                           1.66, 2.31, 1.28, 0.72,  5*0, &
                                47.44,-151.16,-36.06,15.80,-21.07,-1.14, &
                                          -0.86,-0.85,-1.14,-0.30,  6*0, &
                                         13.98,13.60, 6.30,15.28,-13.22, &
                                           0.65, 0.28, 0.71, 0.75,  7*0, &
                                               -64.80,-7.21,13.65, 1.08, &
                                          -0.88,-0.66, 0.31,-0.01,  8*0, &
                                                      9.77,-16.59, 8.82, &
                                           1.88,-0.07, 0.49, 0.76,  9*0, &
                                                            -0.34,-9.23, &
                                           1.44, 1.44,-0.26,-0.05, 10*0, &
                                                                 -11.86/
  DATA (GY1D(I),I=5546,5625) / &
                                          -2.38,-0.59,-0.47, 0.37, 11*0, &
                                          -3.84, 0.18, 0.09, 0.13, 12*0, &
                                                 3.09,-1.13, 0.45, 13*0, &
                                                      -0.33,-0.46, 14*0, &
                                                            -0.40, 16*0/
!          g(n,m) for 2025
  DATA (GY1D(I),I=5626,5770) /0, &
       -29350.00,-2556.20,1360.90,894.70,-232.90,64.30,79.60,23.10, 4.70, &
                                           -1.30, 3.00,-2.00, 0.20,  2*0, &
        -1410.30,2950.90,-2404.20,799.60,369.00,63.80,-76.90,10.90, 8.00, &
                                           -6.40,-1.40,-0.10,-0.90,  3*0, &
                   1648.70,1243.80,55.80,187.20,76.70,-8.80,-17.50, 3.00, &
                                            0.20,-2.50, 0.40, 0.60,  4*0, &
                        453.40,-281.10,-138.70,-115.70,59.30, 2.00,-0.20, &
                                            2.00, 2.40, 1.20, 0.70,  5*0, &
                                 12.00,-141.90,-40.90,15.80,-21.80,-2.50, &
                                           -1.00,-0.60,-1.20,-0.20,  6*0, &
                                          20.90,14.90, 2.50,16.90,-13.10, &
                                           -0.50, 0.00, 0.60, 0.50,  7*0, &
                                               -60.80,-11.20,14.90, 2.40, &
                                           -0.90,-0.60, 0.50, 0.10,  8*0, &
                                                      14.30,-16.80, 8.60, &
                                            1.50,-0.10, 0.50, 0.70,  9*0, &
                                                              1.00,-8.70, &
                                            0.90, 1.10,-0.10, 0.00, 10*0, &
                                                                  -12.80/
  DATA (GY1D(I),I=5771,5850) / &
                                           -2.60,-1.00,-0.50, 0.30, 11*0, &
                                           -3.90,-0.10,-0.20, 0.20, 12*0, &
                                                  2.60,-1.20, 0.40, 13*0, &
                                                       -0.70,-0.50, 14*0, &
                                                             -0.40, 16*0/

!          h(n,m) for 1900
      DATA (HY1D(I),I=1,145) /16*0, &
          5922,  -1061,  -330,   195,  -210,    -9,   -45,     8,   -20, &
                                          2,     0,     0,     0,   3*0, &
                  1121,     3,   -69,    53,    83,   -13,   -14,    14, &
                                          1,     0,     0,     0,   4*0, &
                          523,  -210,   -33,     2,   -10,     7,     5, &
                                          2,     0,     0,     0,   5*0, &
                                 -75,  -124,   -35,    -1,   -13,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                          3,    36,    28,     5,    -2, &
                                         -4,     0,     0,     0,   7*0, &
                                               -69,   -12,    16,     8, &
                                          0,     0,     0,     0,   8*0, &
                                                      -22,    -5,    10, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -18,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=146,225) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1905
      DATA (HY1D(I),I=226,370) /16*0, &
          5909,  -1086,  -357,   203,  -193,    -7,   -46,     8,   -20, &
                                          2,     0,     0,     0,   3*0, &
                  1065,    34,   -77,    56,    86,   -14,   -15,    14, &
                                          1,     0,     0,     0,   4*0, &
                          480,  -201,   -32,     4,   -11,     7,     5, &
                                          2,     0,     0,     0,   5*0, &
                                 -65,  -125,   -32,     0,   -13,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                         11,    32,    28,     5,    -2, &
                                         -4,     0,     0,     0,   7*0, &
                                               -67,   -12,    16,     8, &
                                          0,     0,     0,     0,   8*0, &
                                                      -22,    -5,    10, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -18,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=371,450) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1910
      DATA (HY1D(I),I=451,595) /16*0, &
          5898,  -1128,  -389,   211,  -172,    -5,   -47,     8,   -20, &
                                          2,     0,     0,     0,   3*0, &
                  1000,    62,   -90,    57,    89,   -14,   -15,    14, &
                                          1,     0,     0,     0,   4*0, &
                          425,  -189,   -33,     5,   -12,     6,     5, &
                                          2,     0,     0,     0,   5*0, &
                                 -55,  -126,   -29,     1,   -13,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                         21,    28,    28,     5,    -2, &
                                         -4,     0,     0,     0,   7*0, &
                                               -65,   -13,    16,     8, &
                                          0,     0,     0,     0,   8*0, &
                                                      -22,    -5,    10, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -18,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=596,675) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1915
      DATA (HY1D(I),I=676,820) /16*0, &
          5875,  -1191,  -421,   218,  -148,    -2,   -48,     8,   -20, &
                                          2,     0,     0,     0,   3*0, &
                   917,    84,  -109,    58,    93,   -14,   -15,    14, &
                                          1,     0,     0,     0,   4*0, &
                          360,  -173,   -34,     8,   -12,     6,     5, &
                                          2,     0,     0,     0,   5*0, &
                                 -51,  -126,   -26,     2,   -13,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                         32,    23,    28,     5,    -2, &
                                         -4,     0,     0,     0,   7*0, &
                                               -62,   -15,    16,     8, &
                                          0,     0,     0,     0,   8*0, &
                                                      -22,    -5,    10, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -18,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=821,900) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1920
      DATA (HY1D(I),I=901,1045) /16*0, &
          5845,  -1259,  -445,   220,  -122,     0,   -49,     8,   -20, &
                                          2,     0,     0,     0,   3*0, &
                   823,   103,  -134,    58,    96,   -14,   -15,    14, &
                                          1,     0,     0,     0,   4*0, &
                          293,  -153,   -38,    11,   -13,     6,     5, &
                                          2,     0,     0,     0,   5*0, &
                                 -57,  -125,   -22,     4,   -14,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                         43,    18,    28,     5,    -2, &
                                         -4,     0,     0,     0,   7*0, &
                                               -57,   -16,    17,     9, &
                                          0,     0,     0,     0,   8*0, &
                                                      -22,    -5,    10, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -19,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=1046,1125) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1925
      DATA (HY1D(I),I=1126,1270) /16*0, &
          5817,  -1334,  -462,   216,   -96,     3,   -50,     8,   -20, &
                                          2,     0,     0,     0,   3*0, &
                   728,   119,  -163,    58,    99,   -14,   -15,    14, &
                                          1,     0,     0,     0,   4*0, &
                          229,  -130,   -44,    14,   -14,     6,     5, &
                                          2,     0,     0,     0,   5*0, &
                                 -70,  -122,   -18,     5,   -14,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                         51,    13,    29,     5,    -2, &
                                         -4,     0,     0,     0,   7*0, &
                                               -52,   -17,    17,     9, &
                                          0,     0,     0,     0,   8*0, &
                                                      -21,    -5,    10, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -19,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=1271,1350) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1930
      DATA (HY1D(I),I=1351,1495) /16*0, &
          5808,  -1424,  -480,   205,   -72,     4,   -51,     8,   -20, &
                                          2,     0,     0,     0,   3*0, &
                   644,   133,  -195,    60,   102,   -15,   -15,    14, &
                                          1,     0,     0,     0,   4*0, &
                          166,  -109,   -53,    19,   -14,     5,     5, &
                                          2,     0,     0,     0,   5*0, &
                                 -90,  -118,   -16,     6,   -14,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                         58,     8,    29,     5,    -2, &
                                         -4,     0,     0,     0,   7*0, &
                                               -46,   -18,    18,     9, &
                                          0,     0,     0,     0,   8*0, &
                                                      -20,    -5,    10, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -19,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=1496,1575) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1935
      DATA (HY1D(I),I=1576,1720) /16*0, &
          5812,  -1520,  -494,   188,   -51,     4,   -52,     8,   -20, &
                                          2,     0,     0,     0,   3*0, &
                   586,   146,  -226,    64,   104,   -17,   -15,    15, &
                                          1,     0,     0,     0,   4*0, &
                          101,   -90,   -64,    25,   -14,     5,     5, &
                                          2,     0,     0,     0,   5*0, &
                                -114,  -115,   -15,     7,   -15,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                         64,     4,    29,     5,    -3, &
                                         -4,     0,     0,     0,   7*0, &
                                               -40,   -19,    18,     9, &
                                          0,     0,     0,     0,   8*0, &
                                                      -19,    -5,    11, &
                                         -1,     0,     0,     0,   9*0, &
                                                             -19,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=1721,1800) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1940
      DATA (HY1D(I),I=1801,1945) /16*0, &
          5821,  -1614,  -499,   169,   -33,     4,   -52,     8,   -21, &
                                          2,     0,     0,     0,   3*0, &
                   528,   163,  -252,    71,   105,   -18,   -14,    15, &
                                          1,     0,     0,     0,   4*0, &
                           43,   -72,   -75,    33,   -14,     5,     5, &
                                          2,     0,     0,     0,   5*0, &
                                -141,  -113,   -15,     7,   -15,    -3, &
                                          6,     0,     0,     0,   6*0, &
                                         69,     0,    29,     5,    -3, &
                                         -4,     0,     0,     0,   7*0, &
                                               -33,   -20,    19,     9, &
                                          0,     0,     0,     0,   8*0, &
                                                      -19,    -5,    11, &
                                         -1,     0,     0,     0,   9*0, &
                                                             -19,    -2, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=1946,2025) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1945
      DATA (HY1D(I),I=2026,2170) /16*0, &
          5810,  -1702,  -499,   144,   -12,     6,   -45,    12,   -27, &
                                          5,     0,     0,     0,   3*0, &
                   477,   186,  -276,    95,   100,   -18,   -21,    17, &
                                          1,     0,     0,     0,   4*0, &
                          -11,   -55,   -67,    16,     2,   -12,    29, &
                                        -20,     0,     0,     0,   5*0, &
                                -178,  -119,    -9,     6,    -7,    -9, &
                                         -1,     0,     0,     0,   6*0, &
                                         82,   -16,    28,     2,     4, &
                                         -6,     0,     0,     0,   7*0, &
                                               -39,   -17,    18,     9, &
                                          6,     0,     0,     0,   8*0, &
                                                      -22,     3,     6, &
                                         -4,     0,     0,     0,   9*0, &
                                                             -11,     1, &
                                         -2,     0,     0,     0,  10*0, &
                                                                      8/
      DATA (HY1D(I),I=2171,2250) / &
                                          0,     0,     0,     0,  11*0, &
                                         -2,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1950
      DATA (HY1D(I),I=2251,2395) /16*0, &
          5815,  -1810,  -476,   136,     3,    -1,   -35,     5,   -24, &
                                         13,     0,     0,     0,   3*0, &
                   381,   206,  -278,   103,    99,   -17,   -22,    19, &
                                         -2,     0,     0,     0,   4*0, &
                          -46,   -37,   -87,    33,     0,     0,    12, &
                                        -10,     0,     0,     0,   5*0, &
                                -210,  -122,   -12,    10,   -21,     2, &
                                          2,     0,     0,     0,   6*0, &
                                         80,   -12,    36,    -8,     2, &
                                         -3,     0,     0,     0,   7*0, &
                                               -30,   -18,    17,     8, &
                                          6,     0,     0,     0,   8*0, &
                                                      -16,    -4,     8, &
                                         -3,     0,     0,     0,   9*0, &
                                                             -17,   -11, &
                                          6,     0,     0,     0,  10*0, &
                                                                     -7/
      DATA (HY1D(I),I=2396,2475) / &
                                         11,     0,     0,     0,  11*0, &
                                          8,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1955
      DATA (HY1D(I),I=2476,2620) /16*0, &
          5820,  -1898,  -462,   133,    15,    -9,   -50,    10,   -11, &
                                         -4,     0,     0,     0,   3*0, &
                   291,   216,  -274,   110,    96,   -24,   -15,    12, &
                                          0,     0,     0,     0,   4*0, &
                          -83,   -23,   -98,    48,    -4,     5,     7, &
                                         -8,     0,     0,     0,   5*0, &
                                -230,  -121,   -16,     8,   -23,     6, &
                                         -2,     0,     0,     0,   6*0, &
                                         78,   -12,    28,     3,    -2, &
                                         -4,     0,     0,     0,   7*0, &
                                               -24,   -20,    23,    10, &
                                          1,     0,     0,     0,   8*0, &
                                                      -18,    -4,     7, &
                                         -3,     0,     0,     0,   9*0, &
                                                             -13,    -6, &
                                          7,     0,     0,     0,  10*0, &
                                                                      5/
      DATA (HY1D(I),I=2621,2700) / &
                                         -1,     0,     0,     0,  11*0, &
                                         -3,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1960
      DATA (HY1D(I),I=2701,2845) /16*0, &
          5791,  -1967,  -414,   135,    16,   -10,   -55,    11,   -18, &
                                          4,     0,     0,     0,   3*0, &
                   206,   224,  -278,   125,    99,   -28,   -14,    12, &
                                          1,     0,     0,     0,   4*0, &
                         -130,     3,  -117,    60,    -6,     7,     2, &
                                          0,     0,     0,     0,   5*0, &
                                -255,  -114,   -20,     7,   -18,     0, &
                                          2,     0,     0,     0,   6*0, &
                                         81,   -11,    23,     4,    -3, &
                                         -5,     0,     0,     0,   7*0, &
                                               -17,   -18,    23,     9, &
                                          1,     0,     0,     0,   8*0, &
                                                      -17,     1,     8, &
                                         -1,     0,     0,     0,   9*0, &
                                                             -20,     0, &
                                          6,     0,     0,     0,  10*0, &
                                                                      5/
      DATA (HY1D(I),I=2846,2925) / &
                                          0,     0,     0,     0,  11*0, &
                                         -7,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1965
      DATA (HY1D(I),I=2926,3070) /16*0, &
          5776,  -2016,  -404,   148,    19,   -11,   -61,     7,   -22, &
                                          2,     0,     0,     0,   3*0, &
                   114,   240,  -269,   128,   100,   -27,   -12,    15, &
                                          1,     0,     0,     0,   4*0, &
                         -165,    13,  -126,    68,    -2,     9,     7, &
                                          2,     0,     0,     0,   5*0, &
                                -269,   -97,   -32,     6,   -16,    -4, &
                                          6,     0,     0,     0,   6*0, &
                                         81,    -8,    26,     4,    -5, &
                                         -4,     0,     0,     0,   7*0, &
                                                -7,   -23,    24,    10, &
                                          0,     0,     0,     0,   8*0, &
                                                      -12,    -3,    10, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -17,    -4, &
                                          3,     0,     0,     0,  10*0, &
                                                                      1/
      DATA (HY1D(I),I=3071,3150) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1970
      DATA (HY1D(I),I=3151,3295) /16*0, &
          5737,  -2047,  -366,   167,    26,   -12,   -70,     7,   -21, &
                                          1,     0,     0,     0,   3*0, &
                    25,   251,  -266,   139,   100,   -27,   -15,    16, &
                                          1,     0,     0,     0,   4*0, &
                         -196,    26,  -139,    72,    -4,     6,     6, &
                                          3,     0,     0,     0,   5*0, &
                                -279,   -91,   -37,     8,   -17,    -4, &
                                          4,     0,     0,     0,   6*0, &
                                         83,    -6,    23,     6,    -5, &
                                         -4,     0,     0,     0,   7*0, &
                                                 1,   -23,    21,    10, &
                                          0,     0,     0,     0,   8*0, &
                                                      -11,    -6,    11, &
                                         -1,     0,     0,     0,   9*0, &
                                                             -16,    -2, &
                                          3,     0,     0,     0,  10*0, &
                                                                      1/
      DATA (HY1D(I),I=3296,3375) / &
                                          1,     0,     0,     0,  11*0, &
                                         -4,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1975
      DATA (HY1D(I),I=3376,3520) /16*0, &
          5675,  -2067,  -333,   191,    31,   -13,   -77,     6,   -21, &
                                          1,     0,     0,     0,   3*0, &
                   -68,   262,  -265,   148,    99,   -26,   -16,    16, &
                                          1,     0,     0,     0,   4*0, &
                         -223,    39,  -152,    75,    -5,     4,     7, &
                                          3,     0,     0,     0,   5*0, &
                                -288,   -83,   -41,    10,   -19,    -4, &
                                          4,     0,     0,     0,   6*0, &
                                         88,    -4,    22,     6,    -5, &
                                         -4,     0,     0,     0,   7*0, &
                                                11,   -23,    18,    10, &
                                         -1,     0,     0,     0,   8*0, &
                                                      -12,   -10,    11, &
                                         -1,     0,     0,     0,   9*0, &
                                                             -17,    -3, &
                                          3,     0,     0,     0,  10*0, &
                                                                      1/
      DATA (HY1D(I),I=3521,3600) / &
                                          1,     0,     0,     0,  11*0, &
                                         -5,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1980
      DATA (HY1D(I),I=3601,3745) /16*0, &
          5604,  -2129,  -336,   212,    46,   -15,   -82,     7,   -21, &
                                          1,     0,     0,     0,   3*0, &
                  -200,   271,  -257,   150,    93,   -27,   -18,    16, &
                                          0,     0,     0,     0,   4*0, &
                         -252,    53,  -151,    71,    -5,     4,     9, &
                                          3,     0,     0,     0,   5*0, &
                                -297,   -78,   -43,    16,   -22,    -5, &
                                          6,     0,     0,     0,   6*0, &
                                         92,    -2,    18,     9,    -6, &
                                         -4,     0,     0,     0,   7*0, &
                                                17,   -23,    16,     9, &
                                          0,     0,     0,     0,   8*0, &
                                                      -10,   -13,    10, &
                                         -1,     0,     0,     0,   9*0, &
                                                             -15,    -6, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=3746,3825) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1985
      DATA (HY1D(I),I=3826,3970) /16*0, &
          5500,  -2197,  -310,   232,    47,   -16,   -83,     8,   -21, &
                                          1,     0,     0,     0,   3*0, &
                  -306,   284,  -249,   150,    88,   -27,   -19,    15, &
                                          0,     0,     0,     0,   4*0, &
                         -297,    69,  -154,    69,    -2,     5,     9, &
                                          3,     0,     0,     0,   5*0, &
                                -297,   -75,   -48,    20,   -23,    -6, &
                                          6,     0,     0,     0,   6*0, &
                                         95,    -1,    17,    11,    -6, &
                                         -4,     0,     0,     0,   7*0, &
                                                21,   -23,    14,     9, &
                                          0,     0,     0,     0,   8*0, &
                                                       -7,   -15,     9, &
                                         -1,     0,     0,     0,   9*0, &
                                                             -11,    -7, &
                                          4,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=3971,4050) / &
                                          0,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1990
      DATA (HY1D(I),I=4051,4195) /16*0, &
          5406,  -2279,  -284,   247,    46,   -16,   -80,    10,   -20, &
                                          2,     0,     0,     0,   3*0, &
                  -373,   293,  -240,   154,    82,   -26,   -19,    15, &
                                          1,     0,     0,     0,   4*0, &
                         -352,    84,  -153,    69,     0,     6,    11, &
                                          3,     0,     0,     0,   5*0, &
                                -299,   -69,   -52,    21,   -22,    -7, &
                                          6,     0,     0,     0,   6*0, &
                                         97,     1,    17,    12,    -7, &
                                         -4,     0,     0,     0,   7*0, &
                                                24,   -23,    12,     9, &
                                          0,     0,     0,     0,   8*0, &
                                                       -4,   -16,     8, &
                                         -2,     0,     0,     0,   9*0, &
                                                             -10,    -7, &
                                          3,     0,     0,     0,  10*0, &
                                                                      2/
      DATA (HY1D(I),I=4196,4275) / &
                                         -1,     0,     0,     0,  11*0, &
                                         -6,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 1995
      DATA (HY1D(I),I=4276,4420) /16*0, &
          5306,  -2366,  -262,   262,    46,   -17,   -69,    11,   -20, &
                                          1,     0,     0,     0,   3*0, &
                  -413,   302,  -236,   165,    72,   -25,   -21,    15, &
                                          0,     0,     0,     0,   4*0, &
                         -427,    97,  -143,    67,     4,     8,    12, &
                                          4,     0,     0,     0,   5*0, &
                                -306,   -55,   -58,    24,   -23,    -6, &
                                          5,     0,     0,     0,   6*0, &
                                        107,     1,    17,    15,    -8, &
                                         -5,     0,     0,     0,   7*0, &
                                                36,   -24,    11,     8, &
                                         -1,     0,     0,     0,   8*0, &
                                                       -6,   -16,     5, &
                                         -2,     0,     0,     0,   9*0, &
                                                              -4,    -8, &
                                          1,     0,     0,     0,  10*0, &
                                                                      3/
      DATA (HY1D(I),I=4421,4500) / &
                                         -2,     0,     0,     0,  11*0, &
                                         -7,     0,     0,     0,  12*0, &
                                                 0,     0,     0,  13*0, &
                                                        0,     0,  14*0, &
                                                               0,  16*0/
!          h(n,m) for 2000
      DATA (HY1D(I),I=4501,4645) /16*0, &
        5186.1,-2481.6,-227.6, 272.6,  43.8, -17.4, -64.6,  11.9, -19.7, &
                                        1.7,   0.1,  -0.4,  -0.9,   3*0, &
                -458.0, 293.4,-231.9, 171.9,  63.7, -24.2, -21.5,  13.4, &
                                        0.0,   1.3,   0.3,   0.2,   4*0, &
                       -491.1, 119.8,-133.1,  65.1,   6.2,   8.5,  12.5, &
                                        4.0,  -0.9,   2.5,   1.8,   5*0, &
                              -303.8, -39.3, -61.2,  24.0, -21.5,  -6.2, &
                                        4.9,  -2.6,  -2.6,  -0.4,   6*0, &
                                      106.3,   0.7,  14.8,  15.5,  -8.4, &
                                       -5.9,   0.9,   0.7,  -1.0,   7*0, &
                                              43.8, -25.4,   8.9,   8.4, &
                                       -1.2,  -0.7,   0.3,  -0.1,   8*0, &
                                                     -5.8, -14.9,   3.8, &
                                       -2.9,  -2.8,   0.0,   0.7,   9*0, &
                                                            -2.1,  -8.2, &
                                        0.2,  -0.9,   0.0,   0.3,  10*0, &
                                                                    4.8/
      DATA (HY1D(I),I=4646,4725) / &
                                       -2.2,  -1.2,   0.3,   0.6,  11*0, &
                                       -7.4,  -1.9,  -0.9,   0.3,  12*0, &
                                              -0.9,  -0.4,  -0.2,  13*0, &
                                                      0.8,  -0.5,  14*0, &
                                                             0.1,  16*0/
!          h(n,m) for 2005
      DATA (HY1D(I),I=4726,4870) /16*0, &
       5077.99,-2594.50,-198.86,282.07,42.72,-20.33,-61.14,11.20,-20.11, &
                                        2.19,  0.26, -0.55,-0.76,   3*0, &
               -515.43, 269.72,-225.23,180.25,54.75,-22.57,-20.88,12.69, &
                                        0.10,  1.44,  0.23, 0.33,   4*0, &
                       -524.72, 145.15,-123.45,63.63, 6.82, 9.83, 12.67, &
                                        4.46, -0.77,  2.38, 1.72,   5*0, &
                               -305.36,-19.57,-63.53,25.35,-19.71,-6.72, &
                                        4.76, -2.27, -2.63, -0.54,  6*0, &
                                       103.85,  0.24, 10.93,16.22,-8.16, &
                                       -6.58,  0.90,  0.61, -1.07,  7*0, &
                                              50.94,-26.32,  7.61, 8.10, &
                                       -1.01, -0.58,  0.40, -0.04,  8*0, &
                                                     -4.64,-12.76, 2.92, &
                                       -3.47, -2.69,  0.01,  0.63,  9*0, &
                                                            -0.06,-7.73, &
                                       -0.86, -1.08,  0.02,  0.21, 10*0, &
                                                                   6.01/
      DATA (HY1D(I),I=4871,4950) / &
                                       -2.31, -1.58,  0.28,  0.53, 11*0, &
                                       -7.93, -1.90, -0.87,  0.38, 12*0, &
                                              -1.39, -0.34, -0.22, 13*0, &
                                                      0.88, -0.57, 14*0, &
                                                           -0.82, 16*0/
!          h(n,m) for 2010
      DATA (HY1D(I),I=4951,5095) /16*0, &
       4944.26,-2708.54,-160.40,286.48,44.58,-20.90,-57.80,10.84,-20.54, &
                                           2.73, 0.13,-0.87,-0.87,  3*0, &
                -575.73,251.75,-211.03,189.01,44.18,-21.20,-20.03,11.51, &
                                          -0.10, 1.67, 0.27, 0.30,  4*0, &
                         -537.03,164.46,-118.06,61.54, 6.54,11.83,12.75, &
                                           4.71,-0.66, 2.13, 1.66,  5*0, &
                                -309.72,-0.01,-66.26,24.96,-17.41,-7.14, &
                                           4.44,-1.76,-2.49,-0.59,  6*0, &
                                         101.04, 3.02, 7.03,16.71,-7.42, &
                                          -7.22, 0.85, 0.49,-1.14,  7*0, &
                                               55.40,-27.61, 6.96, 7.97, &
                                          -0.96,-0.39, 0.59,-0.07,  8*0, &
                                                     -3.28,-10.74, 2.14, &
                                          -3.95,-2.51, 0.00, 0.54,  9*0, &
                                                             1.64,-6.08, &
                                          -1.99,-1.27, 0.13, 0.10, 10*0, &
                                                                   7.01/
      DATA (HY1D(I),I=5096,5175) / &
                                          -1.97,-2.11, 0.27, 0.49, 11*0, &
                                          -8.31,-1.94,-0.86, 0.44, 12*0, &
                                                -1.86,-0.23,-0.25, 13*0, &
                                                       0.87,-0.53, 14*0, &
                                                            -0.79, 16*0/
!          h(n,m) for 2015
  DATA (HY1D(I),I=5176,5320) /16*0, &
       4795.99,-2845.41,-115.29,283.54,46.98,-20.61,-54.27,10.04,-21.77, &
                                           3.28, 0.00,-1.08,-0.88,  3*0, &
                -642.17,245.04,-188.43,196.98,33.30,-19.53,-18.26,10.76, &
                                          -0.40, 2.11, 0.37, 0.49,  4*0, &
                         -538.70,180.95,-119.14,58.74, 5.59,13.18,11.74, &
                                           4.55,-0.60, 1.75, 1.56,  5*0, &
                                -329.23,15.98,-66.64,24.45,-14.60,-6.74, &
                                           4.40,-1.05,-2.19,-0.50,  6*0, &
                                         100.12, 7.35, 3.27,16.16,-6.88, &
                                          -7.92, 0.76, 0.27,-1.24,  7*0, &
                                               62.41,-27.50, 5.69, 7.79, &
                                          -0.61,-0.20, 0.72,-0.10,  8*0, &
                                                      -2.32,-9.10, 1.04, &
                                          -4.16,-2.12,-0.09, 0.42,  9*0, &
                                                             2.26,-3.89, &
                                          -2.85,-1.44, 0.29,-0.04, 10*0, &
                                                                   8.44/
  DATA (HY1D(I),I=5321,5400) / &
                                          -1.12,-2.57, 0.23, 0.48, 11*0, &
                                          -8.72,-2.01,-0.89, 0.48, 12*0, &
                                                -2.34,-0.16,-0.30, 13*0, &
                                                       0.72,-0.43, 14*0, &
                                                            -0.71, 16*0/
!          h(n,m) for 2020
  DATA (HY1D(I),I=5401,5545) /16*0, &
        4653.35,-2991.72,-81.96,282.10,47.52,-19.22,-51.50, 8.43,-23.44, &
                                           3.38,-0.02,-1.15,-0.88,  3*0, &
                -734.62,241.80,-158.50,208.36,25.02,-16.85,-15.23,11.04, &
                                          -0.18, 2.50, 0.52, 0.64,  4*0, &
                         -542.52,199.75,-121.43,52.76, 2.36,12.83, 9.86, &
                                           3.50,-0.55, 1.37, 1.40,  5*0, &
                                -350.30,32.09,-64.40,23.56,-11.76,-5.13, &
                                           4.86,-0.39,-1.81,-0.38,  6*0, &
                                          99.14, 8.96,-2.19,14.94,-6.20, &
                                          -8.62, 0.62, 0.08,-1.31,  7*0, &
                                               68.04,-27.19, 3.62, 7.79, &
                                          -0.11,-0.21, 0.71,-0.09,  8*0, &
                                                      -1.90,-6.90, 0.40, &
                                          -4.26,-1.66,-0.15, 0.29,  9*0, &
                                                             2.90,-1.44, &
                                          -3.43,-1.60, 0.55,-0.11, 10*0, &
                                                                   9.60/
  DATA (HY1D(I),I=5546,5625) / &
                                          -0.10,-2.98, 0.16, 0.47, 11*0, &
                                          -8.84,-1.97,-0.93, 0.54, 12*0, &
                                                -2.51,-0.04,-0.41, 13*0, &
                                                       0.52,-0.36, 14*0, &
                                                            -0.60, 16*0/
!          h(n,m) for 2025
  DATA (HY1D(I),I=5626,5770) /16*0, &
        4545.50,-3133.60,-56.90,278.60,45.30,-18.40,-48.90, 7.20,-24.80, &
                                           3.30, 0.00,-1.20,-0.90,  3*0, &
                -814.20,237.60,-134.00,220.00,16.80,-14.40,-12.60,12.10, &
                                           0.10, 2.80, 0.60, 0.70,  4*0, &
                         -549.60,212.00,-122.90,48.90,-1.00,11.50, 8.30, &
                                           2.50,-0.60, 1.00, 1.20,  5*0, &
                                 -375.40,42.90,-59.80,23.50,-9.70,-3.40, &
                                           5.40, 0.10,-1.50,-0.30,  6*0, &
                                         106.20,10.90,-7.40,12.70,-5.30, &
                                          -9.00, 0.50, 0.00,-1.30,  7*0, &
                                               72.80,-25.10, 0.70, 7.20, &
                                           0.40,-0.30, 0.60,-0.10,  8*0, &
                                                      -2.20,-5.20,-0.60, &
                                          -4.20,-1.20,-0.20, 0.20,  9*0, &
                                                             3.90, 0.80, &
                                          -3.80,-1.70, 0.80,-0.20, 10*0, &
                                                                   9.80/
  DATA (HY1D(I),I=5771,5850) / &
                                           0.90,-2.90, 0.10, 0.50, 11*0, &
                                          -9.00,-1.80,-0.90, 0.60, 12*0, &
                                                -2.30, 0.10,-0.60, 13*0, &
                                                       0.20,-0.30, 14*0, &
                                                            -0.50, 16*0/

!          Secular variation rates are nominally okay through 2030
      DATA (GT1D(I),I=1,145) /0, &
                 12.60,-11.20,-1.50,-1.70, 0.60,-0.20,-0.10,-0.10, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  2*0, &
                  10.00,-5.30,-4.40,-2.30, 1.30,-0.30,-0.10, 0.20, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  3*0, &
                        -8.30, 0.40,-5.80, 0.00, 0.80,-0.10, 0.00, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  4*0, &
                             -15.60, 5.40, 0.70, 1.20, 0.50, 0.40, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  5*0, &
                                    -6.80, 2.30,-0.80,-0.10,-0.10, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  6*0, &
                                           1.00, 0.40,-0.80, 0.30, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  7*0, &
                                                 0.90,-0.80, 0.10, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  8*0, &
                                                       0.90, 0.00, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  9*0, &
                                                             0.30, 0.00, &
                                           0.00, 0.00, 0.00, 0.00, 10*0, &
                                                                   0.00/
      DATA (GT1D(I),I=146,225) / &
                                           0.00, 0.00, 0.00, 0.00, 11*0, &
                                           0.00, 0.00, 0.00, 0.00, 12*0, &
                                                 0.00, 0.00, 0.00, 13*0, &
                                                       0.00, 0.00, 14*0, &
                                                             0.00, 16*0/

      DATA (HT1D(I),I=1,145) /16*0, &
                -21.50,-27.30, 3.80,-1.30,-0.50, 0.30, 0.60,-0.30, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  3*0, &
                       -11.10,-0.20, 4.10, 2.10,-1.60, 0.50, 0.40, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  4*0, &
                              -3.90, 1.60, 0.50,-0.40,-0.70,-0.30, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  5*0, &
                                    -4.10, 1.70, 0.80, 0.00, 0.40, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  6*0, &
                                           1.90, 0.70,-0.90,-0.50, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  7*0, &
                                                 0.90, 0.50,-0.60, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  8*0, &
                                                      -0.30, 0.30, 0.00, &
                                           0.00, 0.00, 0.00, 0.00,  9*0, &
                                                             0.20, 0.00, &
                                           0.00, 0.00, 0.00, 0.00, 10*0, &
                                                                   0.00/
      DATA (HT1D(I),I=146,225) / &
                                           0.00, 0.00, 0.00, 0.00, 11*0, &
                                           0.00, 0.00, 0.00, 0.00, 12*0, &
                                                 0.00, 0.00, 0.00, 13*0, &
                                                       0.00, 0.00, 14*0, &
                                                             0.00, 16*0/

!          Do not need to load new coefficients if date has not changed
      ICHG = 0
      IF (DATE .EQ. DATEL) GO TO 300
      DATEL = DATE
      ICHG = 1

!          Trap out of range date:
      IF (DATE .LT. EPOCH(1)) GO TO 9100
      IF (DATE .GT. EPOCH(NEPO)+5.) WRITE(0,9200) DATE, EPOCH(NEPO) + 5.

      DO 100 I=1,NEPO
      IF (DATE .LT. EPOCH(I)) GO TO 110
      IY = I
  100 CONTINUE
  110 CONTINUE

      NMAX  = NMXE(IY)
      TIME  = DATE
      T     = TIME-EPOCH(IY)
      TO5   = T/5.
      IY1   = IY + 1
      GB(1) = 0.0
      GV(1) = 0.0
      I  = 2
      F0 = -1.0D-5
      DO 200 N=1,NMAX
      F0 = F0 * REAL(N)/2.
      F  = F0 / SQRT(2.0)
      NN = N+1
      MM = 1
      IF (IY .LT. NEPO) GB(I) = (GYR(NN,MM,IY) + &
                                (GYR(NN,MM,IY1)-GYR(NN,MM,IY))*TO5) * F0
      IF (IY .EQ. NEPO) GB(I) = (GYR(NN,MM,IY) + GT(NN,MM)    *T  ) * F0
      GV(I) = GB(I) / REAL(NN)
      I = I+1
      DO 200 M=1,N
      F  = F / SQRT( REAL(N-M+1) / REAL(N+M) )
      NN = N+1
      MM = M+1
      I1 = I+1
      IF (IY .LT. NEPO) THEN
      GB(I)  = (GYR(NN,MM,IY) + &
                 (GYR(NN,MM,IY1)-GYR(NN,MM,IY))*TO5) * F
      GB(I1) = (HYR(NN,MM,IY) + &
                 (HYR(NN,MM,IY1)-HYR(NN,MM,IY))*TO5) * F
      ELSE
      GB(I)  = (GYR(NN,MM,IY) +GT (NN,MM)    *T  ) * F
      GB(I1) = (HYR(NN,MM,IY) +HT (NN,MM)    *T  ) * F
      ENDIF
      RNN = REAL(NN)
      GV(I)  = GB(I)  / RNN
      GV(I1) = GB(I1) / RNN
  200 I = I+2

  300 CONTINUE

      RETURN

!          Error trap diagnostics:
 9100 WRITE (0,'(''COFRM:  DATE'',F9.3,'' preceeds earliest available ('',F6.1,'')'')') DATE, EPOCH(1)
      CALL EXIT (1)
 9200 FORMAT('COFRM:  DATE',F9.3,' is after the last recommended for extrapolation (',F6.1,')')
      END

      SUBROUTINE DYPOL (COLAT,ELON,VP)
!          Computes parameters for dipole component of geomagnetic field.
!          COFRM must be called before calling DYPOL!
!          940504 A. D. Richmond
!
!          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed
!
!          RETURNS:
!            COLAT = Geocentric colatitude of geomagnetic dipole north pole
!                    (deg)
!            ELON  = East longitude of geomagnetic dipole north pole (deg)
!            VP    = Magnitude, in T.m, of dipole component of magnetic
!                    potential at geomagnetic pole and geocentric radius
!                    of 6371.0088 km
!
!          HISTORY:
!          Nov 2009: Change definition of earth's mean radius (RE) from 6371.2
!                    to the WGS84 value (6371.0088), by J.T. Emmert, NRL.

      PARAMETER (RTOD = 57.2957795130823, RE = 6371.0088)
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG

!          Compute geographic colatitude and longitude of the north pole of
!          earth centered dipole
      GPL   = SQRT (GB(2)**2 + GB(3)**2 + GB(4)**2)
      CTP   = GB(2) / GPL
      STP   = SQRT (1. - CTP*CTP)
      COLAT = ACOS (CTP) * RTOD
      ELON  = ATAN2 (GB(4),GB(3)) * RTOD

!          Compute magnitude of magnetic potential at pole, radius Re.
      VP = .2*GPL*RE
!          .2 = 2*(10**-4 T/gauss)*(1000 m/km) (2 comes through F0 in COFRM).

      RETURN
      END

      SUBROUTINE FELDG (IENTY,GLAT,GLON,ALT, BNRTH,BEAST,BDOWN,BABS)
!          Compute the DGRF/IGRF field components at the point GLAT,GLON,ALT.
!          COFRM must be called to establish coefficients for correct date
!          prior to calling FELDG.
!
!          IENTY is an input flag controlling the meaning and direction of the
!                remaining formal arguments:
!          IENTY = 1
!            INPUTS:
!              GLAT = Latitude of point (deg)
!              GLON = Longitude (east=+) of point (deg)
!              ALT  = Ht of point (km)
!            RETURNS:
!              BNRTH  north component of field vector (Gauss)
!              BEAST  east component of field vector  (Gauss)
!              BDOWN  downward component of field vector (Gauss)
!              BABS   magnitude of field vector (Gauss)
!
!          IENTY = 2
!            INPUTS:
!              GLAT = X coordinate (in units of earth radii 6371.0088 km )
!              GLON = Y coordinate (in units of earth radii 6371.0088 km )
!              ALT  = Z coordinate (in units of earth radii 6371.0088 km )
!            RETURNS:
!              BNRTH = X component of field vector (Gauss)
!              BEAST = Y component of field vector (Gauss)
!              BDOWN = Z component of field vector (Gauss)
!              BABS  = Magnitude of field vector (Gauss)
!          IENTY = 3
!            INPUTS:
!              GLAT = X coordinate (in units of earth radii 6371.0088 km )
!              GLON = Y coordinate (in units of earth radii 6371.0088 km )
!              ALT  = Z coordinate (in units of earth radii 6371.0088 km )
!            RETURNS:
!              BNRTH = Dummy variable
!              BEAST = Dummy variable
!              BDOWN = Dummy variable
!              BABS  = Magnetic potential (T.m)
!
!          INPUT from COFRM through COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
!            NMAX = Maximum order of spherical harmonic coefficients used
!            GB   = Coefficients for magnetic field calculation
!            GV   = Coefficients for magnetic potential calculation
!            ICHG = Flag indicating when GB,GV have been changed
!
!          HISTORY:
!          Apr 1983: written by Vincent B. Wickwar (Utah State Univ.).
!
!          May 1994 (A.D. Richmond): Added magnetic potential calculation
!
!          Oct 1995 (Barnes): Added ICHG
!
!          Nov 2009: Change definition of earth's mean radius (RE) from 6371.2
!                    to the WGS84 value (6371.0088), by J.T. Emmert, NRL.

      PARAMETER (DTOR = 0.01745329251994330, RE = 6371.0088)
      COMMON /MAGCOF/ NMAX,GB(255),GV(225),ICHG
      DIMENSION G(255), H(255), XI(3)
      SAVE IENTYP, G
      DATA IENTYP/-10000/

      IF (IENTY .EQ. 1) THEN
      IS   = 1
      RLAT = GLAT * DTOR
      CT   = SIN (RLAT)
      ST   = COS (RLAT)
      RLON = GLON * DTOR
      CP   = COS (RLON)
      SP   = SIN (RLON)
        CALL GD2CART (GLAT,GLON,ALT,XXX,YYY,ZZZ)
        XXX = XXX/RE
        YYY = YYY/RE
        ZZZ = ZZZ/RE
      ELSE
        IS   = 2
        XXX  = GLAT
        YYY  = GLON
        ZZZ  = ALT
      ENDIF
      RQ    = 1./(XXX**2+YYY**2+ZZZ**2)
      XI(1) = XXX*RQ
      XI(2) = YYY*RQ
      XI(3) = ZZZ*RQ
      IHMAX = NMAX*NMAX+1
      LAST  = IHMAX+NMAX+NMAX
      IMAX  = NMAX+NMAX-1

      IF (IENTY .NE. IENTYP .OR. ICHG .EQ. 1) THEN
        IENTYP = IENTY
      ICHG = 0
        IF (IENTY .NE. 3) THEN
        DO 10 I=1,LAST
   10     G(I) = GB(I)
        ELSE
        DO 20 I=1,LAST
   20     G(I) = GV(I)
        ENDIF
      ENDIF

      DO 30 I=IHMAX,LAST
   30 H(I) = G(I)

      MK = 3
      IF (IMAX .EQ. 1) MK=1

      DO 100 K=1,MK,2
      I  = IMAX
      IH = IHMAX

   60 IL = IH-I
      F = 2./FLOAT(I-K+2)
      X = XI(1)*F
      Y = XI(2)*F
      Z = XI(3)*(F+F)

      I = I-2
      IF (I .LT. 1) GO TO 90
      IF (I .EQ. 1) GO TO 80

      DO 70 M=3,I,2
      IHM = IH+M
      ILM = IL+M
      H(ILM+1) = G(ILM+1)+ Z*H(IHM+1) + X*(H(IHM+3)-H(IHM-1)) &
                                              -Y*(H(IHM+2)+H(IHM-2))
   70 H(ILM)   = G(ILM)  + Z*H(IHM)   + X*(H(IHM+2)-H(IHM-2)) &
                                              +Y*(H(IHM+3)+H(IHM-1))

   80 H(IL+2) = G(IL+2) + Z*H(IH+2) + X*H(IH+4) - Y*(H(IH+3)+H(IH))
      H(IL+1) = G(IL+1) + Z*H(IH+1) + Y*H(IH+4) + X*(H(IH+3)-H(IH))

   90 H(IL)   = G(IL)   + Z*H(IH)   + 2.*(X*H(IH+1)+Y*H(IH+2))
      IH = IL
      IF (I .GE. K) GO TO 60
  100 CONTINUE

      S = .5*H(1)+2.*(H(2)*XI(3)+H(3)*XI(1)+H(4)*XI(2))
      T = (RQ+RQ)*SQRT(RQ)
      BXXX = T*(H(3)-S*XXX)
      BYYY = T*(H(4)-S*YYY)
      BZZZ = T*(H(2)-S*ZZZ)
      BABS = SQRT(BXXX**2+BYYY**2+BZZZ**2)
      IF (IS .EQ. 1) THEN            ! (convert back to geodetic)
        BEAST = BYYY*CP-BXXX*SP
        BRHO  = BYYY*SP+BXXX*CP
        BNRTH =  BZZZ*ST-BRHO*CT
        BDOWN = -BZZZ*CT-BRHO*ST
      ELSEIF (IS .EQ. 2) THEN        ! (leave in earth centered cartesia
        BNRTH = BXXX
        BEAST = BYYY
        BDOWN = BZZZ
      ENDIF

!          Magnetic potential computation makes use of the fact that the
!          calculation of V is identical to that for r*Br, if coefficients
!          in the latter calculation have been divided by (n+1) (coefficients
!          GV).  Factor .1 converts km to m and gauss to tesla.
      IF (IENTY.EQ.3) BABS = (BXXX*XXX + BYYY*YYY + BZZZ*ZZZ)*RE*.1

      RETURN
      END

      SUBROUTINE SUBSOL (IYR,IDAY,IHR,IMN,SEC,SBSLLAT,SBSLLON)
!          Find subsolar geographic latitude and longitude given the
!          date and time (Universal Time).
!
!          This is based on formulas in Astronomical Almanac for the
!          year 1996, p.  C24. (U.S.  Government Printing Office,
!          1994).  According to the Almanac, results are good to at
!          least 0.01 degree latitude and 0.025 degree longitude
!          between years 1950 and 2050.  Accuracy for other years has
!          not been tested although the algorithm has been designed to
!          accept input dates from 1601 to 2100.  Every day is assumed
!          to have exactly 86400 seconds; thus leap seconds that
!          sometimes occur on June 30 and December 31 are ignored:
!          their effect is below the accuracy threshold of the algorithm.
!
!          961026 A. D. Richmond, NCAR
!
!          INPUTS:
!            IYR  = Year (e.g., 1994). IYR must be in the range: 1601 to 2100.
!            IDAY = Day number of year (e.g., IDAY = 32 for Feb 1)
!            IHR  = Hour of day    (e.g., 13 for 13:49)
!            IMN  = Minute of hour (e.g., 49 for 13:49)
!            SEC  = Second and fraction after the hour/minute.
!          Note:  While IYR is bounds tested, there is no constraint
!                 placed on values: IDAY,IHR,IMN,SEC; e.g., IHR=25 is
!                 properly interpreted.
!          RETURNS:
!            SBSLLAT = geographic latitude of subsolar point (degrees)
!            SBSLLON = geographic longitude of subsolar point (degrees,
!                      between -180 and +180)

      PARAMETER (D2R=0.0174532925199432957692369076847 , &
                 R2D=57.2957795130823208767981548147)
      PARAMETER (MSGUN=6)
      INTEGER IYR,YR,IDAY,IHR,IMN,NLEAP,NCENT,NROT
      REAL SEC,SBSLLAT,SBSLLON,L0,G0,DF,LF,GF,L,G,LAMBDA,EPSILON,N &
         ,GRAD,LAMRAD,SINLAM,EPSRAD,DELTA,UT,ETDEG
!
! Number of years from 2000 to IYR (negative if IYR < 2000):
      YR = IYR - 2000
!
! NLEAP (final) = number of leap days from (2000 January 1) to (IYR January 1)
!                 (negative if IYR is before 1997)
      NLEAP = (IYR-1601)/4
      NLEAP = NLEAP - 99
      IF (IYR.LE.1900) THEN
        IF (IYR.LE.1600) THEN
         WRITE(MSGUN,*) 'SUBSOLR INVALID BEFORE 1601: INPUT YEAR = ',IYR
         STOP
        ENDIF
        NCENT = (IYR-1601)/100
        NCENT = 3 - NCENT
        NLEAP = NLEAP + NCENT
      ENDIF
      IF (IYR.GE.2101) THEN
        WRITE(MSGUN,*) 'SUBSOLR INVALID AFTER 2100:  INPUT YEAR = ',IYR
        STOP
      ENDIF
!
! L0 = Mean longitude of Sun at 12 UT on January 1 of IYR:
!     L0 = 280.461 + .9856474*(365*(YR-NLEAP) + 366*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = 280.461 + .9856474*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = (280.461 - 360.) + (.9856474*365 - 360.)*(YR-4*NLEAP)
!          + (.9856474*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
      L0 = -79.549 + (-.238699*(YR-4*NLEAP) + 3.08514E-2*NLEAP)
!
! G0 = Mean anomaly at 12 UT on January 1 of IYR:
!     G0 = 357.528 + .9856003*(365*(YR-NLEAP) + 366*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = 357.528 + .9856003*(365*(YR-4*NLEAP) + (366+365*3)*NLEAP)
!          - (ARBITRARY INTEGER)*360.
!        = (357.528 - 360.) + (.9856003*365 - 360.)*(YR-4*NLEAP)
!          + (.9856003*(366+365*3) - 4*360.)*NLEAP,
!  where ARBITRARY INTEGER = YR+1.  This gives:
      G0 = -2.472 + (-.2558905*(YR-4*NLEAP) - 3.79617E-2*NLEAP)
!
! Universal time in seconds:
      UT = FLOAT(IHR*3600 + IMN*60) + SEC
!
! Days (including fraction) since 12 UT on January 1 of IYR:
      DF = (UT/86400. - 1.5) + IDAY
!
! Addition to Mean longitude of Sun since January 1 of IYR:
      LF = .9856474*DF
!
! Addition to Mean anomaly since January 1 of IYR:
      GF = .9856003*DF
!
! Mean longitude of Sun:
      L = L0 + LF
!
! Mean anomaly:
      G = G0 + GF
      GRAD = G*D2R
!
! Ecliptic longitude:
      LAMBDA = L + 1.915*SIN(GRAD) + .020*SIN(2.*GRAD)
      LAMRAD = LAMBDA*D2R
      SINLAM = SIN(LAMRAD)
!
! Days (including fraction) since 12 UT on January 1 of 2000:
      N = DF + FLOAT(365*YR + NLEAP)
!
! Obliquity of ecliptic:
      EPSILON = 23.439 - 4.E-7*N
      EPSRAD = EPSILON*D2R
!
! Right ascension:
      ALPHA = ATAN2(COS(EPSRAD)*SINLAM,COS(LAMRAD))*R2D
!
! Declination:
      DELTA = ASIN(SIN(EPSRAD)*SINLAM)*R2D
!
! Subsolar latitude:
      SBSLLAT = DELTA
!
! Equation of time (degrees):
      ETDEG = L - ALPHA
      NROT = NINT(ETDEG/360.)
      ETDEG = ETDEG - FLOAT(360*NROT)
!
! Apparent time (degrees):
      APTIME = UT/240. + ETDEG
!          Earth rotates one degree every 240 s.
!
! Subsolar longitude:
      SBSLLON = 180. - APTIME
      NROT = NINT(SBSLLON/360.)
      SBSLLON = SBSLLON - FLOAT(360*NROT)
!
      RETURN
      END

      SUBROUTINE MAGLOCTM (ALON,SBSLLAT,SBSLLON,CLATP,POLON,MLT)

!  Computes magnetic local time from magnetic longitude, subsolar coordinates,
!   and geomagnetic pole coordinates.
!  950302 A. D. Richmond, NCAR
!  Algorithm:  MLT is calculated from the difference of the apex longitude,
!   alon, and the geomagnetic dipole longitude of the subsolar point.
!
!   Inputs:
!    ALON    = apex magnetic longitude of the point (deg)
!    SBSLLAT = geographic latitude of subsolar point (degrees)
!    SBSLLON = geographic longitude of subsolar point (degrees)
!    CLATP   = Geocentric colatitude of geomagnetic dipole north pole (deg)
!    POLON   = East longitude of geomagnetic dipole north pole (deg)
!
!   Output:
!    mlt (real) = magnetic local time for the apex longitude alon (hours)
!
! To go from mlt to alon (see comments following Entry mlt2alon for definition
!  of variables), use:
!
!     CALL MLT2ALON (MLT,SBSLLAT,SBSLLON,CLATP,POLON,ALON)
!
!  NOTE: If the calling routine uses subroutine magloctm in conjunction with
!   file magfld.f (which is used by subroutine APEX), then clatp and polon can
!   be found by invoking
!
!     CALL DYPOL (CLATP,POLON,VP)
!
!   where vp is an unneeded variable.  (Note that subroutine COFRM must have
!   been called before DYPOL, in order to set up the coefficients for the
!   desired epoch.)  Alternatively, if subroutine apxntrp is
!   used to get alon from previously computed arrays, then
!   clatp and polon can be obtained for use in magloctm by adding
!
!     COMMON /APXDIPL/ CLATP,POLON,DUM1,DUM2,DUM3
!
!   to the calling routine (where DUM1,DUM2,DUM3 are unneeded dummy variables).
!
        REAL MLT
        CALL SOLGMLON (SBSLLAT,SBSLLON,CLATP,POLON,SMLON)
        MLT = (ALON - SMLON)/15.0 + 12.0
        IF (MLT .GE. 24.0) MLT = MLT - 24.0
        IF (MLT .LT.   0.) MLT = MLT + 24.0
        RETURN
!
      ENTRY MLT2ALON (XMLT,SBSLLAT,SBSLLON,CLATP,POLON,ALONX)
!
!   Inputs:
!    XMLT (real) = magnetic local time for the apex longitude alon (hours,
!                 0. to 24.)
!    SBSLLAT     = geographic latitude of subsolar point (degrees)
!    SBSLLON     = geographic longitude of subsolar point (degrees)
!    CLATP       = Geocentric colatitude of geomagnetic dipole north pole (deg)
!    POLON       = East longitude of geomagnetic dipole north pole (deg)
!
!   Output:
!    ALONX       = apex magnetic longitude of the point (deg, -180. to 180.)
!
        CALL SOLGMLON (SBSLLAT,SBSLLON,CLATP,POLON,SMLON)
        ALONX = 15.*(XMLT - 12.0) + SMLON
        IF (ALONX .GT.  180.) ALONX = ALONX - 360.0
        IF (ALONX .LE. -180.) ALONX = ALONX + 360.0
        RETURN
        END

      SUBROUTINE SOLGMLON (XLAT,XLON,COLAT,ELON,MLON)
! Computes geomagnetic longitude of the point with geocentric spherical
!  latitude and longitude of XLAT and XLON, respectively.
! 940719 A. D. Richmond, NCAR
! Inputs:
!   XLAT  = geocentric spherical latitude (deg)
!   XLON  = geocentric spherical longitude (deg)
!   COLAT = Geocentric colatitude of geomagnetic dipole north pole (deg)
!   ELON  = East longitude of geomagnetic dipole north pole (deg)
! Output:
!   MLON  = Geomagnetic dipole longitude of the point (deg, -180. to 180.)

      REAL MLON
      PARAMETER (RTOD=5.72957795130823E1,DTOR=1.745329251994330E-2)

! Algorithm:
!   Use spherical coordinates.
!   Let GP be geographic pole.
!   Let GM be geomagnetic pole (colatitude COLAT, east longitude ELON).
!   Let XLON be longitude of point P.
!   Let TE be colatitude of point P.
!   Let ANG be longitude angle from GM to P.
!   Let TP be colatitude of GM.
!   Let TF be arc length between GM and P.
!   Let PA = MLON be geomagnetic longitude, i.e., Pi minus angle measured
!     counterclockwise from arc GM-P to arc GM-GP.
!   Then, using notation C=cos, S=sin, spherical-trigonometry formulas
!     for the functions of the angles are as shown below.  Note: STFCPA,
!     STFSPA are sin(TF) times cos(PA), sin(PA), respectively.

      CTP = COS(COLAT*DTOR)
      STP = SQRT(1. - CTP*CTP)
      ANG = (XLON-ELON)*DTOR
      CANG = COS(ANG)
      SANG = SIN(ANG)
      CTE = SIN(XLAT*DTOR)
      STE = SQRT(1.-CTE*CTE)
      STFCPA = STE*CTP*CANG - CTE*STP
      STFSPA = SANG*STE
      MLON = ATAN2(STFSPA,STFCPA)*RTOD
      RETURN
      END

!.  SUBROUTINE TMJD ( YEAR, MONTH, DAY, HOUR, MINUTE, SEC, MJD, FLAG )
!.
!.  This subroutine converts between the following time formats:
!.
!.  (1)  YEAR, MONTH, DAY, HOUR, MINUTE (INTEGER*4) and SEC (REAL*8)
!.
!.  (2)  MJD (REAL*8)  Modified Julian Days, referred to 2000-01-01, 0 UT
!.
!.  The value of FLAG (INTEGER*4) controls the direction of conversion:
!.       FLAG.GT.0 :   input (1), output (2)
!.       FLAG.LE.0 :   input (2), output (1)
!.
!.  Subroutines from lib/igm/util:   TJUDA, IDBL, DBL
!
!   Joerg Warnecke  IGM TU BS  1995-03-22   (joerg@geophys.nat.tu-bs.de)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!
      SUBROUTINE TMJD ( YEAR, MONTH, DAY, HOUR, MINUTE, SEC, &
                        MJD, FLAG )
!
      IMPLICIT NONE
!
      DOUBLE PRECISION  JD0
      PARAMETER        (JD0=2451544.5D0)
!                       Julian Date of 2000-01-01, 0 UT
!
      INTEGER           YEAR, MONTH, DAY, HOUR, MINUTE, FLAG, &
                        I, J, K, L, M, N, IDBL
      DOUBLE PRECISION  SEC, MJD, DBL, JD
      REAL              UT
!
      IF (FLAG.GT.0) THEN
!
        IF ((YEAR.GE.1950).AND.(YEAR.LT.2100)) THEN
!         taken from ESA/ESTEC subroutine `JD2000.FOR'
          J = (14 - MONTH)/12
          L = YEAR - J - 1900*(YEAR/1900) + 100*(2000/(YEAR+1951))
          I = DAY - 36496 + (1461*L)/4 + (367*(MONTH-2+J*12))/12
        ELSE
          CALL TJUDA ( YEAR, MONTH, DAY, 0.0, JD, 1 )
          I = NINT(JD-JD0)
        END IF
        MJD = DBLE(I) + SEC/86400.D0 + DBLE(MINUTE)/1440.D0 &
                                     + DBLE(HOUR)/24.D0
!
      ELSE
!
        IF ((MJD.GT.-18262.D0).AND.(MJD.LT.36525.D0)) THEN
!         taken from ESA/ESTEC subroutine `DJ2000.FOR'
          I = INT(MJD+18262.D0)
          L = (4000*(I+18204))/1461001
          N = I - (1461*L)/4 + 18234
          M = (80*N)/2447
          DAY = N - (2447*M)/80
          J = M/11
          MONTH = M + 2 - 12*J
          YEAR = 1900 + L + J
          SEC = (MJD - DBLE(I-18262))*24.D0
        ELSE
          CALL TJUDA ( YEAR, MONTH, DAY, UT, MJD+JD0, 0 )
          CALL TJUDA ( YEAR, MONTH, DAY, 0.0, JD, 1 )
          I = NINT(JD-JD0)
          SEC = (MJD - DBLE(I))*24.D0
        END IF
        HOUR = INT(SEC)
        SEC = (SEC - DBLE(HOUR))*60.D0
        MINUTE = INT(SEC)
        SEC = (SEC - DBLE(MINUTE))*60.D0
!
!       number of significant decimals of a double precision number
        K = IDBL(0.11111111111111111111D0,+1)
!       number of significant decimals of MJD after the dot
        K = K - IDBL(MJD,-1)
!       number of significant decimals of SEC after the dot = K - 6
!       because a day has 86400 or about 10^6 seconds
        SEC = DBL(SEC,K-6)
!
      END IF
!
      RETURN
      END
!==============================================================================
!------------------------------------------------------------------------------
!*******************************************************************C
!
!.  INTEGER*4 FUNCTION IDBL ( D, I )
!.
!.  Bestimmt zu einer REAL*8-Zahl  D  je nach dem Wert des INTEGER*4-
!.  Schalters  I  folgende Parameter:
!.
!.  fuer I > 0:  die Anzahl von Nachkomma-Dezimalstellen
!.  fuer I < 0:  die Anzahl von Dezimalstellen vor dem Komma
!.  fuer I = 0:  die Gesamtzahl von Dezimalstellen, einschliesslich
!.               Dezimalpunkt und (bei negativen Zahlen) Vorzeichen.
!.               Bei Zahlen groesser als -1 und kleiner als +1  wird
!.               ausserdem eine fuehrende Null mitgezaehlt.
!.
!.  Es wird beruecksichtigt, dass unter SunOS nach dem IEEE-Standard
!.  754  eine REAL*8-Zahl etwa 16  signifikante Dezimalstellen hat.
!.
!   Beispiele             |  IFLT(D,+1)  |  IFLT(D,-1)  |  IFLT(D,0)
!   ----------------------+--------------+--------------+------------
!   +0.075                |     3        |      0       |      5
!   -0.075                |     3        |      0       |      6
!   10.249999999999999    |     2        |      2       |      5
!
!   Joerg Warnecke          26.11.1991
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!
      INTEGER*4 FUNCTION IDBL ( D, I )
!
      INTEGER*4          I
      REAL*8             D, REST, TEST
!
      IDBL = 0
!
      IF (I.LE.0) THEN
        REST = ABS(D)
        TEST = 0.999999999999999D0
        DO WHILE (REST.GT.TEST)
          IDBL = IDBL + 1
          REST = REST*0.1D0
        END DO
      END IF
!
      IF (I.EQ.0) THEN
        IF (IDBL.EQ.0) IDBL = 1
        IDBL = IDBL + 1
        IF (D.LT.0.D0) IDBL = IDBL + 1
      END IF
!
      IF (I.GE.0) THEN
        TEST = ABS(D) * 1.D-15
        REST = MOD(ABS(D),1.D0)
        DO WHILE (REST.GT.TEST)
          REST = REST*10.D0
          TEST = TEST*10.D0
          IDBL = IDBL + 1
          REST = ABS(REST-NINT(REST))
        END DO
      END IF
!
      RETURN
      END
!------------------------------------------------------------------------------
!
!.  SUBROUTINE TJUDA ( JAHR, MONAT, TAG, HOUR, JD, FLAG )
!.
!.  Umrechnung einer Zeitangabe in den Formen
!.  (1)  JAHR, MONAT, TAG (INTEGER*4) Datum christlicher Zeitrechnung
!.       HOUR                (REAL*4) Uhrzeit (UT) in Stunden
!.  und
!.  (2)  JD                  (REAL*8) Tage seit dem 1.1.4713 v.Chr.,
!.                                    12:00:00 UT
!.
!.  Man bezeichnet JD als "Julianisches Datum".
!.
!.  FLAG (INTEGER*4) gibt an, in welcher Richtung umgewandelt wird:
!.   FLAG.GT.0:   (1) -> (2)
!.   FLAG.LE.0:   (2) -> (1)
!
!   Joerg Warnecke        26.11.1991 n.Chr.
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!
      SUBROUTINE TJUDA ( JAHR, MONAT, TAG, HOUR, JD, FLAG )
!
!---- Parameter: erster Tag des Gregorianischen Kalenders ----------C
      INTEGER           JGREG, MGREG, TGREG
      PARAMETER         (JGREG=1582, MGREG=10, TGREG=15)
      DOUBLE PRECISION  DGREG
      PARAMETER         (DGREG=2299160.5D0)
!
      INTEGER           JAHR, MONAT, TAG, JA, MO, FLAG, J, N, M(12)
      REAL              HOUR
      DOUBLE PRECISION  JD, D
!
!---- Laenge der Monate des Jahres ---------------------------------C
      DATA  M / 31, 0, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 /
      M(2) = 28
!
      IF (FLAG.LE.0) THEN
!-------------------------------------------------------------------C
!---- Julianisches Datum -> Christliches Datum ---------------------C
!
!---- Umrechnung auf Julianische Jahre zu 365.25 Tagen -------------C
      IF (JD.GE.DGREG) THEN
!------ Anzahl Gregorianischer 400-Jahres-Zyklen -------------------C
        J = INT ( JD - 1721119.5D0 ) / 146097
        D = JD + DBLE (3*J)
!------ Uebrige volle Jahrhunderte (maximal 3, "kurz") -------------C
        J = MIN ( (INT (JD-1721119.5D0) - J*146097) / 36524, 3 )
        D = D + DBLE (J)
!------ Korrektur, so dass 05.10.1582 = 15.10.1582 -----------------C
        D = D - 1.5D0
      ELSE
        D = JD + 0.5D0
      END IF
!
!---- Anzahl Julianischer 4-Jahres-Zyklen seit Anfang 4713 v.Chr. --C
      J = INT ( D ) / (4*365 + 1)
      D = D - DBLE ( J * (4*365 + 1) )
      N = INT ( D )
!
!---- Uebrige volle Jahre (das erste davon ist ein Schaltjahr !) ---C
      IF (N.GE.(366+365+365)) THEN
        D = D - DBLE (366+365+365)
        JA = 4*J + 3 - 4712
      ELSE IF (N.GE.(366+365)) THEN
        D = D - DBLE (366+365)
        JA = 4*J + 2 - 4712
      ELSE IF (N.GE.366) THEN
        D = D - DBLE (366)
        JA = 4*J + 1 - 4712
      ELSE
        JA = 4*J - 4712
        IF ( (JD.GE.DGREG) .AND. &
             (MOD(JA,100).EQ.0).AND.(MOD(JA,400).NE.0) ) THEN
          IF (N.GT.59) D = D - 1.D0
        ELSE
!         JA ist ein Schaltjahr
          M(2) = 29
        END IF
      END IF
!
!---- Nummer des Monats --------------------------------------------C
      MO = 1
      DO WHILE ((INT(D).GE.M(MO)).AND.(MO.LE.12))
        D = D - DBLE ( M(MO) )
        MO = MO + 1
      END DO
!
!---- Christliches Datum (JA = 0 entspricht dem Jahr 1 v.Chr.) ----C
      HOUR = REAL ( D - INT(D) ) * 24.0
      TAG = INT(D) + 1
      MONAT = MO
!
      IF (JA.GT.0) THEN
        JAHR = JA
      ELSE
        JAHR = JA - 1
      END IF
!
      ELSE
!-------------------------------------------------------------------C
!---- Christliches Datum -> Julianisches Datum ---------------------C
!
!---- Jahreszahl ( JA = 0 entspricht dem Jahr 1 v.Chr. ) -----------C
      IF (JAHR.GT.0) THEN
        JA = JAHR
      ELSE
        JA = JAHR + 1
      END IF
!
!---- Nummer des Monats --------------------------------------------C
      MO = MONAT
      DO WHILE (MO.GT.12)
        JA = JA + 1
        MO = MO - 12
      END DO
      DO WHILE (MO.LT.1)
        JA = JA - 1
        MO = MO + 12
      END DO
!
!---- Tage bis zum Ende des vorigen Jahres (Julianisch) ------------C
      J = JA - 1
      N = (J + 4713) * 365  +  (J + 4716) / 4
!
!---- Ueberpruefen ob Schaltjahr (Julianisch) ----------------------C
      IF (MOD(JA,4).EQ.0) M(2) = 29
!
!---- Korrektur fuer Gregorianischen Kalender ----------------------C
      IF ( (JA.GT.JGREG) .OR. &
          ((JA.EQ.JGREG).AND.(MO.GT.MGREG)) .OR. &
          ((JA.EQ.JGREG).AND.(MO.EQ.MGREG).AND.(TAG.GE.TGREG)) ) &
        THEN
        N = N - J/100 + J/400 + 2
        IF ((MOD(JA,100).EQ.0).AND.(MOD(JA,400).NE.0)) M(2) = 28
      END IF
!
!---- Aufsummieren der Tage bis Anfang des Monats ------------------C
      IF (MO.GT.1) THEN
        DO J = 1, MO-1
          N = N + M(J)
        END DO
      END IF
!
!---- Julianisches Datum -------------------------------------------C
      JD = DBLE ( N + TAG )  +  DBLE ( HOUR ) / 24.D0  -  1.5D0
!
      END IF
      RETURN
      END
!------------------------------------------------------------------------------
!
!.  SUBROUTINE TMJD ( YEAR, MONTH, DAY, HOUR, MINUTE, SEC, MJD, FLAG )
!.
!.  This subroutine converts between the following time formats:
!.
!.  (1)  YEAR, MONTH, DAY, HOUR, MINUTE (INTEGER*4) and SEC (REAL*8)
!.
!.  (2)  MJD (REAL*8)  Modified Julian Days, referred to 2000-01-01, 0 UT
!.
!.  The value of FLAG (INTEGER*4) controls the direction of conversion:
!.       FLAG.GT.0 :   input (1), output (2)
!.       FLAG.LE.0 :   input (2), output (1)
!.
!.  Subroutines from lib/igm/util:   TJUDA, IDBL, DBL
!
!   Joerg Warnecke  IGM TU BS  1995-03-22   (joerg@geophys.nat.tu-bs.de)
!
!++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C
!
!*******************************************************************C
!
!.  REAL*8 FUNCTION DBL ( D, I )
!.
!.  Rundet eine REAL*8-Zahl D auf I (INT*4) Nachkommastellen (I > 0)
!.  oder bis vor die I-te Vorkommastelle (I < 0).
!.
!   Joerg Warnecke          26.11.1991
!
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++C

      REAL*8 FUNCTION DBL ( D, I )
!
      REAL*8          IMAX
      PARAMETER      (IMAX=2.147483647D9)
!
      INTEGER*4       I
      REAL*8          D, E
      CHARACTER*32    STRING
!
      IF (I.GT.0) THEN
        E = 10.D0**I
      ELSE IF (I.LT.0) THEN
        E = 0.1D0**ABS(I)
      ELSE
        E = 1.0D0
      END IF
!
      IF (ABS(D*E).LE.IMAX) THEN
        DBL = DBLE ( NINT ( D*E ) ) / E
      ELSE
        DBL = D
        WRITE (STRING,1000,ERR=999) D*E
        READ (STRING,1000,ERR=999) DBL
        DBL = DBL/E
      END IF
!
 999  RETURN
!
 1000 FORMAT ( F31.0 )
      END
