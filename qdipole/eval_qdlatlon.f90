!***********************************************************************
!
!  File Name: eval_qdlatlon.f90
!  Authors: Martin Paces <martin.paces@eox.at>
!
!  Date: 14-12-2024
!  Version: 1.0
!  Description:
!      Evaluate apex magnetic field coordinates.
!      Wrapper around the low-level subroutines.
!      For more details see APXG2Q in apexsh.f90 
!
! Latest changes:
! 13-05-2016: First version of the subroutine.
! 14-12-2024: Correcting documentation.
!
!***********************************************************************
!
!  eval_qdlatlon
!
!  Evaluate quasi-dipole latitude and longitude for given points
!  in the geocentric polar coordinates and dates (decimal year).
!
!  call eval_qdlatlon(sbslat, sbslon, t_mjd2k, n_data)
!
!  INPUT ARGUMENTS
!    t_dy     vector of times (decimal year)
!    gcrad    vector of geocentric radii (km)
!    gclat    vector of geocentric latitudes (degrees)
!    gdlon    vector of geocentric longitudes (degrees)
!    n_data   number of points, i.e., vector size
!    coeff_file  path to the coefficients (string up to 1023 bytes)
!
!  OUTPUT ARGUMENTS
!    qdlat     array of Quasi-Dipole latitudes
!    qdlon     array of Quasi-Dipole longitudes
!
!***********************************************************************
!
!  eval_qdlatlonvb
!
!  Evaluate quasi-dipole latitude, longitude and base vectors for given
!  points in the geocentric polar coordinates and dates (decimal year).
!
!  call eval_qdlatlon(sbslat, sbslon, t_mjd2k, n_data)
!
!  INPUT ARGUMENTS
!    t_dy     vector of times (decimal year)
!    gcrad    vector of geocentric radii (km)
!    gclat    vector of geocentric latitudes (degrees)
!    gdlon    vector of geocentric longitudes (degrees)
!    n_data   number of points, i.e., vector size
!    coeff_file  path to the coefficients (string up to 1023 bytes)
!
!  OUTPUT ARGUMENTS
!    qdlat     array of Quasi-Dipole latitudes
!    qdlon     array of Quasi-Dipole longitudes
!    f11,f12   arrays of components of the Quasi-Dipole coordinate
!               base vectors F1 (eastward).
!    f21,f22   arrays of components of the Quasi-Dipole coordinate
!               base vectors F2 (northward).
!    f         array of F = |F1 x F2| (see apexsh.f90)
!
!***********************************************************************

      subroutine eval_qdlatlon(qdlat, qdlon, t_dy, gcrad, gclat, gclon,&
                               n_data, coeff_file)
      implicit none

      real*8  qdlon(*), qdlat(*)
      real*8  t_dy(*), gcrad(*), gclat(*), gclon(*)
      real*4  gdalt, gdlat, gdlon, xqdlat, xqdlon, xf1(2), xf2(2), xf
      real*4  epoch, epoch_limit, epoch_old
      integer*4  vecflag, n_data, i
      character*1023 coeff_file
      parameter(epoch_limit = 1./365.25, vecflag = 0)

      epoch_old = -9999.9

      do i = 1, n_data
!        Reload coefficients if time has changed by more than epoch_limit
         epoch = t_dy(i)
         if (abs(epoch - epoch_old) .gt. epoch_limit) then
            call loadapxsh(coeff_file, epoch)
            epoch_old = epoch
         endif
!        Convert geocentric coordinates to geodetic coordinates
         call convrt(4, gdlat, gdalt, sngl(gclat(i)), sngl(gcrad(i)))
         gdlon = gclon(i)
!        Calculate QD-latitude and QD-longitude
         call apxg2q(gdlat, gdlon, gdalt, vecflag, xqdlat, xqdlon, &
                     xf1, xf2, xf)
         qdlat(i) = xqdlat
         qdlon(i) = xqdlon
      enddo

      return
      end subroutine eval_qdlatlon

!***********************************************************************

      subroutine eval_qdlatlonvb(qdlat, qdlon, f11, f12, f21, f22, f, &
        t_dy, gcrad, gclat, gclon, n_data, coeff_file)
      implicit none

      real*8  qdlon(*), qdlat(*), f11(*), f12(*), f21(*), f22(*), f(*)
      real*8  t_dy(*), gcrad(*), gclat(*), gclon(*)
      real*4  gdalt, gdlat, gdlon, xqdlat, xqdlon, xf1(2), xf2(2), xf
      real*4  epoch, epoch_limit, epoch_old
      integer*4  vecflag, n_data, i
      character*256 coeff_file
      parameter(epoch_limit = 1./365.25, vecflag = 1)

      epoch_old = -9999.9

      do i = 1, n_data
!        Reload coefficients if time has changed by more than epoch_limit
         epoch = t_dy(i)
         if (abs(epoch - epoch_old) .gt. epoch_limit) then
            call loadapxsh(coeff_file, epoch)
            epoch_old = epoch
         endif
!        Convert geocentric coordinates to geodetic coordinates
         call convrt(4, gdlat, gdalt, sngl(gclat(i)), sngl(gcrad(i)))
         gdlon = gclon(i)
!        Calculate QD-latitude and QD-longitude
         call apxg2q(gdlat, gdlon, gdalt, vecflag, xqdlat, xqdlon, &
                     xf1, xf2, xf)
         qdlat(i) = xqdlat
         qdlon(i) = xqdlon
         f11(i) = xf1(1)
         f12(i) = xf1(2)
         f21(i) = xf2(1)
         f22(i) = xf2(2)
         f(i) = xf
      enddo

      return
      end subroutine eval_qdlatlonvb
