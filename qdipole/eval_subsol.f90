!***********************************************************************
!
!  File Name: eval_subsol.f90
!  Authors: Martin Paces <martin.paces@eox.at>
!
!  Date: 14-12-2024
!  Version: 1.0
!  Description:
!      Evaluate sub-solar geographic latitude and longitude.
!      Wrapper around the low-level subroutines.
!      For more details see SUBSOL in apex.f90
!  Note:
!      The SUBSOL subroutine expects the time in the Universal Time (UT)
!      while the MJD2000 is aligned to UTC. The offset caused by
!      the accumulated leap seconds is not compensated which leads to
!      an approx. 0.15 dg longitude offset for the currently accumulated
!      36 leap seconds.
!
!  References: n/a
!
! Latest changes:
! 13-05-2016: First version of the subroutine.
! 14-12-2024: Correcting documentation.
!
!***********************************************************************
!
!  eval_subsol
!
!  Evaluate sub-solar geographic latitude and longitude for the given
!  MJD2000 values.
!
!  call eval_subsol(sbslat, sbslon, t_mjd2k, n_data)
!
!  INPUT ARGUMENTS
!    t_mjd2k  vector of times (decimal year)
!    n_data   number of points, i.e., vector size
!
!  OUTPUT ARGUMENTS
!    sbsllat    array of sub-solar latitudes
!    sbsllon    array of sub-solar longitudes
!
!***********************************************************************

      subroutine eval_subsol(sbsllat, sbsllon, t_mjd2k, n_data)
      implicit none

      real*8        sbsllat(*), sbsllon(*), t_mjd2k(*)
      real*8        t_0, sec
      real*4        lat, lon, ssec
      integer       n_data, i
      integer       iyr, imon, iday, idoy, ihour, imin

      do i = 1, n_data
!        1. Calculate civil date
         call tmjd(iyr, imon, iday, ihour, imin, sec, t_mjd2k(i), -1)
!        2. Calculate mjd = nint(t_0) of first day of year 'iyr'
         call tmjd(iyr, 1, 1, 0, 0, 0.d0, t_0, 1)
!        3. Get day of year (doy)
         idoy = int(t_mjd2k(i)) - nint(t_0) + 1
!        4. Evaluate the sub-solar lat/lon
         ssec = sec
         call subsol(iyr, idoy, ihour, imin, ssec, lat, lon)
         sbsllat(i) = lat
         sbsllon(i) = lon
      enddo

      return
      end subroutine eval_subsol
