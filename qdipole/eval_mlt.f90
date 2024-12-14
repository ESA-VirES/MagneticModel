!***********************************************************************
!
!  File Name: eval_mlt.f90
!  Authors: Martin Paces <martin.paces@eox.at>
!
!  Date: 14-12-2024
!  Version: 2.0
!  Description:
!      Evaluate magnetic local time.
!      Wrapper around the low-level subroutines.
!      For more details see MAGLOCTM in apex.f90
!
! Latest changes:
! 09-06-2016: First version of the subroutine.
! 14-12-2024: Review + removal loading of QD coefficients.
!
!***********************************************************************
!
!  eval_mlt
!
!  Evaluate Magnetic Local Time for given QD-longitude and MJD2000 time
!
!  call eval_mlt(t_mlt, qdlon, t_mjd2k, n_data)
!
!  INPUT ARGUMENTS
!    qdlon     array of Quasi-Dipole longitudes
!    t_mjd2k   array of MJD2000 times
!    n_data    array size
!
!  OUTPUT ARGUMENTS
!    t_mlt     array of Magnetic Local Time values.
!
!***********************************************************************

      subroutine eval_mlt(t_mlt, qdlon, t_mjd2k, n_data)

      implicit none

      real*8  qdlon(*), t_mjd2k(*), t_mlt(*)
      real*8  t_0, sec
      real*4  t_mlt4, qdlon4
      real*4  dp_colat, dp_elon, dp_vp
      real*4  sbsllat, sbsllon, ssec
      real*4  epoch, epoch_limit, epoch_old
      integer  n_data, i
      integer  iyr, imon, iday, idoy, ihour, imin
      parameter(epoch_limit = 1./365.25)

      epoch_old = -9999.9

      do i = 1, n_data
!        The epoch is approximated by neglecting the leap years
         epoch = t_mjd2k(i) / 365.25 + 2000
!        Reload coefficients and evaluate the dipole pole location.
!        If time has changed by more than epoch_limit
         if (abs(epoch - epoch_old) .gt. epoch_limit) then
            call cofrm (epoch)
            call dypol(dp_colat, dp_elon, dp_vp)
            epoch_old = epoch
         endif

!        1. Ccalculate civil date
         call tmjd(iyr, imon, iday, ihour, imin, sec, t_mjd2k(i), -1)
!        2. Calculate mjd = nint(t_0) of first day of year 'iyr'
         call tmjd(iyr, 1, 1, 0, 0, 0.d0, t_0, 1)
!        3. Get day of year (doy)
         idoy = int(t_mjd2k(i)) - nint(t_0) + 1
!        4. Evaluate the sub-solar lat/lon
         ssec = sec
         call subsol(iyr, idoy, ihour, imin, ssec, sbsllat, sbsllon)
!        5. Evaluate magnetic local time
         qdlon4 = qdlon(i)
         call magloctm(qdlon4, sbsllat, sbsllon, dp_colat, dp_elon, t_mlt4)
         t_mlt(i) = t_mlt4
      enddo

      return
      end subroutine eval_mlt
