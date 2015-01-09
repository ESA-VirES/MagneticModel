!***************************************************************************************************
!
!  File Name: qdipole.f90
!  Authors: unknown (probably Nils Olsen)
!
!  Date: n/a
!  Version: n/a
!  Description: Calculation of the apex magentic field coordinates. 
!               Wrapper around the low-level subroutines.
!               Derived from the Matlab/MEX code.
!  References: n/a
! 
! Latest changes: 
! 18-01-2011: calculation of MLT corrected
! 08-01-2015: Removing the Matlab/MEX interface.
! 08-01-2015: Changing the input time format from MJD2000 to decimal year.
!
!***************************************************************************************************
!
!  make_appex
! Computes APEX coordinates 
!
!  Converts geodetic to Quasi-Dipole coordinates
!  Must call LOADAPXSH first. All arguments are single precision, except for VECFLAG, which is a
!  4-byte integers.
!
!  call make_apex(qdlat,qdlon,xmlt,f11,f12,f21,f22, &
!                t,gc_rad,gc_lat, gd_lon,N_data,filnam) 
!
!  INPUT ARGUMENTS
!    t        vector of times (decimal year)               
!    gc_rad   vector of geocentric radii (km)
!    gc_lat   vector of geocentric latitudes (degrees)
!    gd_lon   vector of geocentric longitudes (degrees)
!    N_data   number of points, i.e., vector size
!    filnam   path to the coeficients (string up to 128 characters) 
!
!  OUTPUT ARGUMENTS
!    qdlat     array of Quasi-Dipole latitudes
!    qdlon     array of Quasi-Dipole longitudes
!    xmlt      array of magnetic local times
!    f11,f12   arrays of components of the Quasi-Dipole coordinate
!               base vectors F1 (eastward).
!    f21,f22   arrays of components of the Quasi-Dipole coordinate
!               base vectors F2 (northward).
!
!***************************************************************************************************

      subroutine make_apex(qdlat,qdlon,xmlt,f11,f12,f21,f22, &
      t,gc_rad,gc_lat, gd_lon,N_data,filnam) 

      real*8        qdlon(*), qdlat(*), xmlt(*), f11(*), f12(*), f21(*), f22(*)
      real*8        t(*), gc_rad(*), gc_lat(*), gd_lon(*)
      real*8        t_0, sec, sec0
      real*4        xf1(2), xf2(2)
      integer       N_data, vecflag
      character*128 filnam

      parameter(epoch_limit = 1./365.25, vecflag = 1)

!! DEBUG PRINT - START 
!      print *, "filename: ", filnam
!      print *, "number of points: ", N_data
!! DEBUG PRINT - STOP 

      epoch_old = -9999.9
      do i = 1, N_data
!! DEBUG PRINT - START 
!         print *, 1, t(i), gc_lat(i), gd_lon(i), gc_rad(i) 
!! DEBUG PRINT - STOP 
         call CONVRT (4,GLAT,ALT,sngl(gc_lat(i)), sngl(gc_rad(i)))
         glon = gd_lon(i)
!         epoch = t(i)/365.25+2000
         epoch = t(i)
         epoch_dif = abs(epoch - epoch_old)
         
! read/update coefficient file if time has changed by more than epoch_limit
         if (epoch_dif .gt. epoch_limit) then       
            call loadapxsh(filnam, epoch)
            CALL COFRM (epoch)
            CALL DYPOL (DP_COLAT,DP_ELON,DP_VP)      
            epoch_old = epoch
         endif
         
!
!        Calculate QD-latitude and -longitude
!
         call apxg2q(glat,glon,alt,vecflag,qlat,qlon,xf1,xf2,f)
         qdlat(i) = qlat
         qdlon(i) = qlon
         f11(i) = xf1(1)
         f12(i) = xf1(2)
         f21(i) = xf2(1)
         f22(i) = xf2(2)
!
!        Calculate MLT
!
!        IYR  = Year (e.g., 1994). IYR must be in the range: 1601 to 2100.
!        IDOY = Day number of year (e.g., IDAY = 32 for Feb 1)
!        IHR  = Hour of day    (e.g., 13 for 13:49)
!        IMN  = Minute of hour (e.g., 49 for 13:49)
!        SEC  = Second and fraction after the hour/minute.
         
!        here we need a reliable conversion mjd.frac -> civil date
!        1. calculate civil date
         call tmjd( iyr, imon,  iday,  ihour,  imin,  sec,  t(i),  -1)

!        2. calculate mjd = nint(t_0) of first day of year 'iyr'
         call tmjd( iyr, 1, 1, 0, 0, 0.d0, t_0,  1)

!        3. get day of year (doy)
         idoy  = int(t(i)) - nint(t_0) + 1
         rsec  = sec

         CALL SUBSOL(iyr, idoy, ihour, imin, rsec, SBSLLAT, SBSLLON)
         CALL MAGLOCTM(qlon, SBSLLAT, SBSLLON, DP_COLAT, DP_ELON, xmlt4) 
         xmlt(i) = xmlt4
      enddo
      
      return
     end subroutine make_apex
