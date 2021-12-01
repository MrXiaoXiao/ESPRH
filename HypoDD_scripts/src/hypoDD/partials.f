	subroutine partials(fn_srcpar,
     &	nsrc, src_cusp, src_lat, src_lon, src_dep,
     &	nsta, sta_lab, sta_lat, sta_lon,
     &	mod_nl, mod_ratio, mod_v, mod_top,
     &	tmp_ttp, tmp_tts,
     &	tmp_xp, tmp_yp, tmp_zp)

	implicit none

	include'hypoDD.inc'

c	Parameters:
	character	fn_srcpar*80	! Source-parameter file
	integer		nsrc		! No. of sources
	integer		src_cusp(MAXEVE)! [1..nsrc]
	doubleprecision	src_lat(MAXEVE)	! [1..nsrc]
	doubleprecision	src_lon(MAXEVE)	! [1..nsrc]
	real		src_dep(MAXEVE)	! [1..nsrc]
	integer		nsta		! No. of stations
	character	sta_lab(MAXSTA)*7! [1..nsta]
	real		sta_lat(MAXSTA)	! [1..nsta]
	real		sta_lon(MAXSTA)	! [1..nsta]
	integer		mod_nl		! No. of layers
	real		mod_ratio	! Vp/Vs
	real		mod_v(MAXLAY)	! [1..mod_nl]
	real		mod_top(MAXLAY)	! [1..mod_nl]
	real		tmp_ttp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_tts(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_xp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_yp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]
	real		tmp_zp(MAXSTA,MAXEVE)! [1..nsta,1..nsrc]

c	Local variables:
	real		ain
	real		az
	real		del
	real		dist
	integer		i, j, k
	integer		iunit		! Output unit number
	real		pi
	integer		trimlen
        real            vs(MAXLAY)

	parameter(pi=3.141593)

      iunit = 0
      if (trimlen(fn_srcpar).gt.1) then
c        Open source-parameter file
         call freeunit(iunit)
         open(iunit,file=fn_srcpar,status='unknown')
      endif

c     Make sure hypocenters don't fall on layer boundaries
      do i=1,nsrc
         do j=1,mod_nl
            if (abs(src_dep(i)-mod_top(j)).lt.0.0001)
     &         src_dep(i) = src_dep(i)-0.001
         enddo
      enddo

c     Get S velocity model
      do i=1,mod_nl
         vs(i) = mod_v(i)/mod_ratio
      enddo

c     Compute epicentral distances, azimuths, angles of incidence,
c     and P/S-travel times from sources to stations
      do i=1,nsta
         do j=1,nsrc
            call delaz2(src_lat(j), src_lon(j), sta_lat(i), sta_lon(i), 
     &                 del, dist, az)

c           1D ray tracing
            call ttime(dist, src_dep(j), mod_nl, mod_v, mod_top, 
     &                 tmp_ttp(i, j), ain)
            call ttime(dist, src_dep(j), mod_nl, vs, mod_top, 
     &                 tmp_tts(i, j), ain)
            
c           Determine wave speed at the hypocenter
            do k=1,mod_nl
               if (src_dep(j).le.mod_top(k)) goto 10	! break
            enddo
10          continue

c           Depth derivative
            tmp_zp(i,j) = cos((ain * pi)/180.0)/mod_v(k-1)
c           Epicentral derivatives
	    tmp_xp(i,j) = (sin((ain * pi)/180.0) *
     &               cos(((az - 90) * pi)/180.0))/mod_v(k-1)
	    tmp_yp(i,j) = (sin((ain * pi)/180.0) *
     &               cos((az * pi)/180.0))/mod_v(k-1)

c           Write to source-parameter file
            if (iunit .ne. 0)
     &         write(iunit,'(i9,2x,f9.4,2x,f9.4,2x,a7,2x,f9.4,
     &         2x,f9.4,2x,f9.4)')
     &         src_cusp(j), src_lat(j), src_lon(j), sta_lab(i), 
     &         dist, az, ain

         enddo
      enddo

      if (iunit .ne. 0) close(iunit)	! Source-parameter file

      end !of subroutine partials
