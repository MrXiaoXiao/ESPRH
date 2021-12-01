	subroutine partials(fn_srcpar,
     &	nsrc, src_cusp, src_lat, src_lon, src_dep,
     &	nsta, sta_lab, sta_lat, sta_lon,
     &	mod_nl, mod_ratio, mod_v, mod_top,
     &	tmp_ttp, tmp_tts,
     &	tmp_xp, tmp_yp, tmp_zp)

c Compute partial derivatives for straight ray paths (option for 1D layered 
c models commented out). Velocity is taken from first entry in the model 
c specification array in hypoDD.inp.

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
	real		vs(20)

	parameter(pi=3.141593)

	character rcsid*150
	data rcsid /"$Header: /home1/crhet/julian/HYPODD/src/hypoDD/RCS/partials.f,v 1.7 2001/02/17 23:54:57 julian Exp julian $"/
	save rcsid

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

c the following loop needs some more work. I replaced the ray tracing 
c with a straight ray path approach, which you can probably get away in 
c your case. However, you need to get the geometry between source and receiver
c straight; i.e. you need to include receiver depth (probably easiest to 
c shift the source by the amount of receiver depth) to get proper incidence
c angle (ain) and distance (dist). I think this is the only place you need to 
c change anything. of course, any ray trace can be plugged in here.  
      do i=1,nsta
         do j=1,nsrc
            call delaz2(src_lat(j), src_lon(j), sta_lat(i), sta_lon(i), 
     &                 del, dist, az)

cc           1D ray tracing
c            call ttime(dist, src_dep(j), mod_nl, mod_v, mod_top, 
c     &                 tmp_ttp(i, j), ain)
c            call ttime(dist, src_dep(j), mod_nl, vs, mod_top, 
c     &                 tmp_tts(i, j), ain)

c	    Straight ray path: 
c here you need to adopt for receiver depth, i.e. src_dep-sta_dep
            tmp_ttp(i,j)= sqrt(dist**2 + src_dep(j)**2)/mod_v(1)
            tmp_tts(i,j)= sqrt(dist**2 + src_dep(j)**2)/vs(1)
            ain= 180-(atan(dist/src_dep(j))*57.2958)


cc           Determine wave speed (k) at the hypocenter
c            do k=1,mod_nl
c               if (src_dep(j).le.mod_top(k)) goto 10	! break
c            enddo
c10          continue

c	    Straight ray path approach, only first entry in velocity 
c           array taken:
            k= 2 

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
