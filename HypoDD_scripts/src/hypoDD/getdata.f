c Read in data

	subroutine getdata(
     &	log, fn_cc, fn_ct, fn_sta, fn_eve, fn_srcpar,
     &	idata, iphase, ncusp, icusp,
     &  maxdist,maxsep_ct,maxsep_cc,
     &	noisef_dt, mod_nl, mod_ratio, mod_v, mod_top,
     &	ev_date, ev_time, ev_cusp, ev_lat, ev_lon, ev_dep,
     &	ev_mag, ev_herr, ev_zerr, ev_res,
     &	sta_lab, sta_lat, sta_lon,
     &	dt_sta, dt_dt, dt_qual, dt_c1, dt_c2, dt_idx,
     &	dt_ista, dt_ic1, dt_ic2,dt_offs,
     &	nev, nsta, ndt, nccp, nccs, nctp, ncts,
     &	tmp_xp, tmp_yp, tmp_zp, tmp_ttp, tmp_tts)

	implicit none

	include 'hypoDD.inc'

c	Parameters:
	doubleprecision	atoangle	! ASCII-to-angle function
	integer		log
	character*80	fn_cc, fn_ct, fn_sta, fn_eve, fn_srcpar
	integer		idata
	integer		iphase
	integer		ncusp		! No. of events to relocate
	integer		icusp(MAXEVE)	! [1..ncusp] Keys of events to relocate
	real		maxdist
	real		maxsep_ct
	real		maxsep_cc
	real		noisef_dt
	integer		mod_nl
	real		mod_ratio
	real		mod_v(MAXLAY)	! [1..MAXLAY]
	real		mod_top(MAXLAY)	! [1..MAXLAY]
	integer		ev_date(MAXEVE)	! [1..MAXEVE]
	integer		ev_time(MAXEVE)	! [1..MAXEVE]
	integer		ev_cusp(MAXEVE)	! [1..nev] Event keys
	real		ev_lat(MAXEVE)	! [1..nev]
	real		ev_lon(MAXEVE)	! [1..nev]
	real		ev_dep(MAXEVE)	! [1..nev]
	real		ev_mag(MAXEVE)	! [1..nev]
	real		ev_herr(MAXEVE)	! [1..nev]
	real		ev_zerr(MAXEVE)	! [1..nev]
	real		ev_res(MAXEVE)	! [1..nev]
	character	sta_lab(MAXSTA)*7	! [1..MAXSTA]
	real		sta_lat(MAXSTA)	! [1..MAXSTA]
	real		sta_lon(MAXSTA)	! [1..MAXSTA]
	character	dt_sta(MAXDATA)*7	! [1..MAXDATA]
	real		dt_dt(MAXDATA)	! [1..MAXDATA]
	real		dt_qual(MAXDATA)	! [1..MAXDATA]
        real            dt_offs(MAXDATA)   	! [1..MAXDATA] 
	integer		dt_c1(MAXDATA)	! [1..MAXDATA]
	integer		dt_c2(MAXDATA)	! [1..MAXDATA]
	integer		dt_idx(MAXDATA)	! [1..MAXDATA]
	integer		dt_ista(MAXDATA)	! [1..MAXDATA]
	integer		dt_ic1(MAXDATA)	! [1..MAXDATA]
	integer		dt_ic2(MAXDATA)	! [1..MAXDATA]
	integer		nev
	integer		nsta
	integer		ndt
	integer		nccp
	integer		nccs
	integer		nctp
	integer		ncts
	integer		sscanf3		! Formatted string-reading function
	real		tmp_xp(MAXSTA,MAXEVE)! [1..MAXSTA, 1..MAXEVE]
	real		tmp_yp(MAXSTA,MAXEVE)! [1..MAXSTA, 1..MAXEVE]
	real		tmp_zp(MAXSTA,MAXEVE)! [1..MAXSTA, 1..MAXEVE]
	real		tmp_ttp(MAXSTA,MAXEVE)! [1..MAXSTA, 1..MAXEVE]
	real		tmp_tts(MAXSTA,MAXEVE)! [1..MAXSTA, 1..MAXEVE]

c	Local variables:

	real		azim
	character	buf1*20		! Input buffer
	character	buf2*20		! Input buffer
	real		clat
	real		clon
	integer		cusperr(34000)	! [1..nerr] Event keys to not locate
	character	dattim*25
	doubleprecision	elon(20), elat(20)
	integer		nerr
	real		del
	real		dist
	real		dt1, dt2
	real		dtn
	logical		ex
	integer		i
	integer		ic1
	integer		ic2
	integer		ifindi
	integer		ii
	integer		iicusp(MAXEVE)
	integer		iiotc
	integer		iskip
	integer		iunit
	integer		j
	integer		k
	integer		l
	character	line*200
	real		otc
	character	pha*1
	integer		sta_itmp(MAXSTA)
	character	str1*1
	real		tmp
	integer		trimlen
        real		offs
        real		dlat	
        real	 	dlon	
	integer		k1
	integer		k2
        real 		PI
        parameter       (PI=3.141593)

      call datetime (dattim)
      write (*,'("Reading data ...   ",a)') dattim
      write (log,'(/,"~ Reading data ...   ",a)') dattim

c--Read file with events not to be considered in the relocations
c--At the moment, hypoDD appears to simply warn you that you are trying
c  to relocate events you had marked as bad in a file called 'cuspid.err'
      inquire (FILE='cuspid.err',exist=ex)
      if (.not. ex) then
         nerr = 0
      else
         call freeunit (iunit)
         open (iunit,file='cuspid.err', status='unknown')
         i = 1
5        read (iunit,*,end=6) cusperr(i)
         i = i+1
         goto 5
6        continue
         nerr = i-1
         close(iunit)
      endif

c     Read event file
      call freeunit(iunit)
      open (iunit,file=fn_eve,status='unknown')
      i = 1

c--Begin earthquake read loop
c----------------------------------------------------------
10    read (iunit,'(a)',end=20) line
      read (line,*,err=1010) ev_date(i), ev_time(i), ev_lat(i),
     &   ev_lon(i), ev_dep(i), ev_mag(i), ev_herr(i), ev_zerr(i),
     &   ev_res(i), ev_cusp(i)

      if (ev_date(i).lt.10000000) ev_date(i) = ev_date(i)+19000000

c--If earthquake shallower than 10m, force it to 10m
c--This may not be necessary
      if (ev_dep(i) .lt. 0.01) ev_dep(i) = 1   ! no 0-depth for ttime!!!

c--Check if event is on error list. This appears to be just a warning.
      do j=1,nerr
         if (ev_cusp(i).eq.cusperr(j)) then
            write(*,*) 'NOTE: event in error list:', ev_cusp(i)
            write(log,*) 'NOTE: event in error list:', ev_cusp(i)
            goto 15
         endif
      enddo
15    continue

      if (ncusp.gt.1) then
c        Read selected events, skip others
         do j=1,ncusp
            if (icusp(j).eq.ev_cusp(i)) then
                i = i+1
                goto 16
            endif
         enddo
c        From now on, icusp free for work space
      else
c        Read all events
         i = i+1
      endif

16    continue
      if (i.gt.MAXEVE) stop'>>> Increase MAXEVE in hypoDD.inc.'
      goto 10
c----------------------------------------------------------
c--End earthquake read loop

20    nev = i-1
      write (*,'("# events = ",i5)') nev
      if (ncusp.gt.0 .and. ncusp.ne.nev) then
         write (*,'(//,">>> Events repeated in selection list '//
     &      'or missing/repeated in event file!")')
         do i=1,ncusp
            k = 0
            do j=1,nev
               if (icusp(i).eq.ev_cusp(j)) k = k+1
            enddo
            if (k.eq.0) write(*,*) icusp (i),' is missing.'
            if (k.ge.2) then
                write(*,*) icusp (i),' is non-unique.'  
                stop'Event ID must be unique!'
            endif
         enddo
      endif
      close(iunit)

c--Get center of event cluster
      clat = 0
      clon = 0
      do i=1,nev
         clat = clat + ev_lat(i)
         clon = clon + ev_lon(i)
      enddo
      clat = clat/nev
      clon = clon/nev

c--Read station list
      call freeunit (iunit)
      open (iunit,file=fn_sta,status='unknown')
      i = 1
      ii = 1

30    read (iunit,'(a)',end=40) line

c--Split into fields separated by white space
         if (sscanf3(line, "%s%s%s", sta_lab(i), buf1, buf2).ne.3) then
            write (6,*) line
            stop '** Bad station line'
         endif
         call rpad(sta_lab(i))

c--Convert strings to numbers, interpreting colons, if any.
         sta_lat(i) = atoangle(buf1)
         sta_lon(i) = atoangle(buf2)

c--Skip stations at distances larger than maxdist:
         call delaz(clat, clon, sta_lat(i), sta_lon(i), del, dist, azim)
         if (dist.le.maxdist) i = i+1
         if (i.gt.MAXSTA) then
            write (*,*)'>>> Increase station array dimension (MAXSTA)'//
     &      'in hypoDD.inc or decrease search radius for stations '//
     &      '(maxdist) in hypoDD.inp.'
            stop
         endif
         ii = ii+1
      goto 30

c--We now have read the entire station file
40    nsta = i-1
      write (log,'("# stations total = ",i6,/,
     & "# stations < maxdist = ",i6)')ii-1,nsta
      write (*,'("# stations < maxdist = ",i6)') nsta
      close(iunit)

c--Check for duplicated stations
      do i=1,nsta-1
         do j=i+1,nsta
            if (sta_lab(i).eq.sta_lab(j)) then
               write (*,*)sta_lab(i)
               stop'>>> This station is listed twice in station file.'
            endif
         enddo
      enddo

      if (idata.eq.0) goto 150   	!synthetics

      nccp = 0
      nccs = 0
      nctp = 0
      ncts = 0

c--Read cross-correlation dtimes
      call indexxi(nev,ev_cusp,iicusp)
      do i=1,nev
         icusp(i) = ev_cusp(iicusp(i)) ! icusp is workspace array here
      enddo
      i = 1
      iiotc = 0
      if ((idata.eq.1.or.idata.eq.3).and.trimlen(fn_cc).gt.1) then
         call freeunit(iunit)
         open (iunit,file=fn_cc,status='unknown')
50       read (iunit,'(a)',end=60) line
         if (line(1:1).eq.'#') then
            read (line,*,err=1051) str1, ic1, ic2, otc
            iskip = 0
c	skip event pairs with no origin time correction:
            if (abs(otc + 999).lt.0.001) then
               write (log,*)'No OTC for ', ic1, ic2, '. Pair skiped'
               iiotc = iiotc+1
               iskip = 1
               goto 50
            endif
c	skip event pairs with events not in event list: 
            k1= ifindi(nev,icusp,ic1)
            k2= ifindi(nev,icusp,ic2)
            if(k1.eq.0.or.k2.eq.0) then
               iskip=1 
               goto 50
            endif
c	skip event pairs with large separation distance: 
            dlat= ev_lat(iicusp(k1)) - ev_lat(iicusp(k2))
            dlon= ev_lon(iicusp(k1)) - ev_lon(iicusp(k2))
            offs= sqrt( (dlat*111)**2 +
     &            (dlon*(cos(ev_lat(iicusp(k1))*PI/180)*111))**2 +
     &            (ev_dep(iicusp(k1))-ev_dep(iicusp(k2)))**2)
            if(maxsep_cc.gt.0 .and. offs.gt.maxsep_cc) iskip= 1
            goto 50
         else
c           New format, body...
            if (iskip.eq.1) goto 50
            read (line,*,err=1051) dt_sta(i), dt_dt(i), dt_qual(i), pha
            dt_dt(i) = dt_dt(i) - otc
            dt_c1(i) = ic1
            dt_c2(i) = ic2
         endif

c--Skip far-away stations
         do j=1,nsta
            if (dt_sta(i).eq.sta_lab(j)) goto 58
         enddo
         goto 50

c--Only accept P or S phase codes
58       if (pha.eq.'P') then
            if (iphase.eq.2) goto 50
            dt_idx(i) = 1
            nccp = nccp+1
         elseif (pha.eq.'S') then
            if (iphase.eq.1) goto 50
            dt_idx(i) = 2
            nccs = nccs+1
         else
            stop '>>> Phase identifier format error.'
         endif
         dt_offs(i)= offs

         i = i+1
         if (i.gt.MAXDATA) stop'>>> Increase MAXDATA in hypoDD.inc.'
         goto 50

60       if (iphase.ne.2) then
            write (*,'("# cross corr P dtimes = ",i7,
     &      " (no OTC for",i7," event pairs)")') nccp, iiotc
            write (log,'("# cross corr P dtimes = ",i7,
     &      " (no org. time corr. for",i7," event pairs)")') nccp, iiotc
         endif
         if (iphase.ne.1) then
            write (*,'("# cross corr S dtimes = ",i7,
     &      " (no OTC for",i7," event pairs)")') nccs, iiotc
            write (log,'("# cross corr S dtimes = ",i7,
     &      " (no org. time corr. for",i7," event pairs)")') nccs, iiotc
         endif
         close(iunit)
      endif

      if (i.gt.MAXDATA) stop'>>> Increase MAXDATA in hypoDD.inc.'

c--Read catalog P and S absolute times
      if ((idata.eq.2.or.idata.eq.3) .and.
     &   trimlen(fn_ct).gt.1) then
         call freeunit(iunit)
         open (iunit,file=fn_ct,status='unknown')

90       read (iunit,'(a)',end=100) line
         if (line(1:1).eq.'#') then 	
            read(line,*,err=1091) str1, ic1, ic2

            iskip= 0
c	skip event pairs with events not in event list: 
            k1= ifindi(nev,icusp,ic1)
            k2= ifindi(nev,icusp,ic2)
            if(k1.eq.0.or.k2.eq.0) then
               iskip=1 
               goto 90
            endif
c	skip event pairs with large separation distance: 
            dlat= ev_lat(iicusp(k1)) - ev_lat(iicusp(k2))
            dlon= ev_lon(iicusp(k1)) - ev_lon(iicusp(k2))
            offs= sqrt( (dlat*111)**2 +
     &            (dlon*(cos(ev_lat(iicusp(k1))*PI/180)*111))**2 +
     &            (ev_dep(iicusp(k1))-ev_dep(iicusp(k2)))**2)
            if(maxsep_ct.gt.0 .and. offs.gt.maxsep_ct) iskip= 1
            goto 90
         else
            if (iskip.eq.1) goto 90
            read (line,*,err=1091)dt_sta(i), dt1,dt2, dt_qual(i),pha
            dt_c1(i) = ic1
            dt_c2(i) = ic2
         endif

c--Skip far-away data
         do j=1,nsta
            if (dt_sta(i).eq.sta_lab(j)) goto 95
         enddo
         goto 90

c--Store time difference
95       dt_dt(i) = dt1 - dt2
         if (pha.eq.'P') then
            if (iphase.eq.2) goto 90
            dt_idx(i) = 3
            nctp = nctp+1
         elseif (pha.eq.'S') then
            if (iphase.eq.1) goto 90
            dt_idx(i) = 4
            ncts = ncts+1
         else
            write(*,*)line
            stop '>>> Phase identifier format error.'
         endif
         dt_offs(i)= offs

         i = i+1
         if (i.gt.MAXDATA) stop'>>> Increase MAXDATA in hypoDD.inc.'
         goto 90

100      if (iphase.ne.2) then
             write (*,'("# catalog P dtimes = ",i7)') nctp
             write (log,'("# catalog P dtimes = ",i7)') nctp
         endif
         if (iphase.ne.1) then
             write (*,'("# catalog S dtimes = ",i7)') ncts
             write (log,'("# catalog S dtimes = ",i7)') ncts
         endif
         close(iunit)
      endif

      goto 160   ! jump over synthetics

150   continue
c--Generate synthetics dtimes
      if (nev.gt.20) stop'>>> Increase elon/elat array in getdata!'
c     Copy into double precision
      do i=1,nev
         elon(i) = ev_lon(i)
         elat(i) = ev_lat(i)
      enddo
      call partials (fn_srcpar,
     & nev, ev_cusp, elat, elon, ev_dep,
     & nsta, sta_lab, sta_lat, sta_lon,
     & mod_nl, mod_ratio, mod_v, mod_top,
     & tmp_ttp, tmp_tts,
     & tmp_xp, tmp_yp, tmp_zp)

c--Open synthetic dtime file
      nccp = 0
      nccs = 0
      nctp = 0
      if (iphase.eq.1.or.iphase.eq.3) then
         open(20,file='dtime.P.syn',status='unknown')
         l = 1
         do i=1,nsta
            do j=1,nev
               do k=j+1,nev
                  dt_sta(l) = sta_lab(i)
                  dt_c1(l) = ev_cusp(j)
                  dt_c2(l) = ev_cusp(k)
                  dt_qual(l) = 100
                  dt_idx(l) = 1
                  dt_dt(l) = (tmp_ttp(i,j)-tmp_ttp(i,k))
                  tmp = noisef_dt
                  call ran(-tmp, tmp, dtn)
                  dt_dt(l) = dt_dt(l)+dtn
                  write (20,'(a7,2f15.7,2i4,f9.1,2a)')
     &            dt_sta(l), dt_dt(l), -dt_dt(l),
     &            dt_c1(l), dt_c2(l), dt_qual(l), ' 0'
                  l = l+1
               enddo
            enddo
         enddo
         close(20)
         nccp = l-1
      endif

      if (iphase.eq.2.or.iphase.eq.3) then
         open (21,file='dtime.S.syn',status='unknown')
         l = 1
         do i=1,nsta
            do j=1,nev
               do k=j+1,nev
                  dt_sta(l) = sta_lab(i)
                  dt_c1(l) = ev_cusp(j)
                  dt_c2(l) = ev_cusp(k)
                  dt_qual(l) = 100
                  dt_idx(l) = 2
                  dt_dt(l) = (tmp_tts(i,j)-tmp_tts(i,k))
                  tmp = noisef_dt
                  call ran(-tmp, tmp, dtn)
                  dt_dt(l) = dt_dt(l)+dtn
                  write (21,'(a7,2f15.7,2i4,f9.1,2a)')
     &            dt_sta(l), dt_dt(l), -dt_dt(l),
     &            dt_c1(l), dt_c2(l), dt_qual(l), ' 0'
                  l = l+1
               enddo
            enddo
         enddo
         close(21)
         nccs = l-1
      endif

      write (*,'("# synthetic P dtimes ",i6)') nccp
      write (*,'("# synthetic S dtimes ",i6)') nccs
      write (log,'("# synthetic P dtimes ",i6)') nccp
      write (log,'("# synthetic S dtimes ",i6)') nccs
      if (nccp+nccs.gt.MAXDATA)
     &   stop '>>> Increase MAXDATA in hypoDD.inc'

160   ndt = nccp+nccs+nctp+ncts
      write (*,'("# dtimes total = ",i8)') ndt
      write (log,'("# dtimes total = ",i8)') ndt
      if (ndt.gt.MAXDATA) stop'>>> Increase MAXDATA in hypoDD.inc.'
      if (ndt.eq.0) stop

c--Clean events: dtime match
      do i=1,ndt
         dt_ic1(i) = dt_c1(i)   !dt_ic1 is just a workspace array here!
         dt_ic2(i) = dt_c2(i)   !dt_ic2 is just a workspace array here!
      enddo
      call sorti (ndt, dt_ic1)
      call sorti (ndt, dt_ic2)
      k = 1
      do i=1,nev
         if (ifindi(ndt, dt_ic1, ev_cusp(i)).gt.0) goto 174
         if (ifindi(ndt, dt_ic2, ev_cusp(i)).eq.0) goto 175
174      ev_date(k) = ev_date(i)
         ev_time(k) = ev_time(i)
         ev_lat(k) = ev_lat(i)
         ev_lon(k) = ev_lon(i)
         ev_dep(k) = ev_dep(i)
         ev_mag(k) = ev_mag(i)
         ev_herr(k) = ev_herr(i)
         ev_zerr(k) = ev_zerr(i)
         ev_res(k) = ev_res(i)
         ev_cusp(k) = ev_cusp(i)
         k = k+1
175      continue
      enddo
      nev = k-1
      write (*,'("# events after dtime match = ",i10)') nev
      write (log,'("# events after dtime match = ",i10)') nev

c--New & fast: clean stations
      do i=1,nsta
         sta_itmp(i) = 0
      enddo
      do j=1,ndt
         do i=1,nsta
            if (dt_sta(j).eq.sta_lab(i)) then
               sta_itmp(i) = 1
               goto 176
            endif
         enddo
176      continue
      enddo
      k = 1

      do i=1,nsta
         if (sta_itmp(i).eq.1) then
            sta_lab(k) = sta_lab(i)
            sta_lat(k) = sta_lat(i)
            sta_lon(k) = sta_lon(i)
            k = k+1
            goto 177
         endif
177      continue
      enddo

      nsta = k-1
      write(*,'("# stations = ",i6)') nsta
      write(log,'("# stations = ",i6)') nsta

c--New & fast indexing station labels and cuspids
      call indexxi (nev, ev_cusp, iicusp)
      do i=1,nev
         icusp(i) = ev_cusp(iicusp(i)) !icusp is just a workspace array here!
      enddo
      do i=1,ndt
        do j=1,nsta
            if (dt_sta(i).eq.sta_lab(j)) then
               dt_ista(i) = j
               dt_ic1(i) = iicusp(ifindi(nev, icusp, dt_c1(i)))
               dt_ic2(i) = iicusp(ifindi(nev, icusp, dt_c2(i)))
               goto 200
            endif
        enddo
        stop'FATAL ERROR (indexing). Please report to felix'
200     continue
      enddo
      return

c--Error processing
1010  write (*,*) '** Bad earthquake data, so stop:'
      write (*,*) line
      stop

1051  write (*,'(">>> Format error in cross data file,",/,
     & "OR no origin time corrections ",
     & "available for combined use of cat and cross data.")')
       write (*,*) line
       stop 'Program aborted.'

1091   write(*,*)'>>> Format error in catalog data file.'
       write(*,*)line
       stop 'Program aborted.'


      end  ! of subroutine getdata
