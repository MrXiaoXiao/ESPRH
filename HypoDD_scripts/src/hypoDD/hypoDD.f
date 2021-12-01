      program hypoDD

c Author: Felix Waldhauser, felixw@ldeo.columbia.edu
c Version 1.3 - 11/2010 - FW
c
c started 03/1999 
c 01-03/2001  clean up & bug fixes by Bruce Julian, Fred Klein, Keith
c             Richards-Dinger, Felix Waldhauser 
c 10/2004     Version 1.1: fixed errors listed in BugList to V 1.0.
c 06/2007     Version 1.2 started. fixed errors listed in Buglist 1.1
c 06/2007     accomodate negative magnitudes in output format.
c             real -> doubleprecision: covar,lsfit_svd,matmult2,matmutl3,svd
c 07/2010     fixed bug in computing rms values for hypoDD.reloc (Zhonhe Zhao)
c 07/2010     lsfit_svd: Fixed apparent bug in getting 95% confidence errors 
c             from standard errors. Factor 2.7955 was used, but it should be 
c             1.96, assuming a t-distribution of the residuals (Hilary Martens)
c 07/2010     version 1.2
c 08/2010     now compiles with gfortran  (rcs removed and mod problem fixed)
c 08/2010     version 1.3
c 11/2010     fixed /0 problem in rms reporting (NaN in hypoDD.reloc) (c101116)
c
c Purpose:
c Program to determine high-resolution hypocenter locations using the
c double-difference algorithm. hypoDD incorporates catalog and/or cross
c correlation P- and/or S-wave relative travel-time measurements.
c Residuals between observed and theoretical travel time differences
c (or double-differences = DD) are minimized for pairs
c of earthquakes at each station while linking together all observed
c event/station pairs. A least squares solution (SVD or LSQR) is found
c by iteratively adjusting the vector difference between hypocentral pairs.
c
c References:
c For a detailed description of the algorithm see:
c    Waldhauser, F. and W.L. Ellsworth, A double-difference earthquake
c    location algorithm: Method and application to the northern Hayward
c    fault, Bull. Seismol. Soc. Am., 90, 1353-1368, 2000.
c
c For a user guide to hypoDD see USGS open-file report: 
c    Waldhauser, F., HypoDD: A computer program to compute double-difference
c    earthquake locations,  U.S. Geol. Surv. open-file report , 01-113,
c    Menlo Park, California, 2001.
c
c The code is continuously being updated and improved, so feel
c free to send me an occasional request for the newest version:
c felixw@ldeo.columbia.edu
c or goto http://www.ldeo.columbia.edu/~felixw/hypoDD.html

	implicit none

	include'hypoDD.inc'

	real		acond
	real		adamp(10)
	real		adep
	integer		aiter(0:10)
	real		alat
	real		alon
	real		amaxdcc(10)
	real		amaxdct(10)
	real		amaxres_cross(10)
	real		amaxres_net(10)
	integer		amcusp(1000)
	real		awt_ccp(10)
	real		awt_ccs(10)
	real		awt_ctp(10)
	real		awt_cts(10)
	integer		clust(MAXCL,MAXEVE)
	real		cohav
	real		damp
	character	dattim*25
	real		dtav
	integer		dt_c1(MAXDATA)
	integer		dt_c2(MAXDATA)
	real		dt_cal(MAXDATA)
	real		dt_dt(MAXDATA)
	integer		dt_ic1(MAXDATA)
	integer		dt_ic2(MAXDATA)
	integer		dt_idx(MAXDATA)
	integer		dt_ista(MAXDATA)
	real		dt_offs(MAXDATA)
	real		dt_qual(MAXDATA)
	real		dt_res(MAXDATA)
	character	dt_sta(MAXDATA)*7
	real		dt_wt(MAXDATA)
	real		dxav
	real		dyav
	real		dzav
	real		etav
	integer		ev_cusp(MAXEVE)
	integer		ev_date(MAXEVE)
	real		ev_dep(MAXEVE)
	real		ev_herr(MAXEVE)
	real		ev_lat(MAXEVE)
	real		ev_lon(MAXEVE)
	real		ev_mag(MAXEVE)
	real		ev_res(MAXEVE)
	integer		ev_time(MAXEVE)
	real		ev_x(MAXEVE)
	real		ev_y(MAXEVE)
	real		ev_zerr(MAXEVE)
	real		ev_z(MAXEVE)
	logical		ex
	real		exav
	real		eyav
	real		ezav
	character	fn_cc*80
	character	fn_ct*80
	character	fn_eve*80
	character	fn_inp*80
	character	fn_loc*80
	character	fn_reloc*80
	character	fn_res*80
	character	fn_srcpar*80
	character	fn_sta*80
	character	fn_stares*80
	integer		fu0
	integer		fu1
	integer		fu3
	integer		i
	integer		iargc
	integer		ibeg
	integer		iclust
	integer		icusp(MAXEVE)
	integer		idata
	integer		idy
	integer		iend
	integer		ihr
	integer		imn
	integer		imo
	integer		ineg
	integer		iphase
	integer		isolv
	integer		istart
	integer		iter
	integer		itf
	integer		iunit
	integer		iyr
	integer		j
	integer		jiter
	integer		juliam
	integer		k
	integer		kiter
	integer		l
	doubleprecision	lat
	integer		log
	doubleprecision	lon
	real		maxdcc
	real		maxdct
	real 		maxdist
	integer		maxiter
	real		maxres_cross
	real		maxres_net
	integer		mbad
	integer		minobs_cc
	integer		minobs_ct
	real		minwght
	integer		mod_nl
	real		mod_ratio
	real		mod_top(MAXLAY)
	real		mod_v(MAXLAY)
	integer		n
	integer		narguments
	integer		ncc
	integer		nccold
	integer		nccp
	integer		nccs
	integer		nclust
	integer		nct
	integer		nctold
	integer		nctp
	integer		ncts
	integer		ncusp
	integer		ndt
	integer		nev
	integer		nevold
	integer		niter
	integer		noclust(MAXEVE)
	real		noisef_dt
	integer		nsrc
	integer		nsta
	real		picav
	real		resvar1
	real		rms_cc
	real		rms_cc0
	real		rms_cc0old
	real		rms_ccold
	real		rms_ct
	real		rms_ct0
	real		rms_ct0old
	real		rms_ctold
	real		sc
	real		sdc0_dep
	real		sdc0_lat
	real		sdc0_lon
	integer		src_cusp(MAXEVE)
	real		src_dep(MAXEVE)
	real		src_dt(MAXEVE)
	real		src_dx(MAXEVE)
	real		src_dy(MAXEVE)
	real		src_dz(MAXEVE)
	real		src_et(MAXEVE)
	real		src_ex(MAXEVE)
	real		src_ey(MAXEVE)
	real		src_ez(MAXEVE)
	real		src_lat0(MAXEVE)
	doubleprecision	src_lat(MAXEVE)
	real		src_lon0(MAXEVE)
	doubleprecision	src_lon(MAXEVE)
	integer		src_nnp(MAXEVE)
	integer		src_nns(MAXEVE)
	integer		src_np(MAXEVE)
	integer		src_ns(MAXEVE)
	real		src_rmsc(MAXEVE)
	real		src_rmsn(MAXEVE)
	real		src_t0(MAXEVE)
	real		src_t(MAXEVE)
	real		src_x0(MAXEVE)
	real		src_x(MAXEVE)
	real		src_y0(MAXEVE)
	real		src_y(MAXEVE)
	real		src_z0(MAXEVE)
	real		src_z(MAXEVE)
	real		sta_az(MAXSTA)
	real		sta_dist(MAXSTA)
	character	sta_lab(MAXSTA)*7
	real		sta_lat(MAXSTA)
	real		sta_lon(MAXSTA)
	integer		sta_nnp(MAXSTA)
	integer		sta_nns(MAXSTA)
	integer		sta_np(MAXSTA)
	integer		sta_ns(MAXSTA)
	real		sta_rmsc(MAXSTA)
	real		sta_rmsn(MAXSTA)
	character	str1*60
	character	str80*80
	character	str3*3
	real		tav
	real		tmpr1
	real		tmpr2
	real		tmp_ttp(MAXSTA,MAXEVE)
	real		tmp_tts(MAXSTA,MAXEVE)
	real		tmp_xp(MAXSTA,MAXEVE)
	real		tmp_yp(MAXSTA,MAXEVE)
	real		tmp_zp(MAXSTA,MAXEVE)
	integer		trimlen
	real		wt_ccp
	real		wt_ccs
	real		wt_ctp
	real		wt_cts
	real		x
	real		xav
	real		y
	real		yav
	real		zav

      minwght= 0.00001
      rms_ccold= 0
      rms_ctold= 0
      rms_cc0old= 0
      rms_ct0old=  0
c--- open log file:
      call freeunit(log)
      open(log,file='hypoDD.log',status='unknown')
      str1= 'starting hypoDD (v1.3 - 11/2010)...'
      call datetime(dattim)
      write(6,'(a45,a)') str1, dattim
      write(log,'(a45,a)') str1, dattim

c--- get input parameter file name:
      narguments = iargc()
      if(narguments.lt.1) then
        write(*,'(/,a)') 'PARAMETER INPUT FILE [hypoDD.inp <ret>]:'
        read(5,'(a)') fn_inp
        if(trimlen(fn_inp).le.1) then
           fn_inp= 'hypoDD.inp'            !default input file name
        else
           fn_inp= fn_inp(1:trimlen(fn_inp))
        endif
      else
          call getarg(1,fn_inp)
      endif
      inquire(FILE= fn_inp,exist=ex)
      if(.not. ex) stop' >>> ERROR OPENING INPUT PARAMETER FILE.'

c--- get input parameters:
      call getinp(MAXEVE,MAXLAY,log,fn_inp,
     & fn_cc,fn_ct,fn_sta,fn_eve,
     & fn_loc,fn_reloc,fn_res,fn_stares,fn_srcpar,
     & idata,iphase,
     & minobs_cc,minobs_ct,
     & amaxres_cross,amaxres_net,amaxdcc,amaxdct,
     & noisef_dt,maxdist,
     & awt_ccp,awt_ccs,awt_ctp,awt_cts,adamp,
     & istart,maxiter,isolv,niter,aiter,
     & mod_nl,mod_ratio,mod_v,mod_top,
     & iclust,ncusp,icusp)

c--- get data:
      call getdata(
     & log,fn_cc,fn_ct,fn_sta,fn_eve,fn_srcpar,
     & idata,iphase,ncusp,icusp,
     & maxdist,amaxdct(1),amaxdcc(1),
     & noisef_dt,mod_nl,mod_ratio,mod_v,mod_top,
     & ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,
     & ev_mag,ev_herr,ev_zerr,ev_res,
     & sta_lab,sta_lat,sta_lon,
     & dt_sta,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
     & dt_ista,dt_ic1,dt_ic2,dt_offs,
     & nev,nsta,ndt,nccp,nccs,nctp,ncts,
     & tmp_xp,tmp_yp,tmp_zp,tmp_ttp,tmp_tts)

c--- clustering:
      if((idata.eq.1.and.minobs_cc.eq.0).or.
     &   (idata.eq.2.and.minobs_ct.eq.0).or.
     &   (idata.eq.3.and.minobs_ct+minobs_cc.eq.0)) then
         nclust= 1
         clust(1,1)= nev
         do i=1,nev
             clust(1,i+1)= ev_cusp(i)
         enddo
          write(*,'(/,"no clustering performed.")')
          write(log,'(/,"no clustering performed.")')
      else

         call cluster1(log, nev, ndt,
     & idata, minobs_cc, minobs_ct,
     & dt_c1, dt_c2, ev_cusp,
     & clust, noclust, nclust)

      endif

c--- open files
      call freeunit(fu0)
      open(fu0,file=fn_loc,status='unknown')
      call freeunit(fu1)
      open(fu1,file=fn_reloc,status='unknown')
      if(trimlen(fn_stares).gt.1) then
         call freeunit(fu3)
         open(fu3,file=fn_stares,status='unknown')
      endif

      jiter = 0  ! counter for iter with no updating (air quakes)
c--- big loop over clusters starts here:
      if(iclust.ne.0) then
        if (iclust.lt.0 .or. iclust.gt.nclust) then
           write(*,*) 'error: invalid cluster number ',iclust
           write(*,*) 'must be between 1 and nclust (',nclust,')'
           stop
        endif
        ibeg= iclust
        iend= iclust
      else
        ibeg= 1
        iend= nclust
      endif
      do iclust= ibeg,iend
      call datetime(dattim)
      write(log,'(/,"RELOCATION OF CLUSTER:",i2,5x,a,/,
     &"----------------------")')iclust,dattim
      write(*,'(/,"RELOCATION OF CLUSTER:",i2,5x,a,/,
     &"----------------------")')iclust,dattim

c--- get data for each cluster if clustering was invoked:
      if((nclust.eq.1.and.clust(iclust,1).eq.nev).or.
     &   (idata.eq.1.and.minobs_cc.eq.0).or.
     &   (idata.eq.2.and.minobs_ct.eq.0).or.
     &   (idata.eq.3.and.minobs_ct+minobs_cc.eq.0)) goto 50

      ncusp= clust(iclust,1)
      do i=1,ncusp
         icusp(i)= clust(iclust,i+1)
      enddo

      if(idata.ne.0) call getdata(
     & log,fn_cc,fn_ct,fn_sta,fn_eve,fn_srcpar,
     & idata,iphase,ncusp,icusp,
     & maxdist,amaxdct(1),amaxdcc(1),
     & noisef_dt,mod_nl,mod_ratio,mod_v,mod_top,
     & ev_date,ev_time,ev_cusp,ev_lat,ev_lon,ev_dep,
     & ev_mag,ev_herr,ev_zerr,ev_res,
     & sta_lab,sta_lat,sta_lon,
     & dt_sta,dt_dt,dt_qual,dt_c1,dt_c2,dt_idx,
     & dt_ista,dt_ic1,dt_ic2,dt_offs,
     & nev,nsta,ndt,nccp,nccs,nctp,ncts,
     & tmp_xp,tmp_yp,tmp_zp,tmp_ttp,tmp_tts)

50    continue
      nccold= nccp+nccs
      nctold= nctp+ncts
      ncc= nccp+nccs
      nct= nctp+ncts
      nevold= nev

c--- get cluster centroid:
      sdc0_lat= 0
      sdc0_lon= 0
      sdc0_dep= 0
      do i=1,nev
         sdc0_lat= sdc0_lat + ev_lat(i)
         sdc0_lon= sdc0_lon + ev_lon(i)
         sdc0_dep= sdc0_dep + ev_dep(i)
      enddo
      sdc0_lat= sdc0_lat/nev
      sdc0_lon= sdc0_lon/nev
      sdc0_dep= sdc0_dep/nev

      write(log,'("Cluster centroid at:",1x,f10.6,2x,f11.6,2x,f9.6)')
     & sdc0_lat,sdc0_lon,sdc0_dep

c--- Set up cartesian coordinate system (build common block for subr. SDC):
      call setorg(sdc0_lat,sdc0_lon,0.0,0)

c--- get cartesian coordinates for epicenters
      do i=1,nev
         lat= ev_lat(i)
         lon= ev_lon(i)
         call SDC2(x,y,lat,lon,-1)
         ev_x(i)= x *1000
         ev_y(i)= y *1000
         ev_z(i)= (ev_dep(i) - sdc0_dep)*1000
      enddo

      write(log,'("# events:",i5)')nev

c--- write output (mdat.loc):
      write(fu0,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     & 1x,f10.1,
     & 1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f5.2,
     & 1x,f4.1,1x,i3)')
cfw     & 1x,f3.1,1x,i3)')	neg mag format
     & (ev_cusp(i),ev_lat(i),ev_lon(i),ev_dep(i),ev_x(i),ev_y(i),
     & ev_z(i),ev_herr(i)*1000,ev_herr(i)*1000,ev_zerr(i)*1000,
     & int(ev_date(i)/10000),
     & int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     & int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     & mod(real(ev_time(i)),10000.)/100,ev_mag(i),iclust,i=1,nev)
cfw100806     & mod(real(ev_time(i)),10000)/100,ev_mag(i),iclust,i=1,nev)

c--- get initial trial sources:
      call trialsrc(istart,sdc0_lat,sdc0_lon,sdc0_dep,
     & nev,ev_cusp,ev_lat,ev_lon,ev_dep,
     & nsrc,src_cusp,src_lat0,src_lon0,
     & src_x0,src_y0,src_z0,src_t0,
     & src_lat,src_lon,src_dep,
     & src_x,src_y,src_z,src_t)

      write(*,'("Initial trial sources =",i6)')nsrc
      write(log,'("# initial trial sources:",i6)')nsrc

c--- loop over iterations starts here:
c define iteration step at which re-weighting starts: this is dynam. since
c it depends on the number of neg depths runs before.

c first reset aiter() and maxiter
      do i=1,niter
         aiter(i) = aiter(i) - jiter
      enddo
      maxiter = maxiter - jiter

      kiter= 0		! counter for iter with data skipping
      jiter= 0		! counter for iter with no updating (air quakes)
      mbad= 0		! counter for air quakes

      iter= 1
55    call datetime(dattim)
      write(log,'(/,"===ITERATION ",i3," (",i3,") ",a)')
     & iter-jiter, iter, dattim

c--- get weighting parameters for this iteration:
      do i=1,niter
        if(iter.le.aiter(i)) goto 75
      enddo
75    maxres_cross= amaxres_cross(i)
      maxres_net= amaxres_net(i)
      maxdcc= amaxdcc(i)
      maxdct= amaxdct(i)
      wt_ccp= awt_ccp(i)
      wt_ccs= awt_ccs(i)
      wt_ctp= awt_ctp(i)
      wt_cts= awt_cts(i)
      damp= adamp(i)

      write(log, '(/,"Weighting parameters for this iteration:",/,
     &"  wt_ccp= ",f7.4,2X,"wt_ccs= ",f7.4,2X,
     &"maxr_cc= ",f7.4,2X,"maxd_cc= ",f7.2,2X,/,
     &"  wt_ctp= ",f7.4,2x,"wt_cts= ",f7.4,2x,"maxr_ct= ",f7.4,2x,
     &"maxd_ct= ",f7.2,/,"  damp= ",f5.1)')
     & wt_ccp,wt_ccs,maxres_cross,
     & maxdcc,wt_ctp,wt_cts,maxres_net,
     & maxdct,damp

c--- calculate travel times  and slowness vectors:
      write(log,'(/,"~ getting partials for ",i5,
     & " stations and ",i5," source(s) ...")') nsta,nsrc
      call partials(fn_srcpar,
     & nsrc,src_cusp,src_lat,src_lon,src_dep,
     & nsta,sta_lab,sta_lat,sta_lon,
     & mod_nl,mod_ratio,mod_v,mod_top,
     & tmp_ttp,tmp_tts,
     & tmp_xp,tmp_yp,tmp_zp)

c--- get double difference vector:
      call dtres(log,ndt,MAXSTA,nsrc,
     & dt_dt,dt_idx,
     & dt_ista,dt_ic1,dt_ic2,
     & src_cusp,src_t,tmp_ttp,tmp_tts,
     & dt_cal,dt_res)

c--- get a priori weights and reweight residuals
      call weighting(log,ndt,mbad,amcusp,idata,kiter,ineg,
     &               maxres_cross,maxres_net,maxdcc,maxdct,minwght,
     &               wt_ccp,wt_ccs,wt_ctp,wt_cts,
     &               dt_c1,dt_c2,dt_idx,dt_qual,dt_res,dt_offs,
     &               dt_wt)

c--- skip outliers and/or air quakes:
      if(ineg.gt.0) then
          call skip(log,kiter,minwght,
     & ndt,nev,nsrc,nsta,
     & ev_cusp,ev_date,ev_time,ev_mag,
     & ev_lat,ev_lon,ev_dep,ev_x,ev_y,ev_z,
     & ev_herr,ev_zerr,ev_res,
     & src_cusp, src_lat, src_lon, src_dep,
     & src_lat0, src_lon0,
     & src_x, src_y, src_z, src_t, src_x0, src_y0, src_z0, src_t0,
     & sta_lab,sta_lat,sta_lon,sta_dist,sta_az,
     & sta_rmsc,sta_rmsn,sta_np,sta_ns,sta_nnp,sta_nns,
     & dt_sta,dt_c1,dt_c2,dt_idx,dt_dt,dt_qual,dt_cal,
     & dt_ista,dt_ic1,dt_ic2,
     & dt_res,dt_wt,dt_offs,
     & tmp_ttp,tmp_tts,tmp_xp,tmp_yp,tmp_zp,nct,ncc)

c--Dont mess anymore with this cluster if we have wiped out all events
        if(nev.lt.2) then 
          write(log,*)' Cluster has less than 2 events.'
          write(*,*)' Cluster has less than 2 events.'
          goto 778
        endif
      else
         write(log,'("no data skipped.")')

      endif

c--- get initial residual statistics (avrg,rms,var..)
      if(iter.eq.1) then
       resvar1= -999
       call resstat(log,idata,ndt,nev,dt_res,dt_wt,dt_idx,
     & rms_cc,rms_ct,rms_cc0,rms_ct0,
     & rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     &              resvar1)
      endif

c--- least square fitting:

      if(isolv.eq.1) then
         call lsfit_SVD(log,iter,ndt,nev,nsrc,damp,mod_ratio,
     & idata,ev_cusp,src_cusp,
     & dt_res,dt_wt,
     & dt_ista,dt_ic1,dt_ic2,   !new
     & src_dx,src_dy,src_dz,src_dt,src_ex,src_ey,src_ez,src_et,
     & exav,eyav,ezav,etav,dxav,dyav,dzav,dtav,
     & rms_cc,rms_ct,rms_cc0,rms_ct0,
     & rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     & tmp_xp,tmp_yp,tmp_zp,dt_idx)

      else
         call lsfit_lsqr(log,iter,ndt,nev,nsrc,damp,mod_ratio,
     & idata,ev_cusp,src_cusp,
     & dt_res,dt_wt,
     & dt_ista,dt_ic1,dt_ic2,   !new
     & src_dx,src_dy,src_dz,src_dt,src_ex,src_ey,src_ez,src_et,
     & exav,eyav,ezav,etav,dxav,dyav,dzav,dtav,
     & rms_cc,rms_ct,rms_cc0,rms_ct0,
     & rms_ccold,rms_ctold,rms_cc0old,rms_ct0old,
     & tmp_xp,tmp_yp,tmp_zp,dt_idx,acond)
      endif

c--- check for air quakes:
      mbad= 0
      k= 1
      do i= 1,nsrc
        if(src_dep(i) + (src_dz(i)/1000).lt.0) then
            write(log,'(">>>Warning: negative depth - ",i12)')ev_cusp(i)
            amcusp(k)= ev_cusp(i)
            k=k+1
            if(k.gt.1000) stop'>>> More than 1000 air quakes. Too many!'
        endif
      enddo
      mbad= k-1    ! number of neg depth events

c update iteration numbers:
      if(mbad.gt.0) then
         do i= 1,niter
              aiter(i)= aiter(i)+1
         enddo
         jiter= jiter+1	  ! iteration with no update
         maxiter= maxiter+1

         write(log,*)'Number of air quakes (AQ) =',mbad
         if(nsrc-mbad .le. 1) then
           write(*,*)'Warning: number of non-airquakes < 2'
                   write(*,*)'   skipping this cluster'
           write(log,*)'Warning: number of non-airquakes < 2'
                   write(log,*)'   skipping this cluster'
                   goto 778
         endif       
         goto 500   ! skip the updating step
      endif


c--- update source parameters:
      xav= 0 ! mean centroid shift
      yav= 0
      zav= 0
      tav= 0
      alon= 0
      alat= 0
      adep= 0
      if(nsrc.eq.1) nsrc= nev
      do i= 1,nsrc
        src_cusp(i)= ev_cusp(i)
c update absolute source parameters (cart)
        src_x(i)= src_x(i) + src_dx(i)
        src_y(i)= src_y(i) + src_dy(i)
        src_z(i)= src_z(i) + src_dz(i)
        src_t(i)= src_t(i) + src_dt(i)

c update absolute source locations (geogr)
        src_dep(i)= src_dep(i) + (src_dz(i)/1000)
        call SDC2(src_x(i)/1000,src_y(i)/1000,lat,lon,1)
        src_lon(i)= lon
        src_lat(i)= lat
        alon= lon+alon	
        alat= lat+alat
        adep= adep+src_dep(i)

c get mean centroid shift
        xav= xav + (src_x(i) - src_x0(i))	
        yav= yav + (src_y(i) - src_y0(i))
        zav= zav + (src_z(i) - src_z0(i))
        tav= tav + (src_t(i) - src_t0(i))
      enddo
      xav= xav/nsrc
      yav= yav/nsrc
      zav= zav/nsrc
      tav= tav/nsrc
      alon= alon/nsrc
      alat= alat/nsrc
      adep= adep/nsrc

      write(log,'("  cluster centroid at:",1x,f10.6,2x,f11.6,2x,f9.6)')
     & alat,alon,adep
      write(log,'("  mean centroid (origin) shift in x,y,z,t [m,ms]: ",/
     & f7.1,f7.1,f7.1,f7.1)'),xav,yav,zav,tav
      write(log,'("  (OS in std output gives maximum value.)")')

c--- get interevent distance for each observation and average signal coherency:
      cohav= 0
      picav= 0
      j= nct
      k= ncc
      ncc= 0
      nct= 0
      do i= 1,ndt
         dt_offs(i)= sqrt((src_x(dt_ic1(i))-src_x(dt_ic2(i)))**2 +
     &                  (src_y(dt_ic1(i))-src_y(dt_ic2(i)))**2 +
     &                  (src_z(dt_ic1(i))-src_z(dt_ic2(i)))**2)

         if(dt_idx(i).le.2) then
            cohav= cohav + sqrt(dt_qual(i))
            ncc= ncc+1
         else
            picav= picav + dt_qual(i)
            nct= nct+1
         endif

      enddo
      cohav= cohav/ncc
      picav= picav/nct
      write(log,'(/,"More:")')
      write(log,'("  mean phase coherency = ",f5.3)')cohav
      write(log,'("  mean pick quality = ",f5.3)')picav

c--- get number of observations and mean residual at each station
      tmpr1= 0
      tmpr2= 0
      do i= 1,nsta
         sta_np(i)= 0
         sta_ns(i)= 0
         sta_nnp(i)= 0
         sta_nns(i)= 0
         sta_rmsc(i)= 0
         sta_rmsn(i)= 0
         do j= 1,ndt
            if(i.eq.dt_ista(j)) then
               if(dt_idx(j).le.2) then
                 sta_rmsc(i)= sta_rmsc(i)+dt_res(j)**2
                 if(dt_idx(j).eq.1) then
                    sta_np(i)= sta_np(i)+1
                 else
                    sta_ns(i)= sta_ns(i)+1
                 endif
               else
                 sta_rmsn(i)= sta_rmsn(i)+dt_res(j)**2
                 if(dt_idx(j).eq.3) then
                   sta_nnp(i)= sta_nnp(i)+1
                 else
                   sta_nns(i)= sta_nns(i)+1
                 endif
               endif
            endif
         enddo

         if(sta_np(i)+sta_ns(i).gt.0)
     &     sta_rmsc(i)= sqrt(sta_rmsc(i)/(sta_np(i)+sta_ns(i)))
         if(sta_nnp(i)+sta_nns(i).gt.0)
     &     sta_rmsn(i)= sqrt(sta_rmsn(i)/(sta_nnp(i)+sta_nns(i)))
         if(sta_rmsc(i).gt.tmpr1) then
            tmpr1= sta_rmsc(i)
            k= i
         endif
         if(sta_rmsn(i).gt.tmpr2) then
            tmpr2= sta_rmsn(i)
            l= i
         endif
      enddo
      tmpr1= tmpr1*1000
      tmpr2= tmpr2*1000
      if(idata.eq.1.or.idata.eq.3) then
         write(log,'("  station with largest cc rms: ",a7,"=",
     & f7.0," ms (RMSST)")')
     &    sta_lab(k),tmpr1
      endif
      if(idata.eq.2.or.idata.eq.3) then
         write(log,'("  station with largest ct rms: ",a7,"=",
     & f7.0," ms (RMSST)")')
     &    sta_lab(l),tmpr2
      endif

c--- write output scratch mdat.reloc:
      n= trimlen(fn_reloc)
      i=iter-jiter
      write(str80,'(a,".",i3.3,".",i3.3)')fn_reloc(1:n),iclust,i
      call freeunit(iunit)
      open(iunit,file=str80,status='unknown')
      write(iunit,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     & 1x,f10.1,
     & 1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,
     & 1x,f4.1,1x,i3)')
cfw     & 1x,f3.1,1x,i3)')		! negative mag
     & (src_cusp(i),src_lat(i),src_lon(i),src_dep(i),src_x(i),src_y(i),
     & src_z(i),src_ex(i),src_ey(i),src_ez(i),int(ev_date(i)/10000),
     & int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     & int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     & mod(real(ev_time(i)),10000.)/100,ev_mag(i),iclust,i=1,nev)
cfw100806     & mod(real(ev_time(i)),10000)/100,ev_mag(i),iclust,i=1,nev)
      close(iunit)
c WARNING: variable "str" is set to zero value by default
      write(log,'(/,"Relocation results for this iteration are"
     & " stored in ",a)')str80(1:trimlen(str80))

500   continue  ! case of air quakes

c standard output:
      if(mbad.gt.0) then
         str3='   '
      else
         n= iter-jiter
         if(n.lt.1000) write(str3,'(i3)')n
         if(n.lt.100) write(str3,'(1x,i2)')n
         if(n.lt.10) write(str3,'(2x,i1)')n
      endif
      if(isolv.eq.1.and.idata.eq.3) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CT  CC",
     & "    RMSCT      RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ",/,
     &                           "        %   %   %",
     & "   ms     %   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,3(1x,i3),i5,f6.1,i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),nint(nct*100.0/nctold),
     & nint(ncc*100.0/nccold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),mbad
      endif
      if(isolv.eq.1.and.idata.eq.1) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CC",
     & "    RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ",/,
     &                          "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(ncc*100.0/nccold),
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),mbad
      endif
      if(isolv.eq.1.and.idata.eq.2) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CT",
     & "    RMSCT     RST   DX   DY   DZ   DT   OS  AQ",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(nct*100.0/nctold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),mbad
      endif

      if(isolv.eq.2.and.idata.eq.3) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CT  CC",
     & "    RMSCT      RMSCC   RMSST   DX   DY   DZ   DT   ",
     & "OS  AQ  CND",/,
     &                           "        %   %   %",
     & "   ms     %   ms     %    ms    m    m    m   ms   ",
     & " m     ")')
       write(*,'(i2,a3,3(1x,i3),i5,f6.1,i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),nint(nct*100.0/nctold),
     & nint(ncc*100.0/nccold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),mbad,nint(acond)
      endif
      if(isolv.eq.2.and.idata.eq.1) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CC",
     & "    RMSCC   RMSST   DX   DY   DZ   DT   OS  AQ  CND",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(ncc*100.0/nccold),
     & nint(rms_cc*1000),(rms_cc-rms_ccold)*100/rms_ccold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),mbad,nint(acond)
      endif
      if(isolv.eq.2.and.idata.eq.2) then
       if(iter.eq.1) write(*,'(/,"  IT   EV  CT",
     & "    RMSCT   RMSST   DX   DY   DZ   DT   OS  AQ  CND",/,
     &                           "        %   %",
     & "   ms     %    ms    m    m    m   ms    m ")')
       write(*,'(i2,a3,2(1x,i3),i5,f6.1,i6,4i5,i5,i4,i5)')
     & iter,str3,
     & nint(nev*100./nevold),
     & nint(nct*100.0/nctold),
     & nint(rms_ct*1000),(rms_ct-rms_ctold)*100/rms_ctold,
     & nint(max(tmpr1,tmpr2)),
     & nint(dxav),nint(dyav),nint(dzav),nint(dtav),
     & nint(max(abs(xav),abs(yav),abs(zav))),mbad,nint(acond)
      endif

      call datetime(dattim)
      write(log,'("Iteration ",i2," finished ",a)') iter, dattim

      if(iter.eq.maxiter) goto 600	! all iterations done.
      iter= iter+1
      goto 55	! next iteration

c--- update origin time (this is only done for final output!!)
600   continue
      write(*,'(/,"writing out results ...")')
      do i= 1,nev
         src_t(i)= src_t(i)/1000	!from here on src_t in sec!!
         if(src_t(i).gt.5) then
            write(*,*)'WARNING: org time diff > 5s for ',src_cusp(i)
         endif
         iyr= int(ev_date(i)/10000)
         imo= int(mod(ev_date(i),10000)/100)
         idy= int(mod(ev_date(i),100))
         ihr= int(ev_time(i)/1000000)
         imn= int(mod(ev_time(i),1000000)/10000)
         itf= JULIAM(iyr,imo,idy,ihr,imn)

cfw         sc= (mod(real(ev_time(i)),10000)/100) + src_t(i)
         sc= (mod(real(ev_time(i)),10000.)/100) - src_t(i)
cfw100806         sc= (mod(real(ev_time(i)),10000)/100) - src_t(i)
         itf= itf + int(sc / 60.)
         sc=  sc  - int(sc / 60.)*60.
         if(sc.lt.0) then
            itf= itf-1
            sc= 60. + sc
         endif
         call DATUM(itf,iyr,imo,idy,ihr,imn)
         ev_date(i)= iyr*10000 + imo*100 + idy
         ev_time(i)= ihr*1000000 + imn*10000 + nint(sc*100)
      enddo

c--- get # of obs per event:
      do i=1,nev
         src_np(i)= 0
         src_ns(i)= 0
         src_nnp(i)= 0
         src_nns(i)= 0
         src_rmsc(i)= 0
         src_rmsn(i)= 0
      enddo
      do i=1,ndt
         if(dt_idx(i).eq.1) then
             src_np(dt_ic1(i))= src_np(dt_ic1(i))+1
             src_np(dt_ic2(i))= src_np(dt_ic2(i))+1
         endif
         if(dt_idx(i).eq.2) then
             src_ns(dt_ic1(i))= src_ns(dt_ic1(i))+1
             src_ns(dt_ic2(i))= src_ns(dt_ic2(i))+1
         endif
         if(dt_idx(i).le.2) then
             src_rmsc(dt_ic1(i))= src_rmsc(dt_ic1(i))+dt_res(i)**2
             src_rmsc(dt_ic2(i))= src_rmsc(dt_ic2(i))+dt_res(i)**2
         endif
         if(dt_idx(i).eq.3) then
             src_nnp(dt_ic1(i))= src_nnp(dt_ic1(i))+1
             src_nnp(dt_ic2(i))= src_nnp(dt_ic2(i))+1
         endif
         if(dt_idx(i).eq.4) then
             src_nns(dt_ic1(i))= src_nns(dt_ic1(i))+1
             src_nns(dt_ic2(i))= src_nns(dt_ic2(i))+1
         endif
         if(dt_idx(i).ge.3) then
             src_rmsn(dt_ic1(i))= src_rmsn(dt_ic1(i))+dt_res(i)**2
             src_rmsn(dt_ic2(i))= src_rmsn(dt_ic2(i))+dt_res(i)**2
         endif
      enddo
      do i=1,nev
c100710         src_rmsc(i)= sqrt(src_rmsc(i)/nev)
c100710         src_rmsn(i)= sqrt(src_rmsn(i)/nev)
c101116         src_rmsc(i)= sqrt(src_rmsc(i)/(src_np(i)+src_ns(i)))
c101116         src_rmsn(i)= sqrt(src_rmsn(i)/(src_nnp(i)+src_nns(i)))
         if(src_np(i)+src_ns(i).gt.0) then
            src_rmsc(i)= sqrt(src_rmsc(i)/(src_np(i)+src_ns(i)))
         else
            src_rmsc(i)= -9 
         endif
         if(src_nnp(i)+src_nns(i).gt.0) then
            src_rmsn(i)= sqrt(src_rmsn(i)/(src_nnp(i)+src_nns(i)))
         else
            src_rmsn(i)= -9 
         endif
      enddo

c--- output final residuals: mdat.res
      if(trimlen(fn_res).gt.1) then
         call freeunit(iunit)
         open(iunit,file=fn_res,status='unknown')
         write(iunit,'("STA",11x,"DT",8x,
     &"C1",8x,"C2",4x,"IDX",5x,"QUAL",4x,"RES [ms]",3x,"WT",9x,
     &"OFFS")')
         write(iunit,'(a7,1x,f12.7,1x,i9,1x,i9,1x,i1,1x,
     & f9.4,1x,f12.6,1x,f11.6,1x,f8.1)')
     & (dt_sta(j),dt_dt(j),dt_c1(j),dt_c2(j),dt_idx(j),dt_qual(j),
     & dt_res(j)*1000,dt_wt(j),dt_offs(j),j=1,ndt)
         close(iunit)
      endif

c--- output final locations (mdat.reloc):
      write(fu1,'(i9,1x,f10.6,1x,f11.6,1x,f9.3,1x,f10.1,1x,f10.1,
     & 1x,f10.1,
     & 1x,f8.1,1x,f8.1,1x,f8.1,1x,i4,1x,i2,1x,i2,1x,i2,1x,i2,1x,f6.3,
     & 1x,f4.1,1x,i5,1x,i5,1x,i5,1x,i5,1x,f6.3,1x,f6.3,1x,i3)')
cfw     & 1x,f3.1,1x,i5,1x,i5,1x,i5,1x,i5,1x,f6.3,1x,f6.3,1x,i3)') !neg mag
     & (src_cusp(i),src_lat(i),src_lon(i),src_dep(i),src_x(i),src_y(i),
     & src_z(i),src_ex(i),src_ey(i),src_ez(i),int(ev_date(i)/10000),
     & int(mod(ev_date(i),10000)/100),mod(ev_date(i),100),
     & int(ev_time(i)/1000000),int(mod(ev_time(i),1000000)/10000),
     & mod(real(ev_time(i)),10000.)/100,ev_mag(i),
cfw100806     & mod(real(ev_time(i)),10000)/100,ev_mag(i),
     & src_np(i),src_ns(i),src_nnp(i),src_nns(i),
     & src_rmsc(i),src_rmsn(i), iclust,i=1,nev)

c--- output stations (mdat.station):
      if(trimlen(fn_stares).gt.1) then
cfw         write(fu3,'(a5,1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4,1x,i7,1x,
         write(fu3,'(a7,1x,f9.4,1x,f9.4,1x,f9.4,1x,f9.4,1x,i7,1x,
     & i7,1x,i7,1x,i7,1x,f9.4,1x,f9.4,1x,i3)')
     & (sta_lab(i),sta_lat(i),sta_lon(i),sta_dist(i),sta_az(i),
     & sta_np(i),sta_ns(i),sta_nnp(i),sta_nns(i),
     & sta_rmsc(i),sta_rmsn(i),iclust,i=1,nsta)
      endif

778   continue
      enddo  ! loop over clusters (iclust)

      close(fu0)
      if(trimlen(fn_stares).gt.1) close(fu3)
      close(fu1)

      end !of main routine
