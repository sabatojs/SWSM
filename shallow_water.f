      program shallow_water

      implicit none

      integer N,M,trunc,itype,Tf,rhom
      parameter(N=128,M=64,trunc=42,rhom=0,Tf=17280*2)
		! spectral truncation/lat-lon resolution and 
		! final model time step (+1 for t=0, i.e. complete Tf timesteps)
		! spectral transform after Durran

      integer i,j,k,tt,recno,wrtflg,howoften,pwr2,fac3

      !!!! real variables
      double precision u(N,M),v(N,M),h(N,M),z(N,M),d(N,M),hb(N,M)
      double precision y(M),x(N), lats(M), lons(N)
      double precision gauzer(M),gauwts(M)
      double precision eps,fac2,fac4,fac5

      double precision e(N,M),vz(N,M),uz(N,M)
      double precision uphi(N,M),vphi(N,M)
      double precision muz(N,M),v0(N,M),u0(N,M),eta(N,M),etab(N,M)
      double precision mu(M),jjp1,jjm1,jp1jp2,fac1

      double precision H0,f,dt,t,rob,nue,nunm,tau,kf
      double precision A,g,OMEGA,nrg,nu
      integer randseed

      !!!! spectral variables
      double complex esp(N,M),usp(N,M),vsp(N,M)
      double complex phi(N,M,3),psi(N,M,3),chi(N,M,3)
      double complex phib(N,M)
      double complex gba(N,M),gab(N,M),gcd(N,M)
      double complex zsp(N,M),etasp(N,M),etabsp(N,M),dsp(N,M)
      double complex A1(N,M),B1(N,M),C1(N,M)

      logical dofilter,dofrctn,donrg,data_init

      character(17) FILE_OUT
      character(20) FILE_IN

      double complex im,pie

      integer fu,rc
      character(13) namelistfile
      namelist /RUNPARAMS/ A,g,OMEGA,H0,dofrctn,tau,donrg,nrg,nu
      namelistfile = 'runparams.nml'

      im = dcmplx(0.d0,1.d0)
      pie = 4.d0*datan(1.d0)

      FILE_OUT = 'shallow_water.cdf'	! name of output file

      dofilter = .true.			! perform robert filter?
      data_init = .false.		! initialize with some dataset?

      open(action='read', file=namelistfile, iostat=rc, unit=fu)
      read(unit=fu,iostat=rc,nml=RUNPARAMS)
      close(fu)

      write(*,*) ' '
      write(*,*) ' '
      write(*,*) '          PLANET RADIUS: ', A ,' m'
      write(*,*) '   PLANET ROTATION RATE: ', OMEGA ,' 1/s'
      write(*,*) 'ACCELERATION OF GRAVITY: ', g, ' m/s2'
      write(*,*) '       MEAN FLUID DEPTH: ', H0, ' m'
      write(*,*) '  FRICTION DAMPING TIME: ', tau, ' days'
      write(*,*) ' '
      write(*,*) ' '

      !!! choose file name, re-scale CFL time-step/hyper-diffusion
      if (trunc.eq.85) then
	 FILE_IN = 'initial_data_T85.cdf'
	 fac2 = 4.d0
         fac3 = 4
      else if (trunc.eq.42) then
	 FILE_IN = 'initial_data_T42.cdf'
         fac2 = 2.d0
         fac3 = 2
      else
	 FILE_IN = 'initial_data_T21.cdf'
         fac2 = 1.d0
         fac3 = 1
      endif
      !!! large-scale friction?
      if (dofrctn) then
	 fac4 = 1.d0
      else
         fac4 = 0.d0
      endif
      !!! random energy production?
      if (donrg) then
	 fac5 = 1.d0
      else
         fac5 = 0.d0
      endif

      !!! namelist parameters !!!
      !A = 6.378d6			! radius of planet
      !g = 9.81d0			! acceleration due to gravity
      !f = 2.d0*7.292d-5		! coriolis parameter = twice rotation rate in radians per second
      !H0 = 4000.d0			! mean fluid depth
      !dt = 900.d0*2.d0/fac2		! time step 
      !tau =  30.d0 			! frictional time scale in days

      !!! derived parameters
      f = 2.d0*OMEGA		        ! coriolis parameter
      kf = fac4/(tau*86400.d0)		! frictional damping rate in s^-1

      !!! numerical parameters
      rob = 0.05d0			! robert filter
      nue = nu/((fac2/2.d0)**3)  	! sub-grid scale diffusion
      dt = nint((A/(35.774824021263     ! scale CFL time-step for resolution and parameters
     +     *dsqrt(g*H0)))*2.d0 / fac2) 	! default 900s time-step for T42
				      	! for g,H,A = 9.81,4000,6.37e6
      write(*,*) ' '
      write(*,'(A,F8.1,A)') '               TIME-STEP: ', dt, ' s'
      write(*,*) ' '
      write(*,*) ' '

      !!! output file parameters
      howoften = 12*fac3		! write out history (instantaneously) every "howoften" timesteps
      wrtflg=0				! integer flag for output

      pwr2 = 1
      do while (2**pwr2.lt.N)
	 pwr2 = pwr2+1
      enddo
      
      !!!!! define gaussian zeros/weights !!!!!
      call gauss_legendre(gauzer,gauwts,M)
      do j=1,M/2
         y(j) = -dasin(gauzer(j))
         mu(j) = -gauzer(j)
	 lats(j) = 180.d0*y(j)/pie
      enddo
      do j=M/2+1,M
         y(j) =  dasin(gauzer(j))
         mu(j) =  gauzer(j)
	 lats(j) = 180.d0*y(j)/pie
      enddo
      do i=1,N
	 x(i)=dble(i-1)*2.d0*pie/dble(N)
	 lons(i) = 180.d0*x(i)/pie
      enddo

      t = 0.d0                  ! model time in seconds
      tt = 1                    ! model iteration number
      recno = 1                 ! number of netcdf record
      randseed = 35268
      call srand(randseed)

      !!!! initialize physical/spectral variables
      call initialize(N,  M,  trunc,  rhom,  pwr2, x,  y,
     +                FILE_IN,  data_init,
     +                A, H0,  g,  f, gauzer, gauwts,
     +                hb, h, etab, eta, u, v, d, z,
     +                psi, chi, phi, phib,
     +                etabsp, etasp, usp, vsp, dsp, zsp)

C      call nrgprod( N, M, trunc, rhom, pwr2,
C     +              x, y, gauzer, gauwts,
C     +              g, H0, 1.d0, v,h, 436)

      !!! write out initial state
      call ncinit(FILE_OUT,N,M,lats,lons,A,f,g,H0,hb)
      call ncwritehist(h,u,v,z,d,t,N,M,FILE_OUT,recno)

C      stop

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! first time-step: forward
      dt = dt*0.5d0
      call timestepper(N, M, trunc, rhom, tt, pwr2,
     +                 dt, A, g, H0, f, t, kf, nue, nunm, rob,
     +                 x, y, gauzer, gauwts,
     +                 hb, h, etab, eta, u, v, d, z,
     +                 etabsp, etasp, usp, vsp, esp, dsp, zsp,
     +                 phi, psi, chi, donrg, nrg, randseed)
      dt = dt*2.d0

      !!! cycle spectral variables !!!
      do i=1,N
         do j=1,M
	    psi(i,j,1)=psi(i,j,2)
	    psi(i,j,2)=psi(i,j,3)
            psi(i,j,3)=dcmplx(0.d0,0.d0)
	    phi(i,j,1)=phi(i,j,2)
	    phi(i,j,2)=phi(i,j,3)
            phi(i,j,3)=dcmplx(0.d0,0.d0)
	    chi(i,j,1)=chi(i,j,2)
	    chi(i,j,2)=chi(i,j,3)
            chi(i,j,3)=dcmplx(0.d0,0.d0)
         enddo
      enddo

      !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      !!! remaining time-steps: leap-frog
      do tt=2,Tf*fac3/2 + 1
	 randseed = idnint(100000.d0*rand())

      call timestepper(N, M, trunc, rhom, tt, pwr2,
     +                 dt, A, g, H0, f, t, kf, nue, nunm, rob,
     +                 x, y, gauzer, gauwts,
     +                 hb, h, etab, eta, u, v, d, z,
     +                 etabsp, etasp, usp, vsp, esp, dsp, zsp,
     +                 phi, psi, chi, donrg, nrg, randseed)

      !!! write history !!!
      t = dt*dble(tt-1)
      wrtflg = wrtflg+1
      if (wrtflg.eq.howoften) then
      	wrtflg=0
      	recno=recno+1
      	call ncwritehist(h,u,v,z,d,t,N,M,FILE_OUT,recno) 
      endif

      !!! cycle spectral variables !!!
      do i=1,N
         do j=1,M
	    psi(i,j,1)=psi(i,j,2)
	    psi(i,j,2)=psi(i,j,3)
            psi(i,j,3)=dcmplx(0.d0,0.d0)
	    phi(i,j,1)=phi(i,j,2)
	    phi(i,j,2)=phi(i,j,3)
            phi(i,j,3)=dcmplx(0.d0,0.d0)
	    chi(i,j,1)=chi(i,j,2)
	    chi(i,j,2)=chi(i,j,3)
            chi(i,j,3)=dcmplx(0.d0,0.d0)
         enddo
      enddo

      enddo !!!!!!!!!! time loop

      return
      end

