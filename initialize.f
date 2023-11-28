      subroutine initialize(N,	M,  trunc,  rhom,  pwr2, x,  y, 
     +			    FILE_IN,  data_init,  
     + 			    A, H0,  g,  f, gauzer, gauwts, 
     + 			    hb, h, etab, eta, u, v, d, z, 
     + 			    psi, chi, phi, phib, 
     +                      etabsp, etasp, usp, vsp, dsp, zsp)

      implicit none

      integer N,M,trunc,rhom
      integer i,j,k,tt,pwr2,recno


      !!  real variables
      double precision u(N,M),v(N,M),h(N,M),z(N,M),d(N,M),hb(N,M)
      double precision y(M),x(N),lats(M),lons(N)
      double precision gauzer(M),gauwts(M)
      double precision eps, jjp1
      double precision v0(N,M),u0(N,M),eta(N,M),etab(N,M)

      !!  spectral variables
      double complex esp(N,M),usp(N,M),vsp(N,M)
      double complex phi(N,M,3),psi(N,M,3),chi(N,M,3)
      double complex phib(N,M)
      double complex gba(N,M),gab(N,M),gcd(N,M)
      double complex zsp(N,M),etasp(N,M),etabsp(N,M),dsp(N,M)
      double complex A1(N,M),B1(N,M),C1(N,M)

      double precision pie,A,g,H0,f,dt,t,UU
      double complex im

      logical data_init

      character(20) FILE_IN

      im = dcmplx(0.d0,1.d0)
      pie = 4.d0*datan(1.d0)

      t = 0.d0			! model time in seconds
      tt = 1			! model iteration number
      recno = 1			! number of netcdf record

      !!! initialization of data
      if (data_init) then !***************** read in initial data ******************! 
      call ncreadhist(h,u,v,z,d,t,tt,N,M,FILE_IN,recno)
      do i=1,N
	 do j=1,M
	    hb(i,j)= 0.d0
            eta(i,j) = g*(h(i,j) - H0)
            etab(i,j) = g*hb(i,j)
            u0(i,j)=u(i,j)*dcos(y(j))
            v0(i,j)=v(i,j)*dcos(y(j))
	    do k=1,3
	       psi(i,j,k)=dcmplx(0.d0,0.d0)
	       phi(i,j,k)=dcmplx(0.d0,0.d0)
	       chi(i,j,k)=dcmplx(0.d0,0.d0)
	    enddo
	 enddo
      enddo
      else  !***************** or define i.c.'s here *****************! 
	    !!!!! x,y are lon,lat in radians
	    !!!!!
	    !!!!! enter your own formulae for input variables
	    !!!!!
	    !!!!! default: 45^o latitude mountain with
	    !!!!! westerly flow: stationary mountain lee Rossby waves
      do i=1,N
	 do j=1,M
	    UU = 10.d0
	    ! mountain at 45 latitude and 180 longitude
	    hb(i,j)= 0.d0
     +		+(H0/5.d0)*dexp(-A*A*((x(i)-pie)*(x(i)-pie)*dcos(y(j))
     +		+(y(j)-1.d0*pie/4.d0)*(y(j)-1.d0*pie/4.d0) )
     +		/(2.d0*(1.d6*1.d6)))

	    ! westerly midlatitude jet stream
	    ! flow over the mountain
	    u(i,j) = 0.d0
C	    u(i,j)=UU*( dsin(2.d0*y(j))*
C     +			dsin(2.d0*y(j)) )
	    !if (i.eq.1) then
	    !	print *,u(i,j)
            !endif

	    v(i,j)=0.d0

	    ! approximately geostrophically balanced thickness
	    ! average depth of H0
	    h(i,j)=H0 !-(f*UU*A/(3.d0*g))*(dsin(y(j))**3.d0)

	    ! vorticity and divergence of westerly flow
	    z(i,j)=4.d0*UU*dsin(2.d0*y(j))*dcos(2.d0*y(j))
C(UU/A)*(dsin(2.d0*y(j))/dcos(y(j)))*
C     +			(4.d0*cos(y(j)) - dsin(y(j)))
	    d(i,j)=0.d0

	    ! do not change these
            eta(i,j) = g*(h(i,j) - H0)
            etab(i,j) = g*hb(i,j)

            u0(i,j)=u(i,j)*dcos(y(j))
            v0(i,j)=v(i,j)*dcos(y(j))

	    do k=1,3
	       psi(i,j,k)=dcmplx(0.d0,0.d0)
	       phi(i,j,k)=dcmplx(0.d0,0.d0)
	       chi(i,j,k)=dcmplx(0.d0,0.d0)
	    enddo
	    ! do not change these

	 enddo
      enddo
      endif ! data_init

      !!!  calculate spectral coefficients for zeroth time step
      call spectral_transform(etasp,eta,N,M,trunc,rhom,gauzer,gauwts,
     + 				1,0,pwr2)
      call spectral_transform(etabsp,etab,N,M,trunc,rhom,gauzer,gauwts,
     +				1,0,pwr2)
      call spectral_transform(dsp,d,N,M,trunc,rhom,gauzer,gauwts,
     +				1,0,pwr2)
      call spectral_transform(zsp,z,N,M,trunc,rhom,gauzer,gauwts,
     +				1,0,pwr2)
      do i=1,trunc+1
         do j=i,i*rhom + trunc+1
            if (j.eq.1) then
	    phi(i,j,1) = etasp(i,j)
	    phib(i,j) = etabsp(i,j)
	    chi(i,j,1) = dcmplx(0.d0,0.d0)
	    psi(i,j,1) = dcmplx(0.d0,0.d0)
	    else
            jjp1 = dble(j)*dble(j-1)
	    phi(i,j,1) = etasp(i,j)
	    phib(i,j) = etabsp(i,j)
	    chi(i,j,1) = -dsp(i,j)/jjp1
	    psi(i,j,1) = -zsp(i,j)/jjp1
	    endif
         enddo
      enddo

      !!!! update spectral coefficients for U,V,eta,d,z
      do i=1,trunc+1
	 do j=i,i*rhom + trunc+1
            jjp1=dble(j)*dble(j-1)
	    usp(i,j) = im*dble(i-1)*chi(i,j,1)
     +			+dble(j-2)*eps(i-1,j-1)*psi(i,j-1,1)
     +			-dble(j+1)*eps(i-1,j)*psi(i,j+1,1)
	    vsp(i,j) = im*dble(i-1)*psi(i,j,1)
     +			-dble(j-2)*eps(i-1,j-1)*chi(i,j-1,1)
     +			+dble(j+1)*eps(i-1,j)*chi(i,j+1,1)
	    etasp(i,j) = phi(i,j,1)
	    dsp(i,j) = -jjp1*chi(i,j,1)
	    zsp(i,j) = -jjp1*psi(i,j,1)
	 enddo
	 !!!!!!! carry one more harmonic in j for U,V
	    usp(i,trunc+2) = 
     +			dble(trunc)*eps(i-1,trunc+1)*psi(i,trunc+1,2)
	    vsp(i,trunc+2) =
     +			-dble(trunc)*eps(i-1,trunc+1)*chi(i,trunc+1,2)
      enddo

      !!!! update variables in physical space 
C      call spectral_transform(etasp,eta,N,M,trunc,rhom,gauzer,gauwts,
C     +				-1,0,pwr2)
C      call spectral_transform(dsp,d,N,M,trunc,rhom,gauzer,gauwts,
C     +				-1,0,pwr2)
C      call spectral_transform(zsp,z,N,M,trunc,rhom,gauzer,gauwts,
C     +				-1,0,pwr2)
C      call spectral_transform(vsp,v0,N,M,trunc,rhom,gauzer,gauwts,
C     +				-1,1,pwr2)
C      call spectral_transform(usp,u0,N,M,trunc,rhom,gauzer,gauwts,
C     +				-1,1,pwr2)
      do i=1,N
	 do j=1,M
C	    u(i,j)=u0(i,j)*A/dcos(y(j))
C	    v(i,j)=v0(i,j)*A/dcos(y(j))
C            h(i,j) = eta(i,j)/g + H0
	    ! set spectral components time 1 = time 2
	       phi(i,j,2) = phi(i,j,1)
	       psi(i,j,2) = psi(i,j,1)
	       chi(i,j,2) = chi(i,j,1)
	 enddo
      enddo

      return
      end


