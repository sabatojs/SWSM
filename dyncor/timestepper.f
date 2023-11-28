      subroutine timestepper(N, M, trunc, rhom, tt, pwr2, 
     +			     dt, A, g, H0, f, t, kf, nue, nunm, rob, 
     +			     x, y, gauzer, gauwts,
     + 			     hb, h, etab, eta, u, v, d, z,
     +		             etabsp, etasp, usp, vsp, esp, dsp, zsp,
     +			     phi, psi, chi, donrg, nrg, randseed)

      implicit none

      integer N,M,trunc,rhom

      integer i,j,k,tt,pwr2,randseed

      logical donrg

      double precision u(N,M),v(N,M),h(N,M),z(N,M),d(N,M),hb(N,M)
      double precision y(M),x(N)
      double precision gauzer(M),gauwts(M)
      double precision eps

      !! real variables
      double precision e(N,M),vz(N,M),uz(N,M)
      double precision uphi(N,M),vphi(N,M)
      double precision muz(N,M),v0(N,M),u0(N,M),eta(N,M),etab(N,M)
      double precision mu(M),jjp1,jjm1,jp1jp2,fac1,nrg

      !! spectral variables
      double complex esp(N,M),usp(N,M),vsp(N,M)
      double complex phi(N,M,3),psi(N,M,3),chi(N,M,3)
      double complex gba(N,M),gab(N,M),gcd(N,M)
      double complex zsp(N,M),etasp(N,M),etabsp(N,M),dsp(N,M)
      double complex A1(N,M),B1(N,M),C1(N,M)
      double complex enrg(N,M)

      double precision A,g,H0,f,dt,t,rob,nue,nunm,kf
      double complex im,pie

      im = dcmplx(0.d0,1.d0)
      pie = 4.d0*datan(1.d0)

      !!! grid multiplication of nonlinear terms
      do i=1,N
	 do j=1,M
	    e(i,j)=0.5d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            u0(i,j) = u(i,j)*dcos(y(j))/A
            v0(i,j) = v(i,j)*dcos(y(j))/A
            uz(i,j) = u(i,j)*dcos(y(j))*z(i,j)
            vz(i,j) = v(i,j)*dcos(y(j))*z(i,j)
            muz(i,j) = -uz(i,j)
            uphi(i,j) = u(i,j)*dcos(y(j))*(eta(i,j)-etab(i,j))
            vphi(i,j) = v(i,j)*dcos(y(j))*(eta(i,j)-etab(i,j))
	 enddo
      enddo

      !!! spectral transforms
      call spectral_transform(esp,e,N,M,trunc,rhom,gauzer,gauwts,
     +				1,1,pwr2)
      call Gnm_R_S(gba,vz,muz,N,M,A,trunc,rhom,gauzer,gauwts,pwr2)
      call Gnm_R_S(gab,uz,vz,N,M,A,trunc,rhom,gauzer,gauwts,pwr2)
      call Gnm_R_S(gcd,uphi,vphi,N,M,A,trunc,rhom,gauzer,gauwts,pwr2)

      !!! compute energy production !!!
      if (donrg) then
         call nrgprod( N, M, trunc, rhom, pwr2,
     +                 x, y, gauzer, gauwts,
     +                 g, H0, nrg, enrg, randseed)
      else
	 do i=1,N
	    do j=1,M
	       enrg(i,j)=0.d0
	    enddo
         enddo
      endif

      !!!! time-step chi,phi,psi
      do i=1,trunc+1
	 do j=i,i*rhom + trunc+1

         if (j.eq.1) then !!! enforce mass conservation in spectral space

            phi(i,j,3) = phi(i,j,2)
            chi(i,j,3) = dcmplx(0.d0,0.d0)
            psi(i,j,3) = dcmplx(0.d0,0.d0)

         else

            jjp1 = dble(j)*dble(j-1)
            jjm1 = dble(j-2)*dble(j-1)               
            jp1jp2 = dble(j)*dble(j+1)               

	    nunm = nue*(jjp1*jjp1*jjp1 - 4.d0*jjp1)/(A**4.d0)

            fac1 = 1.d0 + 4.d0*dt*dt*jjp1*g*H0/(A*A)

            usp(i,j) = dcmplx(0.d0,1.d0)*dble(i-1)*chi(i,j,2)
     +                  +dble(j-2)*eps(i-1,j-1)*psi(i,j-1,2)
     +                  -dble(j+1)*eps(i-1,j)*psi(i,j+1,2)
            vsp(i,j) = dcmplx(0.d0,1.d0)*dble(i-1)*psi(i,j,2)
     +                  -dble(j-2)*eps(i-1,j-1)*chi(i,j-1,2)
     +                  +dble(j+1)*eps(i-1,j)*chi(i,j+1,2)

	    !!!! divergence equation !!!!
            A1(i,j)=chi(i,j,1) -(2.d0*dt/jjp1)*(
     +               gba(i,j)					! nonlinear terms
     +              +jjp1*esp(i,j)/(A*A)           
     +              +jjp1*enrg(i,j)/(A*A)               ! forcing term     
     +              -f*(usp(i,j)				! coriolis terms
     +                  +jjm1*eps(i-1,j-1)*psi(i,j-1,2)
     +                  +jp1jp2*eps(i-1,j)*psi(i,j+1,2))
     +		    +nunm*chi(i,j,2) 				! sub-grid diffusion 
     +		    +jjp1*kf*chi(i,j,2) 			! large-scale frictional damping
     +                           )
	    !!!! vorticity equation !!!!
            B1(i,j)=(-2.d0*dt/jjp1)*(
     +              -gab(i,j)					! nonlinear terms
     +              -f*(vsp(i,j)				! coriolis terms
     +                  +jjm1*eps(i-1,j-1)*chi(i,j-1,2)
     +                  -jp1jp2*eps(i-1,j)*chi(i,j+1,2))
     +		    +nunm*psi(i,j,2) 				! sub-grid diffusion 
     +		    +jjp1*kf*psi(i,j,2) 			! large-scale frictional damping
     +                           )
	    !!!! thickness equation !!!!
            C1(i,j)=phi(i,j,1)+2.d0*dt*(
     +			-gcd(i,j)				! nonlinear terms
     +			-nue*jjp1*jjp1*phi(i,j,2)/(A**4.d0)	! sub-grid diffusion
     +					)
	    !!! update !!!
            chi(i,j,3) = (A1(i,j) - 2.d0*dt*C1(i,j)/(A*A))/fac1
            psi(i,j,3) = psi(i,j,1) + B1(i,j)
            phi(i,j,3) = C1(i,j)+2.d0*dt*jjp1*g*H0*chi(i,j,3)
	    !!! energy production term !!!
C	    if (donrg) then
C	       phi(i,j,3) = phi(i,j,3) + dt*enrg(i,j)
C	       psi(i,j,3) = psi(i,j,3) !+ dt*enrg(i,j)*(f/g)**2.d0
C	    endif

         endif

      !!!!!! time (Robert–Asselin) filtering !!!!!!
            psi(i,j,2) = psi(i,j,2) + rob*
     +			(psi(i,j,1)-2.d0*psi(i,j,2)+psi(i,j,3))
            phi(i,j,2) = phi(i,j,2) + rob*
     +			(phi(i,j,1)-2.d0*phi(i,j,2)+phi(i,j,3))
            chi(i,j,2) = chi(i,j,2) + rob*
     +			(chi(i,j,1)-2.d0*chi(i,j,2)+chi(i,j,3))
	 enddo ! meridional index
      enddo ! zonal index
      !!!! update spectral coefficients for U,V,eta,d,z
      do i=1,trunc+1
	 do j=i,i*rhom + trunc+1
            jjp1=dble(j)*dble(j-1)
	    usp(i,j) = im*dble(i-1)*chi(i,j,2)
     +			+dble(j-2)*eps(i-1,j-1)*psi(i,j-1,2)
     +			-dble(j+1)*eps(i-1,j)*psi(i,j+1,2)
	    vsp(i,j) = im*dble(i-1)*psi(i,j,2)
     +			-dble(j-2)*eps(i-1,j-1)*chi(i,j-1,2)
     +			+dble(j+1)*eps(i-1,j)*chi(i,j+1,2)
	    etasp(i,j) = phi(i,j,2)
	    dsp(i,j) = -jjp1*chi(i,j,2)
	    zsp(i,j) = -jjp1*psi(i,j,2)
	 enddo
	 !!!!!!! carry one extra harmonic in j for U,V to avoid aliasing
	    usp(i,trunc+2) = 
     +			dble(trunc)*eps(i-1,trunc+1)*psi(i,trunc+1,2)
	    vsp(i,trunc+2) =
     +			-dble(trunc)*eps(i-1,trunc+1)*chi(i,trunc+1,2)
      enddo
      !!!! inverse spectral transforms
      call spectral_transform(etasp,eta,N,M,trunc,rhom,gauzer,gauwts,
     +				-1,0,pwr2)
      call spectral_transform(dsp,d,N,M,trunc,rhom,gauzer,gauwts,
     +				-1,0,pwr2)
      call spectral_transform(zsp,z,N,M,trunc,rhom,gauzer,gauwts,
     +				-1,0,pwr2)
      call spectral_transform(vsp,v0,N,M,trunc,rhom,gauzer,gauwts,
     +				-1,1,pwr2)
      call spectral_transform(usp,u0,N,M,trunc,rhom,gauzer,gauwts,
     +				-1,1,pwr2)

      !!! update variables in physical space 
      do i=1,N
	 do j=1,M
	    u(i,j)=u0(i,j)*A/dcos(y(j))
	    v(i,j)=v0(i,j)*A/dcos(y(j))
	    e(i,j)=0.5d0*(u(i,j)*u(i,j)+v(i,j)*v(i,j))
            h(i,j) = eta(i,j)/g + H0
	 enddo
      enddo

      return
      end

