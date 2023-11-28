      subroutine spectral_transform(a,u,N,M,trunc,rhom,gauzer,
     +				    gauwts,itype,jtype,pwr2)

      implicit none

      double complex gaussquad
      double precision plgndr

      integer N,M,trunc,pwr2,itype,jtype,rhom
      integer i,j,k

      double precision pie
      double precision gauwts(M),gauzer(M)
      double complex a(N,M),amu(N,M)
      double complex u1(N),a1(M)
      double precision u(N,M)

      !!!!! jtype tells spectral_transform if it is for U or V
      !!!!! in which case, we need to carry one more legendre polynomial
      !!!!! in the expansion trunc ==> trunc + 1 for j loops
      !!!!!   jtype=0 ==> standard    jtype=1 ==> U/V 

      if (itype.eq.1) then
      !!!!!!!!!!!! forward transform physical-to-spectral !!!!!!!!!!
      !!!!!!!!! FFT !!!!!!!!!
      do j=1,M
         do i=1,N
            u1(i)=dcmplx(u(i,j),0.d0)
	 enddo
      call fft(u1,N,pwr2,1)
         do i=1,trunc+1
            amu(i,j)=u1(i)
	 enddo
      enddo
      !!!!!!!!! Gaussian Quadrature !!!!!!!!
      do i=1,trunc+1
	 do j=i,i*rhom + trunc + 1 + jtype*(1-rhom)
            do k=1,M
	       a1(k)=amu(i,k)
	    enddo
            a(i,j)=gaussquad(gauwts,gauzer,a1,M,i-1,j-1)
	    if (i.gt.1) then
               a(N+2-i,j)=dconjg(a(i,j))!!!! the factor of (-1)**n is 
					!!!! wrapped in the plgndr routine
	    endif
	 enddo
      enddo

      else !!! itype
      !!!!!!!!!!!! inverse transform spectral-to-physical !!!!!!!!!!!
      !!!!!!! sum over Legendre polynomials !!!!!!!!
      do i=1,trunc+1
	 do k=1,M/2
	    amu(i,k)=dcmplx(0.d0,0.d0)
	    amu(i,M+1-k)=dcmplx(0.d0,0.d0)
	    do j=i,i*rhom + trunc + 1 + jtype*(1-rhom)
	       amu(i,k)=amu(i,k)+a(i,j)
     +		        *plgndr(i-1,j-1,-gauzer(k))
	       amu(i,M+1-k)=amu(i,M+1-k)+a(i,j)
     +		        *plgndr(i-1,j-1,gauzer(k))
	    enddo
            if (i.gt.1) then
	       amu(N+2-i,k)=dconjg(amu(i,k))
	       amu(N+2-i,M+1-k)=dconjg(amu(i,M+1-k))
	    endif
	 enddo
      enddo

      !!!!!!!!! IFFT !!!!!!!!!
      do j=1,M
         do i=1,N
            u1(i)=amu(i,j)
	 enddo
      call fft(u1,N,pwr2,-1)
         do i=1,N
            u(i,j)=dreal(u1(i))
	 enddo
      enddo
     
      endif !!!!! itype

      return
      end
