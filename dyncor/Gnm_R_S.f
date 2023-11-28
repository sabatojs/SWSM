      subroutine Gnm_R_S(Gnm,R,S,N,M,A,trunc,rhom,gauzer,gauwts,pwr2)

      implicit none

      double complex gaussquad
      double precision epsln

      integer N,M,trunc,rhom,pwr2
      integer i,j,k

      double precision gauwts(M),gauzer(M),ommu2
      double complex Gnm(N,M),Rn(N,M),Sn(N,M)
      double complex r1(N),s1(N),r2(M),s2(M),s3(M)
      double precision R(N,M),S(N,M),A

      !!!!!!!! FFT !!!!!!!!
      do j=1,M
         do i=1,N
            r1(i)=dcmplx(R(i,j),0.d0)
            s1(i)=dcmplx(S(i,j),0.d0)
	 enddo
      call fft(r1,N,pwr2,1)
      call fft(s1,N,pwr2,1)
         do i=1,trunc+1
            Rn(i,j)=r1(i)
            Sn(i,j)=s1(i)
	 enddo
      enddo

      !!!!!!!!! Gaussian Quadrature !!!!!!!!!
      do i=1,trunc+1
	do j=i,i*rhom + trunc
	   epsln = dsqrt(dble(j*j-(i-1)*(i-1)) / 
     +			(4.d0*dble(j*j)-1.d0))
           do k=1,M/2
	   ommu2 = (1.d0 - gauzer(k)*gauzer(k))
           r2(k) = dcmplx(0.d0,1.d0)*dcmplx(dble(i-1)/ommu2,0.d0)
     +		   *Rn(i,k)
           s2(k) = Sn(i,k)*dcmplx(dble(j)*gauzer(k)/ommu2,0.d0)
           s3(k) = dcmplx((2.d0*dble(j-1)+1.d0)*epsln/ommu2
     +			,0.d0)*Sn(i,k)
	   enddo
           do k=M/2+1,M
	   ommu2 = (1.d0 - gauzer(k)*gauzer(k))
           r2(k) = dcmplx(0.d0,1.d0)*dcmplx(dble(i-1)/ommu2,0.d0)
     +		   *Rn(i,k)
           s2(k) = Sn(i,k)*dcmplx(-dble(j)*gauzer(k)/ommu2,0.d0)
           s3(k) = dcmplx((2.d0*dble(j-1)+1.d0)*epsln/ommu2
     +			,0.d0)*Sn(i,k)
	   enddo
        Gnm(i,j)=(
     +		 gaussquad(gauwts,gauzer,s2,M,i-1,j-1)
     +		+gaussquad(gauwts,gauzer,s3,M,i-1,j)
     +		+gaussquad(gauwts,gauzer,r2,M,i-1,j-1)
     +		 )/A
           if (i.ne.1) Gnm(N+2-i,j)=dconjg(Gnm(i,j))
        enddo
      enddo

      return
      end

