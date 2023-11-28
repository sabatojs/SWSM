      subroutine nrgprod( N, M, trunc, rhom, pwr2,
     +			  x, y, gauzer, gauwts,
     + 		   	  g, H0, nrg, esp, randseed)

      implicit none

      integer N,M,trunc,rhom

      integer i,j,k,pwr2,randseed

      !! real variables
      double precision e(N,M),e2(N,M)
      double precision y(M),x(N)
      double precision gauzer(M),gauwts(M)
      double precision jjp1
      double precision g,H0,nrg, randnum

      !! spectral variables
      double complex esp(N,M),esp2(N,M)
      double precision amp, phase, a, b

      double precision pie

      pie = 4.d0*datan(1.d0)

      call srand(randseed)
C      print *, rand(), rand(), rand(), rand()
C      print *, rand(86456), rand(), rand(), rand()
C      stop

      do i=1,trunc+1
	 do j=i,i*rhom + trunc+1
	    jjp1 = dble(j*(j-1))
	    amp = exp(-((jjp1-3*(trunc**2)/4)**2)/(2*trunc**2))    
	    randnum =  rand()
	    phase = 2.d0*pie*randnum
	    esp(i,j)=nrg*dcmplx(amp*dcos(phase),amp*dsin(phase))
	 enddo ! meridional index
      enddo ! zonal index

      if (.false.) then
      !!!! inverse spectral transform
      call spectral_transform(esp2,e2,N,M,trunc,rhom,gauzer,gauwts,
     +				-1,0,pwr2)
      !!! modify so there's less stirring in the tropics
      do i=1,N
         do j=1,M
	    e2(i,j) = e2(i,j)*dsin(y(j))**2.d0
         enddo
      enddo
      !!!! spectral transform
      call spectral_transform(esp,e2,N,M,trunc,rhom,gauzer,gauwts,
     +				1,0,pwr2)
C      !!!! inverse spectral transform
C      call spectral_transform(esp,e,N,M,trunc,rhom,gauzer,gauwts,
C     +				-1,0,pwr2)
      endif

      return
      end

