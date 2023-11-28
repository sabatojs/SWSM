      subroutine gauss_legendre(x,w,n)

      implicit none
      
      integer n
      double precision x1,x2,pie
      double precision x(n),w(n)
      double precision eps
      parameter (eps=1.5d-16)
   
      integer i,j,m
      double precision p1,p2,p3,pp,xl,xm,z,z1

      pie=4.d0*datan(1.d0)

      x1 = -1.d0
      x2 =  1.d0

      m=(n+1)/2
      xm=0.5d0*(x2+x1)
      xl=0.5d0*(x2-x1)
      do i=1,m
         z=dcos(pie*(dble(i)-.25d0)/(dble(n)+.5d0))

  1      continue
	    p1=1.d0
	    p2=0.d0
	    do j=1,n
	          p3=p2
		  p2=p1
		  p1=((2.d0*dble(j)-1.d0)*z*p2
     +		    -(dble(j)-1.d0)*p3)/dble(j)
            enddo

	    pp=dble(n)*(z*p1-p2)/(z*z-1.d0)
	    z1=z
	    z=z1-p1/pp
	 if(abs(z-z1).gt.eps)goto 1
	 x(i)=xm-x1*z
	 x(n+1-i)=xm+xl*z
	 w(i)=2.d0*xl/((1.d0-z*z)*pp*pp)
	 w(n+1-i)=w(i)
      enddo

      return
      end 


