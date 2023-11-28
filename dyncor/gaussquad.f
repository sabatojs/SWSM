      double complex function gaussquad(w,x,a,M1,n,m) 
     $			result (int)
      implicit none

      integer M1,i,n,m
      double precision x(M1),w(M1)
      double complex a(M1)
      double precision nhplgndr,shplgndr,plgndr

      int = dcmplx(0.d0,0.d0)
      do i=1,M1/2
	 nhplgndr = plgndr(n,m,x(i))
	 shplgndr = plgndr(n,m,-x(i))
	 int = int + dcmplx(w(i),0.d0)*(  
     +			   a(M1+1-i)*nhplgndr 
      !			       NH          
     +			  +   a(i)*shplgndr  )
      ! 		       SH
      enddo

      end function gaussquad
