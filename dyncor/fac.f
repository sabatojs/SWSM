CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function fac(n,m)
C
C
      implicit none

      integer n,m
      double precision n1,m1,mpn,mmn,a,b,c

      n1=dble(n)
      m1=dble(m)

      mmn = m1-n1
      mpn = m1+n1

      a=dsqrt((2.d0*m1+1.d0)/2.d0)

C     'b' factorial = (m+n)!/(m-n)!
         b=dsqrt(mpn)
         mpn=mpn-1.d0
         do while (mpn.gt.mmn)
	       b=b*dsqrt(mpn)
         mpn=mpn-1.d0
	 enddo
	 if (n.eq.0.and.m.eq.0) b=1.d0

      fac=a/b

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
