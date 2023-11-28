CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      double precision function eps(n1,m1)

      implicit none
                                                                                                                                                 
      integer n1,m1
      double precision n,m
                                                                                                                                                 
      n=dble(n1)
      m=dble(m1)
                                                                                                                                                 
      eps = dsqrt( (m*m - n*n) / (4.d0*m*m - 1.d0) )
                                                                                                                                                 
      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

