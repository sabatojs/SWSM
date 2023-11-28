CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      subroutine fft(u,N,pwr2,itype)
      implicit none

      integer pwr2,itype,N,i
      double complex u(N)
      double precision x(N),y(N)

      do i=1,N
         x(i)=dreal(u(i))
         y(i)=dimag(u(i))
      enddo
      
      call sffteu( x, y, N, pwr2, itype )

      do i=1,N
         u(i)=dcmplx(x(i),y(i))
      enddo

      return
      end
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
