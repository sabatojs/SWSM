      double precision function plgndr(n,m,x)
      !!!!!!  n is the zonal wavenumber
      !!!!!!  m is meridional wavenumber

      implicit none

      integer n,m
      double precision x
      double precision fac
    
      integer i,ll
      double precision fact,pll,pmm,pmmp1,somx2
      
      pmm=1.d0
      if(n.gt.0)then
	somx2=dsqrt((1.d0-x)*(1.d0+x))
	fact=1.d0
	do i=1,n
	   pmm=-pmm*fact*somx2
	   fact=fact+2.d0
	enddo
      endif
      if(m.eq.n)then
	plgndr=pmm
      else
	pmmp1=x*(2.d0*dble(n)+1.d0)*pmm
	if(m.eq.n+1) then
	  plgndr=pmmp1
	else
	  do ll=n+2,m
	    pll=(x*(2.d0*dble(ll)-1.d0)
     +		*pmmp1-(dble(ll)+dble(n)-1.d0)*pmm)/
     +		(dble(ll)-dble(n))
	    pmm=pmmp1
	    pmmp1=pll
	  enddo
	  plgndr=pll
	endif
      endif
          plgndr = plgndr*fac(n,m)
      if (n.eq.0.and.m.ne.0) plgndr=plgndr*dsqrt(dble(m))
      return
      end





