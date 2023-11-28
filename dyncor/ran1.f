      FUNCTION ran1(idum)
      INTEGER idum,IA,IM,IQ,IR,NTAB,NDIV
      REAL ran1,AM,EPS,RNMX
      PARAMETER (IA=16807,IM=2147483647,AM=1./IM,IQ=127773,IR=2836,
     *		NTAB=32,NDIV=1+(IM-1)/NTAB,EPS=1.2e-7,RNMX=1.-EPS)

CCCC “Minimal” random number generator of Park and Miller with 
CCCC Bays-Durham shuffle and added safeguards. Returns a uniform 
CCCC random deviate between 0.0 and 1.0 (exclusive of the endpoint 
CCCC values). Call with idum a negative integer to initialize; 
CCCC thereafter, do not alter idum between successive deviates in a 
CCCC sequence. RNMX should approximate the largest floating value 
CCCC that is less than 1.

      INTEGER j,k,iv(NTAB),iy
      SAVE iv,iy
      DATA iv /NTAB*0/, iy /0/

      if (idum.le.0.or.iy.eq.0) then 	!!!!! Initialize.
          idum=max(-idum,1) 		!!!!! Be sure to prevent idum = 0.
          do j=NTAB+8,1,-1 		!!!!! Load the shuffle table 
					!!!!! (after 8 warm-ups).
               k=idum/IQ
               idum=IA*(idum-k*IQ)-IR*k
               if (idum.lt.0) idum=idum+IM
               if (j.le.NTAB) iv(j)=idum
          enddo
	  iy=iv(1)
      endif
      k=idum/IQ 			!!!!! Start here when not 
					!!!!! initializing.
      idum=IA*(idum-k*IQ)-IR*k 		!!!!! Compute idum=mod(IA*idum,IM) 
					!!!!! without overflows by
      if (idum.lt.0) idum=idum+IM 	!!!!! Schrage’s method.
      j=1+iy/NDIV 			!!!!! Will be in the range 
					!!!!! 1:NTAB.
      iy=iv(j) 				!!!!! Output previously stored 
					!!!!! value and refill the 
					!!!!! shuffle table.
      iv(j)=idum 
      ran1=min(AM*iy,RNMX) 		!!!!! Because users don’t 
					!!!!! expect endpoint values.
      return
      END
