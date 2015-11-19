      subroutine dgFun(numAtoms,x,f,g,lam,task,
     +           numDs,iInd,jInd,lb,ub)
      character*(*) task
      integer numAtoms,numDs
      integer iInd(*),jInd(*)
      double precision f,lam
      double precision x(*),g(*),lb(*),ub(*)
c
c     **********
c
c     Subroutine dgFun
c
c     This subroutine computes the function and the gradient 
c     of the distance geometry problem.

c     Argonne National Laboratory.
c     Jorge More' and Zhijun Wu.
c     03/01/97
c
c     **********
c
      logical fEval, gEval
      integer l, k, i1, i2, i3, j1, j2, j3, n
      parameter (n=5)
      double precision s, s2, r1, r2, r3, eps, sqrt
      parameter (eps=1.0d-4)
      double precision hl, hl2, hu, hu2, hr, dhl, dhu, dhr, hlu
      double precision tl, tl2, tu, tu2, lu, lu2
      double precision u(5), w(5)
c
      u(1)=3.4361591188377376033
      u(2)=2.5327316742327897964
      u(3)=1.7566836492998817734
      u(4)=1.0366108297895136541
      u(5)=0.3429013272237046088
      w(1)=0.76404328552326206291d-5
      w(2)=0.13436457467812326922d-2
      w(3)=0.33874394455481063136d-1
      w(4)=0.24013861108231468641
      w(5)=0.61086263373532579878
c
      fEval = task .eq. 'F' .or. task .eq. 'FG'
      gEval = task .eq. 'G' .or. task .eq. 'FG' 

c
      if (fEval) f = 0.0d0
      if (gEval) then
         do 10 k = 1, 3*numAtoms
            g(k) = 0.0d0
   10    continue
      end if
c
      do 80 k = 1, numDs
	 i1 = 3*(iInd(k)-1) + 1
	 i2 = i1 + 1
	 i3 = i2 + 1
	 j1 = 3*(jInd(k)-1) + 1
	 j2 = j1 + 1
	 j3 = j2 + 1
	 r1 = x(i1) - x(j1)
	 r2 = x(i2) - x(j2)
	 r3 = x(i3) - x(j3)
	 s2 = r1*r1 + r2*r2 + r3*r3
	 s = sqrt(s2)
	 hr = 0.0d0
	 dhr = 0.0d0
	 if (s .gt. eps) then
	    do 20 l = 1, n
	       lu = lam*u(l)
	       tu = s + lu
	       tl = s - lu
	       tu2 = tu*tu
	       tl2 = tl*tl
	       hl = 0.0d0
	       if (tl2 .lt. lb(k)) then
	          hl = (tl2-lb(k))/lb(k)
               endif
	       hl2 = hl*hl
	       hu = 0.0d0
	       if (tu2 .lt. lb(k)) then
	          hu = (tu2-lb(k))/lb(k) 
               endif
	       hu2 = hu*hu
	       if (fEval) then
	          hr = hr + (tl*hl2 + tu*hu2)*w(l)
               endif
	       if (gEval) then
		  dhl = 4.0d0*tl2*hl/lb(k)
		  dhu = 4.0d0*tu2*hu/lb(k)
                  dhr = dhr + s*(dhl + dhu)*w(l)
		  dhr = dhr + lu*(hl2 - hu2)*w(l)
               endif
	       hl = 0.0d0
	       if (tl2 .gt. ub(k)) then
	          hl = (tl2-ub(k))/ub(k)
               endif
	       hl2 = hl*hl
	       hu = 0.0d0
	       if (tu2 .gt. ub(k)) then
	          hu = (tu2-ub(k))/ub(k)
               endif
	       hu2 = hu*hu
	       if (fEval) then
	          hr = hr + (tl*hl2 + tu*hu2)*w(l)
               endif
	       if (gEval) then
		  dhl = 4.0d0*tl2*hl/ub(k)
		  dhu = 4.0d0*tu2*hu/ub(k)
                  dhr = dhr + s*(dhl + dhu)*w(l)
		  dhr = dhr + lu*(hl2 - hu2)*w(l)
               endif
   20       continue
	    hr = hr/s
	    dhr = dhr/s2
         else
	    do 30 l = 1, n
	       lu = lam*u(l)
	       tu = s + lu
	       tl = s - lu
	       lu2 = lu*lu
	       tu2 = tu*tu
	       tl2 = tl*tl
	       hl = 0.0d0
	       hu = 0.0d0
               if (lu2 .lt. lb(k)) then
		  hl = (tl2-lb(k))/lb(k)
		  hl2 = hl*hl
		  hu = (tu2-lb(k))/lb(k)
		  hu2 = hu*hu
		  if (fEval) then
		     hr = hr + (hl2 + hu2)*w(l)
		     hlu = s2 + lu2
                     hr = hr + 8.0d0*lu2*(hlu-lb(k))*w(l)/lb(k)/lb(k)
                  endif
		  if (gEval) then
		     hlu = s2 + 5.0d0*lu2
		     dhr = dhr + 8.0d0*(hlu -lb(k))*w(l)/lb(k)/lb(k)
                  endif
               endif
	       hl = 0.0d0
	       hu = 0.0d0
               if (lu2 .gt. ub(k)) then
		  hl = (tl2-ub(k))/ub(k)
		  hl2 = hl*hl
		  hu = (tu2-ub(k))/ub(k)
		  hu2 = hu*hu
		  if (fEval) then
		     hr = hr + (hl2 + hu2)*w(l)
		     hlu = s2 + lu2
                     hr = hr + 8.0d0*lu2*(hlu-ub(k))*w(l)/ub(k)/ub(k)
                  endif
		  if (gEval) then
		     hlu = s2 + 5.0d0*lu2
		     dhr = dhr + 8.0d0*(hlu-ub(k))*w(l)/ub(k)/ub(k)
                  endif
               endif
   30       continue
         endif
	 if (fEval) then
            f = f + hr
         endif
	 if (gEval) then
	    if (s .gt. eps) then
	       r1 = dhr*r1/s
	       r2 = dhr*r2/s
	       r3 = dhr*r3/s
	    else
	       r1 = dhr*r1
	       r2 = dhr*r2
	       r3 = dhr*r3
            endif
	    g(i1) = g(i1) + r1
	    g(i2) = g(i2) + r2
	    g(i3) = g(i3) + r3
	    g(j1) = g(j1) - r1
	    g(j2) = g(j2) - r2
	    g(j3) = g(j3) - r3
	 end if
   80 continue
c
      end
