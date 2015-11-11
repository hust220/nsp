      subroutine dgErr
     +   (numAtoms,xOpt,fOpt,numDs,iInd,jInd,lbd,ubd,dErr)
      integer numAtoms,numDs
      integer iInd(*),jInd(*)
      double precision fOpt, dErr(*), xOpt(*), lbd(*), ubd(*)
c
c     **********
c
c     Subroutine dgErr
c
c     This subroutine calculates the distance errors
c     for the given set of distance pairs.
c
c     Argonne National Laboratory.
c     Jorge More' and Zhijun Wu.
c     06/01/1997

c     **********
c

      integer k, i1, i2, i3, j1, j2, j3
      double precision d, r1, r2, r3, lErr, uErr, sqrt

      dErr (1) = 1.0e10
      dErr (2) = 0.0e00
      dErr (3) = 0.0e00

      do k = 1, numDs

         i1 = 3 * (iInd (k) - 1) + 1
         i2 = i1 + 1
         i3 = i2 + 1
         j1 = 3 * (jInd (k) - 1) + 1
         j2 = j1 + 1
         j3 = j2 + 1

         r1 = xOpt (i1) - xOpt (j1)
         r2 = xOpt (i2) - xOpt (j2)
         r3 = xOpt (i3) - xOpt (j3)
         d = r1*r1 + r2*r2 + r3*r3

         uErr = (sqrt (d) - sqrt (ubd (k))) / sqrt (ubd (k))
         if (uErr .lt. 0.0e0) uErr = 0.0e0
         lErr = (sqrt (lbd (k)) - sqrt (d)) / sqrt (lbd (k))
         if (lErr .lt. 0.0e0) lErr = 0.0e0

         if (uErr .lt. lErr) uErr = lErr

         if (dErr (3) .lt. uErr) dErr (3) = uErr
         if (dErr (1) .gt. uErr) dErr (1) = uErr

         dErr (2) = dErr (2) + uErr

      end do 

      dErr (2) = dErr (2) / numDs

      end
