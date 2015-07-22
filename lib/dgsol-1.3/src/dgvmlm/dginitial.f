      subroutine dgInitial(numAtoms,x,numDs,iInd,jInd,lb,ub,
     +                     seed,intWork)
      integer numAtoms, numDs
      integer seed
      integer iInd(*), jInd(*), intWork(*)
      double precision x(*), lb(*), ub(*)
c     **********
c
c     Subroutine dgInitial
c
c     This subroutine provides an initial point by choosing the 
c     positions of the atoms such that
c
c        || x_i - x_j || = (l_{i,j} + u_{i,j})/2
c
c     A random atom is placed at the origin and other atoms are
c     slected by generating a spanning tree for the graph.
c
c     The algorithm fails if the graph of distance data is not connected.
c
c     Jorge More' and Zhijun Wu.
c     Argonne National Laboratory.
c
c     **********

      integer i, j, k, i1, i2, i3, j1, j2, j3
      integer listStart, listEnd
      double precision r, s, r1, r2, r3

      real surn01

      do i = 1, numAtoms
         intWork(i) = 0
      end do

c     Place a random atom at the origin.

      i = max(int(surn01(seed)*numAtoms),1)

      i1 = 3*(i - 1) + 1
      i2 = i1 + 1
      i3 = i2 + 1

      x(i1) = 0.0d0
      x(i2) = 0.0d0
      x(i3) = 0.0d0

c     Mark the selected atom as chosen.

      intWork(i) = 1
      intWork(numAtoms+1) = i

c     Initialize the markers to the list of atoms to be examined.

      listEnd       = 1
      listStart     = 1

      do while (listStart .le. listEnd)

c        Select an atom i from the list of atoms to be examined.

         i = intWork(numAtoms+listStart)
         listStart = listStart + 1

         i1 = 3*(i - 1) + 1
         i2 = i1 + 1
         i3 = i2 + 1

         do k = 1, numDs

c           Look for an atom j adjacent to the selected atom.

            j = 0
            if (i .eq. iInd(k)) j = jInd(k)
            if (i .eq. jInd(k)) j = iInd(k)

c           Check that there is an atom j adjacent to i.

            if (j .ne. 0) then

c              If atom j has not been selected before, then
c              assign coordinates that satisfy the bounds.

               if (intWork(j) .eq. 0) then

                  j1 = 3*(j - 1) + 1
                  j2 = j1 + 1
                  j3 = j2 + 1

                  r1 = 2.0d0*(surn01(seed) - 0.5d0)
                  r2 = 2.0d0*(surn01(seed) - 0.5d0)
                  r3 = 2.0d0*(surn01(seed) - 0.5d0)

                  r = sqrt(r1*r1 + r2*r2 + r3*r3)
                  s = (sqrt(lb(k)) + sqrt(ub(k))) / 2.0d0

                  x(j1) = x(i1) + s*r1/r
                  x(j2) = x(i2) + s*r2/r
                  x(j3) = x(i3) + s*r3/r

c                 Mark the atom j as selected and place atom j
c                 on the list of atoms to be examined.

                  intWork(j) = 1
                  listEnd = listEnd + 1
                  intWork(numAtoms+listEnd) = j

               end if

            end if

         end do

      end do

      end
