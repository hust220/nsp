      subroutine dgsel(n,k,x,ind)
      integer n, k
      integer ind(n)
      double precision x(n)
c     **********
c
c     Subroutine dgsel
c
c     Given an array x of length n, this subroutine permutes
c     the elements of the array so that the kth element is the
c     kth smallest element of the array x.
c
c     The subroutine statement is
c
c       subroutine dgsel(k,n,x)
c
c     where
c
c       k is a positive integer input variable not less than n.
c
c       n is a positive integer input variable.
c
c       x is an array of length n. On output the array x is
c         modified so that x(k) is the kth smallest element of x.
c
c     At each stage of the algorithm there are integers l and u
c     with l .le. u such that if j is an integer in l ,..., u then
c
c           x(i) .le. x(j) if i is less than l,
c
c     and
c
c           x(j) .lt. x(i) if i is greater than u.
c
c     The integers l and u are updated by choosing an element of
c     the subarray x(l) ,..., x(u), and partitioning the subarray
c     so that x(m) is the chosen element, x(i) .lt. x(m) for i in
c     l ,..., m, and x(m) .le. x(i) for i in m+1 ,..., u.
c     The integers l and u are updated so that the kth smallest
c     element of x remains in the subarray.
c
c     **********
      integer i, l,lc,lp,m,p,p1,p2,p3,u
      integer swap

      u = n
      l = 1
      lc = n
      lp = 2*n

c     Start of iteration loop.

   10 continue

c        Choose the partition as the median of the elements in
c        positions l+s*(u-l) for s = 0, 0.25, 0.5, 0.75, 1.
c        Move the partition element into position l.

         p1 = (u+3*l)/4
         p2 = (u+l)/2
         p3 = (3*u+l)/4

c        Order the elements in positions l and p1.

         if (dabs(x(ind(l))) .gt. dabs(x(ind(p1)))) then
            swap = ind(l)
            ind(l) = ind(p1)
            ind(p1) = swap
            end if

c        Order the elements in positions p2 and p3.

         if (dabs(x(ind(p2))) .gt. dabs(x(ind(p3)))) then
            swap = ind(p2)
            ind(p2) = ind(p3)
            ind(p3) = swap
            end if

c        Swap the larger of the elements in positions p1
c        and p3, with the element in position u, and reorder
c        the first two pairs of elements as necessary.

         if (dabs(x(ind(p3))) .gt. dabs(x(ind(p1)))) then
            swap = ind(p3)
            ind(p3) = ind(u)
            ind(u) = swap
            if (dabs(x(ind(p2))) .gt. dabs(x(ind(p3)))) then
               swap = ind(p2)
               ind(p2) = ind(p3)
               ind(p3) = swap
               end if
         else
            swap = ind(p1)
            ind(p1) = ind(u)
            ind(u) = swap
            if (dabs(x(ind(l))) .gt. dabs(x(ind(p1)))) then
               swap = ind(l)
               ind(l) = ind(p1)
               ind(p1) = swap
               end if
            end if

c        Find the third largest element of the four remaining
c        elements (the median), and place in position l.

         if (dabs(x(ind(p1))) .gt. dabs(x(ind(p3)))) then
            if (dabs(x(ind(l))) .le. dabs(x(ind(p3)))) then
               swap = ind(l)
               ind(l) = ind(p3)
               ind(p3) = swap
               end if
         else
            if (dabs(x(ind(p2))) .le. dabs(x(ind(p1)))) then
               swap = ind(l)
               ind(l) = ind(p1)
               ind(p1) = swap
            else
               swap = ind(l)
               ind(l) = ind(p2)
               ind(p2) = swap
               end if
            end if

c        Partition the array about the element in position l.

         m = l
         do 20 i = l+1, u
            if (dabs(x(ind(i))) .lt. dabs(x(ind(l))))then
               m = m + 1
               swap = ind(m)
               ind(m) = ind(i)
               ind(i) = swap
               end if
   20       continue

c        Move the partition element into position m.

         swap = ind(l)
         ind(l) = ind(m)
         ind(m) = swap

c        Adjust the values of l and u.

         if (k .ge. m) l = m+1
         if (k .le. m) u = m-1

c        Check for multiple medians if the length of the subarray
c        has not decreased by 1/3 after two consecutive iterations.

         if (3*(u-l) .gt. 2*lp .and. k .gt. m) then

c           Partition the remaining elements into those elements
c           equal to x(m), and those greater than x(m). Adjust
c           the values of l and u.

            p = m
            do 30 i = m+1, u
               if (dabs(x(ind(i))) .eq. dabs(x(ind(m)))) then
                  p = p + 1
                  swap = ind(p)
                  ind(p) = ind(i)
                  ind(i) = swap
                  end if
   30          continue
            l = p + 1
            if (k .le. p) u = p - 1
            end if

c        Update the length indicators for the subarray.

         lp = lc
         lc = u-l

c        Termination test.

         if (l .lt. u) go to 10
      return

      end
