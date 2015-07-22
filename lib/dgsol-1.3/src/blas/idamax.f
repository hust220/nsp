      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(*), dmax
      integer i, incx, ix, n
c
      idamax = 0
      if (n .lt. 1 .or. incx .le. 0) return
      idamax = 1
      if (n .eq. 1) return
      if (incx .eq. 1) go to 30
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 20 i = 2, n
         if (dabs(dx(ix)) .le. dmax) go to 10
         idamax = i
         dmax = dabs(dx(ix))
   10    continue
         ix = ix + incx
   20 continue

      return
c
c        code for increment equal to 1
c
   30 continue
      dmax = dabs(dx(1))
      do 40 i = 2, n
         if (dabs(dx(i)) .le. dmax) go to 40
         idamax = i
         dmax = dabs(dx(i))
   40 continue

      return

      end
