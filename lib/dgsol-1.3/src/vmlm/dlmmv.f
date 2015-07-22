      subroutine dlmmv(n,m,s,y,rho,scale,mark,v,wa)
      integer n, m, mark
      double precision scale
      double precision s(n,m), y(n,m), rho(m), v(n), wa(m)
c     **********
c
c     This subroutine computes the matrix-vector product H*v
c     where H is the inverse BFGS approximation.
c
c     The matrix H depends on an initial matrix H0, m steps
c     s(1),...,s(m), and m gradient differences y(1),...,y(m).
c     These vectors are stored in the arrays s and y.
c     The most recent step and gradient difference are stored
c     in columns s(1,mark) and y(1,mark), respectively.
c     The initial matrix H0 is assumed to be scale*I.
c
c     The subroutine statement is
c
c       subroutine dlmmv(n,m,s,y,rho,scale,mark,v,wa)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       m is an integer variable.
c         On entry m specifies the number of steps and gradient
c            differences that are stored.
c         On exit m is unchanged.
c
c       s is a double precision array of dimension (n,m)
c         On entry s contains the m steps.
c         On exit s is unchanged.
c
c       y is a double precision array of dimension (n,m)
c         On entry y contains the m gradient differences.
c         On exit y is unchanged.
c
c       rho is a double precision array of dimension m
c         On entry rho contains the m innerproducts (s(i),y(i)).
c         On exit rho is unchanged.
c
c       scale is a double precision variable
c         On entry scale specifies the initial matrix H0 = scale*I.
c         On exit scale is unchanged.
c
c       mark is an integer variable.
c         On entry mark points to the current s(i) and y(i).
c         On exit mark is unchanged.
c
c       v is a double precision array of dimension n.
c         On entry v contains the vector v.
c         On exit v contains the matrix-vector product H*v.
c
c       wa is a double precision work array of dimension m.
c
c     Subprograms called
c
c       Level 1 BLAS ... daxpy, ddot, dscal
c
c     MINPACK-2 project. November 1993.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
c
c     **********
      integer i, k
      double precision beta

      double precision ddot
      external daxpy, ddot, dscal

      k = mark + 1
      do 10 i = 1, m
         k = k - 1
         if (k .eq. 0) k = m
         wa(k) = ddot(n,s(1,k),1,v,1)/rho(k)
         call daxpy(n,-wa(k),y(1,k),1,v,1)
   10 continue
      call dscal(n,scale,v,1)
      do 20 i = 1, m
         beta = wa(k) - ddot(n,y(1,k),1,v,1)/rho(k)
         call daxpy(n,beta,s(1,k),1,v,1)
         k = k + 1
         if (k .eq. m+1) k = 1
   20 continue

      end
