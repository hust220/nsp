      subroutine dvmlm(n,x,f,fgrad,frtol,fatol,fmin,task,m,s,y,rho,
     +                 isave,dsave,wa1,wa2)
      character*(*) task
      integer n, m
      integer isave(5)
      double precision f, frtol, fatol, fmin
      double precision x(n), fgrad(n), s(n,m), y(n,m), rho(m), wa1(n),
     +                 wa2(m), dsave(24)
c     **********
c
c     Subroutine dvmlm
c
c     This subroutine computes a local minimizer of a function
c     of n variables by a limited memory variable metric method.
c     The user must evaluate the function and the gradient.
c
c     This subroutine uses reverse communication.
c     The user must choose an initial approximation x to the
c     minimizer, evaluate the function and the gradient at x,
c     and make the initial call with task set to 'START'.
c     On exit task indicates the required action.
c
c     A typical invocation of dvmlm has the following outline:
c
c     Choose a starting vector x.
c     Evaluate the function at x; store in f.
c     Evaluate the gradient at x; store in fgrad.
c
c     task = 'START'
c  10 continue
c        call dvmlm(n,x,f,fgrad,frtol,fatol,fmin,task,m,s,y,rho,
c                   isave,dsave,wa1,wa2)
c        if (task .eq. 'FG') then
c           Evaluate the function at x; store in f.
c           Evaluate the gradient at x; store in fgrad.
c           go to 10
c        else if (task .eq. 'NEWX') then
c           The approximation x, function f, and gradient fgrad
c           are available for inspection.
c           go to 10
c        end if
c
c     NOTE: The user must not alter work arrays between calls.
c
c     The subroutine statement is
c
c       subroutine dvmlm(n,x,f,fgrad,frtol,fatol,fmin,task,m,s,y,rho,
c                        isave,dsave,wa1,wa2)
c
c     where
c
c       n is an integer variable.
c         On entry n is the number of variables.
c         On exit n is unchanged.
c
c       x is a double precision array of dimension n.
c         On entry x is an approximation to the solution.
c         On exit x is the current approximation.
c
c       f is a double precision variable.
c         On entry f is the value of the function at x.
c         On final exit f is the value of the function at x.
c
c       fgrad is a double precision array of dimension n.
c         On entry fgrad is the value of the gradient at x.
c         On final exit fgrad is the value of the gradient at x.
c
c       frtol is a double precision variable.
c         On entry frtol specifies the relative error desired in the
c            function. Convergence occurs if the estimate of the
c            relative error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than frtol.
c         On exit frtol is unchanged.
c
c       fatol is a double precision variable.
c         On entry fatol specifies the absolute error desired in the
c            function. Convergence occurs if the estimate of the
c            absolute error between f(x) and f(xsol), where xsol
c            is a local minimizer, is less than fatol.
c         On exit fatol is unchanged.
c
c       fmin is a double precision variable.
c         On entry fmin specifies a lower bound for the function.
c            The subroutine exits with a warning if f < fmin.
c         On exit fmin is unchanged.
c
c       task is a character variable of length at least 60.
c         On initial entry task must be set to 'START'.
c         On exit task indicates the required action:
c
c            If task(1:2) = 'FG' then evaluate the function and
c            gradient at x and call dvmlm again.
c
c            If task(1:4) = 'NEWX' then a new iterate has been
c            computed. The approximation x, function f, and
c            gradient fgrad are available for examination.
c
c            If task(1:4) = 'CONV' then the search is successful.
c
c            If task(1:4) = 'WARN' then the subroutine is not able
c            to satisfy the convergence conditions. The exit value
c            of x contains the best approximation found.
c
c            If task(1:5) = 'ERROR' then there is an error in the
c            input arguments.
c
c         On exit with convergence, a warning or an error, the
c            variable task contains additional information.
c
c       m is an integer variable.
c          On entry m specifies the amount of storage.
c          On exit m is unchanged.
c
c       s is a double precision work array of dimension (n,m).
c
c       y is a double precision work array of dimension (n,m).
c
c       rho is a double precision work array of dimension m.
c
c       isave is an integer work array of dimension 5.
c
c       dsave is a double precision work array of dimension 24.
c
c       wa1 is a double precision work array of dimension n.
c
c       wa2 is a double precision work array of dimension m.
c
c     Subprograms called
c
c       MINPACK-2 ... dcsrch, dlmmv
c
c       Level 1 BLAS ... daxpy, dcopy, ddot, dnrm2, dscal
c
c     MINPACK-2 Project. April 1995.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick, Richard G. Carter, and Jorge J. More'.
c
c     **********
      double precision zero, one
      parameter (zero=0.0d0,one=1.0d0)

      character*30 work
      integer i, iter, mark
      double precision gd, gd0, f0, scale, sftol, sgtol, stp, stpmin,
     +                 stpmax, sxtol

      double precision dnrm2, ddot
      external dcsrch, dlmmv, daxpy, dcopy, ddot, dnrm2, dscal

      if (task(1:5) .eq. 'START') then

c        Check the input arguments for errors.

         if (n .le. 0) task = 'ERROR: N .LE. 0'
         if (m .le. 0) task = 'ERROR: M .LE. 0'
         if (frtol .le. zero) task = 'ERROR: FRTOL .LE. 0'
         if (fatol .le. zero) task = 'ERROR: FATOL .LE. 0'
         if (f .le. fmin) task = 'ERROR: INITIAL F .LE. FMIN'

c        Exit if there are errors on input.

         if (task(1:5) .eq. 'ERROR') return

c        Initialize local variables.

         iter = 1
         mark = 1

c        Initialize step information.

         scale = dnrm2(n,fgrad,1)
         do 10 i = 1, n
            s(i,1) = fgrad(i)/scale
   10    continue

c        Initialize line search parameters.

         sftol = 1.0d-3
         sgtol = 0.9d0
         sxtol = 0.1d0

c        Set work to start the search.

         work = 'START SEARCH'

      else

c        Restore local variables.

         if (isave(1) .eq. 1) work = 'SEARCH'
         if (isave(1) .eq. 2) work = 'SEARCH DIRECTION'
         iter = isave(2)
         mark = isave(3)
         sftol = dsave(1)
         sgtol = dsave(2)
         sxtol = dsave(3)
         f0 = dsave(4)
         gd = dsave(5)
         gd0 = dsave(6)
         stp = dsave(7)
         stpmin = dsave(8)
         stpmax = dsave(9)
         scale = dsave(10)
      end if

   20 continue

      if (work .eq. 'START SEARCH') then

c        Initialize the line search subroutine.

         f0 = f
         stp = one
         gd0 = -ddot(n,fgrad,1,s(1,mark),1)
         stpmin = zero
         stpmax = (fmin-f0)/(sgtol*gd0)
         stp = min(stp,stpmax)
         call dcopy(n,x,1,wa1,1)
         call dcopy(n,fgrad,1,y(1,mark),1)
         task = 'START SEARCH'
         work = 'SEARCH'

      end if

      if (work .eq. 'SEARCH') then

c        Determine the line search parameter.

         if (f .lt. fmin) then
            task = 'WARNING: F .LT. FMIN'
            go to 30
         end if
         gd = -ddot(n,fgrad,1,s(1,mark),1)

         call dcsrch(stp,f,gd,sftol,sgtol,sxtol,task,stpmin,stpmax,
     +               isave(4),dsave(11))

c        Compute the new iterate.

         call dcopy(n,wa1,1,x,1)
         call daxpy(n,-stp,s(1,mark),1,x,1)

c        Continue if the line search has converged.

         if (task(1:4) .ne. 'CONV' .and.
     +       task .ne. 'WARNING: XTOL TEST SATISFIED') go to 30

c        Compute the step and gradient change.

         iter = iter + 1
         call daxpy(n,-one,fgrad,1,y(1,mark),1)
         call dscal(n,stp,s(1,mark),1)
         rho(mark) = ddot(n,y(1,mark),1,s(1,mark),1)

c        Compute the scale.

         if (rho(mark) .gt. zero) then
            scale = rho(mark)/ddot(n,y(1,mark),1,y(1,mark),1)
         else
            scale = one
         end if

c        Set task to signal a new iterate.
c        Set work to compute a new search direction.

         task = 'NEWX'
         work = 'SEARCH DIRECTION'

c        Test for convergence.

         if (abs(f-f0) .le. fatol .and. stp*abs(gd0) .le. fatol)
     +       task = 'CONVERGENCE: FATOL TEST SATISFIED'
         if (abs(f-f0) .le. frtol*abs(f0) .and.
     +       stp*abs(gd0) .le. frtol*abs(f0))
     +       task = 'CONVERGENCE: FRTOL TEST SATISFIED'

         go to 30

      end if

      if (work .eq. 'SEARCH DIRECTION') then

c        Compute -H*g.

         call dcopy(n,fgrad,1,wa1,1)

         call dlmmv(n,min(m,iter-1),s,y,rho,scale,mark,wa1,wa2)

         mark = mark + 1
         if (mark .eq. m+1) mark = 1
         call dcopy(n,wa1,1,s(1,mark),1)

c        Set task and work to initialize the line search.

         task = 'START SEARCH'
         work = 'START SEARCH'

         go to 20

      end if

   30 continue

c     Save local variables.

      if (work .eq. 'SEARCH') isave(1) = 1
      if (work .eq. 'SEARCH DIRECTION') isave(1) = 2
      isave(2) = iter
      isave(3) = mark
      dsave(1) = sftol
      dsave(2) = sgtol
      dsave(3) = sxtol
      dsave(4) = f0
      dsave(5) = gd
      dsave(6) = gd0
      dsave(7) = stp
      dsave(8) = stpmin
      dsave(9) = stpmax
      dsave(10) = scale

      end
