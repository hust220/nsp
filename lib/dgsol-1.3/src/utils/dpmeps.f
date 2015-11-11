      double precision function dpmeps()
c     **********
c
c     Subroutine dpmeps
c
c     This subroutine computes the machine precision parameter
c     dpmeps as the smallest floating point number such that
c     1 + dpmeps differs from 1.
c
c     This subroutine is based on the subroutine machar described in
c
c     W. J. Cody, MACHAR: A subroutine to dynamically determine
c     machine parameters, ACM Trans. Math. Software 14 (1988), 303-311.
c
c     The function statement is
c
c       double precision dpmeps()
c
c     MINPACK-2 Project. February 1991.
c     Argonne National Laboratory and University of Minnesota.
c     Brett M. Averick.
c
c     **********
      integer i, ibeta, irnd, it, itemp, negep
      double precision a, b, beta, betain, betah, temp, tempa, temp1
      double precision zero, one, two
      data zero, one, two/0.0d0, 1.0d0, 2.0d0/

c     Determine ibeta, beta ala malcolm.

      a = one
      b = one
   10 continue
      a = a + a
      temp = a + one
      temp1 = temp - a
      if (temp1-one .eq. zero) go to 10
   20 continue
      b = b + b
      temp = a + b
      itemp = int(temp-a)
      if (itemp .eq. 0) go to 20
      ibeta = itemp
      beta = dble(ibeta)

c     Determine it, irnd.

      it = 0
      b = one
   30 continue
      it = it + 1
      b = b*beta
      temp = b + one
      temp1 = temp - b
      if (temp1-one .eq. zero) go to 30
      irnd = 0
      betah = beta/two
      temp = a + betah
      if (temp-a .ne. zero) irnd = 1
      tempa = a + beta
      temp = tempa + betah
      if ((irnd .eq. 0) .and. (temp-tempa .ne. zero)) irnd = 2

c     Determine dpmeps.

      negep = it + 3
      betain = one/beta
      a = one
      do 40 i = 1, negep
         a = a*betain
   40 continue
   50 continue
      temp = one + a
      if (temp-one .ne. zero) go to 60
      a = a*beta
      go to 50
   60 continue
      dpmeps = a
      if ((ibeta .eq. 2) .or. (irnd .eq. 0)) go to 70
      a = (a*(one+a))/two
      temp = one + a
      if (temp-one .ne. zero) dpmeps = a

   70 continue

      end
