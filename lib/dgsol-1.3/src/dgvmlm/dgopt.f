      subroutine dgopt 
     +(numAtoms,xOpt,fOpt,numDs,iInd,jInd,lbd,ubd,intWkSp,dblWkSp)
      
      integer numAtoms,numDs,iInd(*),jInd(*)
      double precision fOpt,xOpt(*),lbd(*),ubd(*)
      integer intWkSp(*)
      double precision dblWkSp(*)
      
      character*60 task

      integer k,n,m,maxFev,numFev,conSteps,smlDs,intSave(7)
      double precision faTol,fMin,frTol,dblSave(23)
      double precision alpha,beta,lam,lam0,sqrt,sqrt2,sqrt5

c     Parameters for VMLM:

      parameter (m=10)
      parameter (faTol=1.0e-08,frTol=1.0e-08,fMin=-1.0e32) 

      n = 3 * numAtoms

c     Determine the initial lambda value:
c     lam = (1 / sqrt5 - sqrt2) * u * (l - 1) / (u - 1) + sqrt2 * u

      sqrt2 = sqrt (2.0d0)
      sqrt5 = sqrt (5.0d0)

      do k = 1, numDs
         intWkSp (k) = k
      end do

      smlDs = int (0.50d0*numDs)
      call dgsel (numDs,smlDs,ubd,intWkSp)
     
      alpha = 1.0d0 / sqrt5 - sqrt2
      beta  = sqrt (ubd (intWkSp (smlDs)))

      alpha = alpha * beta / (beta - 1.0d0)
      beta  = sqrt2 * beta

      lam0  = sqrt (lbd (intWkSp (smlDs)))
      lam0  = alpha * (lam0 - 1.0d0) + beta

      conSteps = int (1.0d1*lam0)
      lam0 = 0.1d0 * conSteps
      conSteps = 2 * conSteps

c     Start optimization

      do 80 k = 1, conSteps

         lam = (conSteps-k) * lam0 / conSteps

c        Initial function evaluation

         task = 'FG'
         call dgFun
     +   (numAtoms,xOpt,fOpt,dblWkSp(1),sqrt2*lam,task,
     +   numDs,iInd,jInd,lbd,ubd)

c        Initialization

         numFev = 0
         maxFev = 10000

         if (fOpt .lt. faTol) then
            task = 'CONVERGENCE: FATOL SATISFIED AT THE BEGINNING'
            goto 80
         endif

c        Start the vmlm algorithm

         task = 'START'
 
   50    continue

         call dvmlm(n,xOpt,fOpt,dblWkSp(1),frTol,faTol,fMin,task,m,
     +   dblWkSp(n+1),dblWkSp(n*(m+1)+1),dblWkSp(n*(2*m+1)+1),
     +   intSave,dblSave,dblWkSp(n*(2*m+1)+m+1),dblWkSp(n*(2*m+2)+m+1))

c        Do not alter values in dblWkSp.

         if (task .eq. 'FG') then
            call dgFun
     +      (numAtoms,xOpt,fOpt,dblWkSp(1),sqrt2*lam,task,
     +      numDs,iInd,jInd,lbd,ubd)
            numFev = numFev + 1
            if (numFev .lt. maxFev) go to 50
            task = 'STOP: MAXFEV LIMIT REACHED'
         end if
         if (task .eq. 'NEWX') go to 50

   80 continue

      end

