c----------------------------------------------------------------------- 
c 
c     This file is part of the Test Set for IVP solvers 
c     http://www.dm.uniba.it/~testset/ 
c 
c        generic DDASSL driver 
c 
c     DISCLAIMER: see 
c     http://www.dm.uniba.it/~testset/disclaimer.php 
c 
c     The most recent version of this source file can be found at 
c     http://www.dm.uniba.it/~testset/src/drivers/dassld.f 
c 
c     This is revision 
c     $Id: dassld.f,v 1.6 2006/10/02 10:19:08 testset Exp $ 
c 
c----------------------------------------------------------------------- 
 
      program dassld 
 
      integer maxneq,maxord 
      parameter (maxneq=400,maxord=5) 
      integer lrw, liw 
      parameter (lrw = 40+(maxord+4)*maxneq+maxneq**2, liw = 20+maxneq) 
 
      integer neqn,ndisc,mljac,mujac,mlmas,mumas,ind(maxneq), 
     +        info(15),iwork(liw),ierr,ipar(2),idid 
      double precision y(maxneq),yprime(maxneq),t(0:100), 
     +                 h0,rtol(maxneq),atol(maxneq), 
     +                 rwork(lrw),rpar(1) 
      logical numjac,nummas,consis,tolvec 
 
      double precision yref(maxneq) 
      character fullnm*40, problm*8, type*3 
      character driver*8, solver*8 
      parameter (driver = 'dassld', solver='DDASSL') 
      external ideres,idejac 
      external oderes,odejac 
      external daeres,daejac 
 
      double precision solmax 
      real gettim, timtim, cputim 
      external gettim 
      double precision scd,mescd 
 
      integer ldjac,ldmas 
      double precision mas(maxneq*maxneq) 
      common /neqcom/ neqn 
      common /jaccom/ ldjac,mljac,mujac 
      common /mascom/ mas,ldmas,mlmas,mumas 
 
 
      character  fileout*140,  filepath*100 
      character formatout*30, namefile*100 
      logical printsolout, solref, printout 
      integer nindsol, indsol(maxneq) 
   
 
      integer i, j, icount(11:15), it 
 
      do 10 i=1,11 
         info(i) = 0 
   10 continue 
c 
c     we use DASSL's intermediate-output mode (restart every step) 
c     instead of stopping/restarting after DASSL's IDID=-1 
c     ('500 STEPS TAKEN ON THIS CALL BEFORE REACHING TOUT') 
c 
      info(3) = 1 
c 
c     for an accurate solution in the Test Set endpoint 
c     (the scd-value) we do not want the solution to be interpolated 
c     so not only for discontinuities (as it should) but also 
c     for the Test Set endpoint (see also rwork(1)): 
c 
      info(4) = 1 
 
c----------------------------------------------------------------------- 
c     check the problem definition interface date 
c----------------------------------------------------------------------- 
 
      call chkdat(driver,solver,20060828) 
 
c----------------------------------------------------------------------- 
c     get the problem dependent parameters 
c----------------------------------------------------------------------- 
 
      call prob(fullnm,problm,type, 
     +          neqn,ndisc,t, 
     +          numjac,mljac,mujac, 
     +          nummas,mlmas,mumas, 
     +          ind) 
      if (type.eq.'IDE') then 
         do 20 i=1,neqn 
            if (problm.ne.'wheel') then 
               if (ind(i).gt.1) then 
                  print *, 'DASSLD: ERROR: ', 
     +                     'DASSL can not solve higher index problems' 
                  stop 
               endif 
            endif 
   20    continue 
      elseif (type.eq.'ODE') then 
         nummas = .false. 
         mlmas    = 0 
         mumas    = 0 
      elseif (type.eq.'DAE') then 
         nummas = .false. 
         if (problm.ne.'pump') then 
            do 30 i=1,neqn 
               if (ind(i).gt.1) then 
                  print *, 'DASSLD: ERROR: ', 
     +                     'DASSL can not solve higher index problems' 
                  stop 
               endif 
   30       continue 
         endif 
         if (mlmas.lt.neqn) then 
            ldmas = mlmas + mumas + 1 
         else 
            ldmas = neqn 
         endif 
         call meval(ldmas,neqn,t(0),y,yprime,mas,ierr,rpar,ipar) 
      else 
         print *, 'DASSLD: ERROR: ', 
     +            'unknown Test Set problem type', type 
         stop 
      endif 
 
      if (.not.numjac) info(5) = 1 
      if (mljac.lt.neqn .and. mlmas.lt.neqn) then 
         info(6) = 1 
         iwork(1)=mljac 
         iwork(2)=mujac 
      endif 
 
      if (mljac.lt.neqn) then 
         ldjac = 2*mljac + mujac + 1 
      else 
         ldjac = neqn 
      endif 
 
c----------------------------------------------------------------------- 
c     get the initial values 
c----------------------------------------------------------------------- 
 
      call init(neqn,t(0),y,yprime,consis) 
      if (type.eq.'ODE') then 
         ierr =  0 
         call feval(neqn,t(0),y,yprime,yprime,ierr,rpar,ipar) 
         if (ierr.ne.0) then 
            print *, 'DASSLD: ERROR: ', 
     +               'DASSLD could not evaluate f(t(0))' 
            stop 
         endif 
         consis = .true. 
      endif 
      if (.not.consis) info(11)=1 
 
c----------------------------------------------------------------------- 
c     read the tolerances 
c----------------------------------------------------------------------- 
 
      call getinp(driver,problm,solver,fullnm, 
     +            tolvec,rtol,atol,h0,solmax) 
 
      call settolerances(neqn,rtol,atol,tolvec) 
 
      if (problm .eq. 'pump') then 
c 
c        no error control for the index-2-variable current: 
c 
         rtol(9) = 1000d0 
         atol(9) = 1000d0 
      endif 
      if (problm .eq. 'wheel') then 
c 
c        the index 2 case: only the Lagrange multipliers have to 
c        be excluded from error control. 
c 
         do 40 i=2,15 
            atol(i) = atol(1) 
            rtol(i) = rtol(1) 
 40      continue 
         do 50 i=16,17 
            atol(i) = 1d10 
            rtol(i) = 1d10 
 50      continue 
         tolvec = .true. 
      endif 
 
      if (tolvec) info(2)=1 
 
 
      call  setoutput(neqn,solref,printsolout, 
     +                        nindsol,indsol) 
 
 
      if (printsolout) then 
  
 
          call getarg(1,filepath) 
          call getarg(2,namefile) 
           
          if (lnblnk(namefile) .gt. 0) then 
            write(fileout,'(a,a,a,a)') filepath(1:lnblnk(filepath)), 
     +     namefile(1:lnblnk(namefile)), solver(1:lnblnk(solver)),'.txt' 
                 
             open(UNIT=90,FILE=fileout) 
  
             call mfileout(namefile,solver,filepath,nindsol,indsol) 
   
          else 
            write(fileout,'(a,a,a,a)') filepath(1:lnblnk(filepath)), 
     +      problm(1:lnblnk(problm)), solver(1:lnblnk(solver)),'.txt' 
                 
             open(UNIT=90,FILE=fileout) 
  
             call mfileout(problm,solver,filepath,nindsol,indsol) 
          end if   
      end if 
 
c----------------------------------------------------------------------- 
c     call of the subroutine DASSL 
c----------------------------------------------------------------------- 
 
      do 60 i=11,15 
         icount(i) = 0 
   60 continue 
 
      if (printsolout) then 
c the initial condition is printed in the oputput file 
        write(formatout,'(a,i5,a)') '(e23.15,',nindsol,'(e23.15))' 
        write(90,formatout)  
     +       t(0), (y(indsol(it)),it=1,nindsol) 
        info(3) = 1 
      end if 
       
 
      timtim = gettim() 
      timtim = gettim() - timtim 
 
      cputim = gettim() 
 
      do 90 i=0,ndisc 
         
         rwork(1) = t(i+1) 
 
   70    continue 
 
            if (type.eq.'IDE') then 
               call ddassl(ideres,neqn,t(i),y,yprime,t(i+1), 
     +                     info,rtol,atol,idid, 
     +                     rwork,lrw,iwork,liw,rpar,ipar,idejac) 
            elseif (type.eq.'ODE') then 
               call ddassl(oderes,neqn,t(i),y,yprime,t(i+1), 
     +                     info,rtol,atol,idid, 
     +                     rwork,lrw,iwork,liw,rpar,ipar,odejac) 
            elseif (type.eq.'DAE') then 
               call ddassl(daeres,neqn,t(i),y,yprime,t(i+1), 
     +                     info,rtol,atol,idid, 
     +                     rwork,lrw,iwork,liw,rpar,ipar,daejac) 
            endif 
 
         if (idid.eq.1) then 
            if (printsolout) then 
             write(90,formatout)  
     +         rwork(4), (y(indsol(it)),it=1,nindsol) 
            end if 
            info(1) = 1 
            goto 70 
         elseif (idid.lt.0) then 
            print *, 'DASSLD: ERROR: ', 
     +               'DASSL returned IDID = ', idid 
            stop 
         endif 
 
         if (printsolout) then 
             write(90,formatout)  
     +         rwork(4), (y(indsol(it)),it=1,nindsol) 
         end if 
         info(1) = 0 
 
         do 80 j=11,15 
            icount(j) = icount(j) + iwork(j) 
   80    continue 
 
   90 continue 
 
       
 
      cputim = gettim() - cputim - timtim 
 
      do 100 i=11,15 
         iwork(i) = icount(i) 
  100 continue 
 
c----------------------------------------------------------------------- 
c     print numerical solution in endpoint and integration statistics 
c----------------------------------------------------------------------- 
      printout = .true. 
      if (solref) then  
        call solut(neqn,t(ndisc+1),yref) 
       
        call getscd(mescd,scd,neqn,yref,y,problm,tolvec,atol,rtol, 
     +             printout) 
      else 
        call printsol(neqn,y,problm) 
      end if 
      call report( 
     +   driver,problm,solver, 
     +   rtol,atol,h0,solmax, 
     +   iwork,cputim,scd,mescd 
     +) 
      end 
 
C======================================================================= 
C     `Test Set for IVP Solvers' IDE wrappers for DASSL 
C======================================================================= 
c 
c     since in DASSL the format of the subroutines for the 
c     function G and the derivative PD differ from the format 
c     in the testset, we transform them 
c 
c         G = f                   -> ideres 
c        PD = dG/dy + cj * dG/dy' -> idejac 
c 
c----------------------------------------------------------------------- 
      subroutine ideres(t,y,yprime,delta,ires,rpar,ipar) 
      integer ires,ipar(*) 
      double precision t,y(*),yprime(*),delta(*),rpar(*) 
      integer neqn 
      common /neqcom/ neqn 
      call feval(neqn,t,y,yprime,delta,ires,rpar,ipar) 
      return 
      end 
 
      subroutine idejac(t,y,yprime,pd,cj,rpar,ipar) 
      integer ipar(*) 
      double precision t,y(*),yprime(*),pd(*),cj,rpar(*) 
c TODO: this one, however, we do not have IDEs with analytical Jacobians 
      stop 'DASSLD: ERROR: idejac not yet implemented' 
      return 
      end 
C======================================================================= 
C     `Test Set for IVP Solvers' ODE wrappers for DASSL 
C======================================================================= 
c 
c     since in DASSL the format of the subroutines for the 
c     residual G and the derivative PD differ from the format 
c     in the testset, we transform them 
c 
c         G = f - y'                             -> oderes 
c        PD = dG/dy + cj * dG/dy' = df/dy - cj   -> odejac 
c 
c----------------------------------------------------------------------- 
      subroutine oderes(t,y,yprime,delta,ires,rpar,ipar) 
      integer ires,ipar(*) 
      double precision t,y(*),yprime(*),delta(*),rpar(*) 
      integer neqn 
      common /neqcom/ neqn 
      integer i 
c compute delta = f - y' 
      call feval(neqn,t,y,yprime,delta,ires,rpar,ipar) 
      do 10 i=1,neqn 
         delta(i) = delta(i) - yprime(i) 
   10 continue 
      return 
      end 
 
      subroutine odejac(t,y,yprime,pd,cj,rpar,ipar) 
      integer ipar(*) 
      double precision t,y(*),yprime(*),pd(*),cj,rpar(*) 
      integer ierr 
      integer neqn,ldjac,mljac,mujac 
      common /neqcom/ neqn 
      common /jaccom/ ldjac,mljac,mujac 
      integer j 
c compute pd = df/dy 
      ierr = 0 
      if (mljac.lt.neqn) then 
         call jeval(ldjac,neqn,t,y,yprime,pd(mljac+1),ierr,rpar,ipar) 
      else 
         call jeval(ldjac,neqn,t,y,yprime,pd(1),ierr,rpar,ipar) 
      endif 
      if (ierr.ne.0) then 
         print *, 'DASSLD: ERROR: ', 
     +            'DASSL can not handle JEVAL IERR' 
         stop 
      endif 
c compute pd = df/dy - cj 
      if (mljac.lt.neqn) then 
         do 10 j=1,neqn 
            pd((j-1)*ldjac+mljac+mujac+1) = 
     +      pd((j-1)*ldjac+mljac+mujac+1) - cj 
   10    continue 
      else 
         do 20 j=1,neqn 
            pd((j-1)*ldjac+j) = 
     +      pd((j-1)*ldjac+j) - cj 
   20    continue 
      endif 
      return 
      end 
C======================================================================= 
C     `Test Set for IVP Solvers' DAE wrappers for DASSL 
C======================================================================= 
c 
c     since in DASSL the format of the subroutines for the 
c     function G and the derivative PD differ from the format 
c     in the testset, we transform them 
c 
c         G = f - My'                              -> daeres 
c        PD = dG/dy + cj * dG/dy' = df/dy - cj*M   -> daejac 
c 
c----------------------------------------------------------------------- 
      subroutine daeres(t,y,yprime,delta,ires,rpar,ipar) 
      integer ires,ipar(*) 
      double precision t,y(*),yprime(*),delta(*),rpar(*) 
      integer maxneq,neqn,ldmas,mlmas,mumas 
      parameter (maxneq=400) 
      double precision mas(maxneq*maxneq) 
      common /neqcom/ neqn 
      common /mascom/ mas,ldmas,mlmas,mumas 
      integer i,j 
c compute f in delta 
      call feval(neqn,t,y,yprime,delta,ires,rpar,ipar) 
c compute delta := f - M*y' 
      if (mlmas.lt.neqn) then 
         do 20 j=1,neqn 
            do 10 i=max(1,j-mumas),min(neqn,j+mlmas) 
               delta(i) = 
     +         delta(i) - mas((j-1)*ldmas+mumas+1-j+i)*yprime(j) 
   10       continue 
   20    continue 
      else 
         do 40 j=1,neqn 
            do 30 i=1,neqn 
               delta(i) = 
     +         delta(i) - mas((j-1)*ldmas+i)*yprime(j) 
   30       continue 
   40    continue 
      endif 
      return 
      end 
 
      subroutine daejac(t,y,yprime,pd,cj,rpar,ipar) 
      integer ipar(*) 
      double precision t,y(*),yprime(*),pd(*),cj,rpar(*) 
      integer ierr 
      integer maxneq,neqn,ldjac,mljac,mujac,ldmas,mlmas,mumas 
      parameter (maxneq=400) 
      double precision mas(maxneq*maxneq) 
      common /neqcom/ neqn 
      common /jaccom/ ldjac,mljac,mujac 
      common /mascom/ mas,ldmas,mlmas,mumas 
      integer i,j 
c compute pd = df/dy 
      ierr = 0 
      if (mljac.lt.neqn) then 
         call jeval(ldjac,neqn,t,y,yprime,pd(mljac+1),ierr,rpar,ipar) 
      else 
         call jeval(ldjac,neqn,t,y,yprime,pd(1),ierr,rpar,ipar) 
      endif 
      if (ierr.ne.0) then 
         print *, 'DASSLD: ERROR: ', 
     +            'DASSL can not handle JEVAL IERR' 
         stop 
      endif 
c compute pd = df/dy - cj*M 
      if (mljac.lt.neqn) then 
c df/dy banded and M banded 
         do 20 j=1,neqn 
            do 10 i=max(1,             mumas+1 +    1-j), 
     +              min(mlmas+mumas+1, mumas+1 + neqn-j) 
               pd((j-1)*ldjac+mljac+mujac-mumas+i) = 
     +         pd((j-1)*ldjac+mljac+mujac-mumas+i) - 
     +         cj*mas((j-1)*ldmas+i) 
   10       continue 
   20    continue 
      elseif (mlmas.lt.neqn) then 
c df/dy full and M banded 
         do 40 j=1,neqn 
            do 30 i=max(1,             mumas+1 +    1-j), 
     +              min(mlmas+mumas+1, mumas+1 + neqn-j) 
               pd((j-1)*ldjac+j-1+i-mumas) = 
     +         pd((j-1)*ldjac+j-1+i-mumas) - 
     +         cj*mas((j-1)*ldmas+i) 
   30       continue 
   40    continue 
      else 
c df/dy full and M full 
         do 60 j=1,neqn 
            do 50 i=1,neqn 
               pd((j-1)*ldjac+i) = 
     +         pd((j-1)*ldjac+i) - cj*mas((j-1)*ldmas+i) 
   50       continue 
   60    continue 
      endif 
      return 
      end 
