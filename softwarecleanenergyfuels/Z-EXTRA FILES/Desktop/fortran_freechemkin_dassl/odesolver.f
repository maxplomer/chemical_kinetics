      program odesolver
      EXTERNAL RES, JACDUM
      INTEGER NEQ,INFO(15),IDID,LRW,IWORK(100000),LIW,IPAR(1000),i
      DOUBLE PRECISION T, Y(54), YPRIME(54), TOUT, RTOL, ATOL,
     * RWORK(100000), RPAR(1000)

      NEQ = 54
      T=0
      Y(1:53) = 0.0
      Y(54) = 900.0
      Y(4) =  2.5742557936E-6
      Y(14) = 1.2871278968E-6
      Y(48) = 9.6792017840E-6
      TOUT=0.01
      info(1)=0
      info(2)=0
      info(3)=0
      info(4)=0
      info(5)=0
      info(6)=0
      info(7)=0
      info(8)=0
      info(9)=0
      info(10)=0
      info(11)=0
      info(12)=0
      info(13)=0
      info(14)=0
      info(15)=0
      RTOL = 0.00000001D0
      ATOL = 0.00000001D0
      LRW = 100000
      LIW = 100000
      idid=0

      DO 10 i = 1,801
        CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JACDUM)
        PRINT *, T, Y(54), Y(4), Y(48)
        TOUT = TOUT+0.01
  10  continue


      STOP
      end program odesolver



      SUBROUTINE JACDUM ()
      RETURN
      END



      SUBROUTINE  RES (T,Y,YDOT,DELTA,IRES,RPAR,IPAR)
      INTEGER IRES, IPAR
      DOUBLE PRECISION  T, Y(54), YDOT(54), DELTA(54), RPAR
      DOUBLE PRECISION  CV(53), U(53), C(53), WDOT(53)
      DOUBLE PRECISION Temp
      C=Y(1:53)
      Temp=Y(54)
      CALL GETCV(53,Temp,CV)
      CALL GETU(53,Temp,U)
      CALL GETWC(53,325,Temp,C,WDOT)
      DELTA(1:53) = WDOT - YDOT(1:53)
      DELTA(54) = YDOT(54) + DOT_PRODUCT(U,WDOT)/DOT_PRODUCT(C,CV)
      RETURN
      END






      SUBROUTINE GETH (KK,T,H)
C     enthalpy, h=f(T) for ideal gas, h=f(P,T) for real gas
C     units: erg/mol      erg=g*cm^2/s^2=1e-7 J
C     note: h=integral of cp dT from 298 to T    +     h_0,f,298
      INTEGER KK, i
      DOUBLE PRECISION T, T2, T3, T4, T5, RU
      DOUBLE PRECISION, dimension(KK) :: H, H_OVER_RU, COMMONT
      DOUBLE PRECISION, dimension(14,KK) :: N
      CALL GETRU (RU)
      CALL GETCOMMONT(KK, COMMONT)
      CALL GETTHERMODATA(KK,N)
      T2=T*T
      T3=T2*T
      T4=T3*T
      T5=T4*T
      do 10 i = 1, KK
        IF (T .GT. COMMONT(i)) THEN
          H_OVER_RU(i) = N(1,i)*T + N(2,i)*T2/2 + N(3,i)*T3/3 + 
     *                   N(4,i)*T4/4 + N(5,i)*T5/5 + N(6,i)
        ELSE
          H_OVER_RU(i) = N(8,i)*T + N(9,i)*T2/2 + N(10,i)*T3/3 + 
     *                   N(11,i)*T4/4 + N(12,i)*T5/5 + N(13,i)
        END IF
  10  continue
      H=H_OVER_RU*RU
      RETURN
      END

      SUBROUTINE GETG (KK,T,G)
C     gibbs free energy at standard state(p=1atm)
C     units: erg/mol
C     note: g=h - Ts
      INTEGER KK
      DOUBLE PRECISION RU, T
      DOUBLE PRECISION, dimension(KK) :: G, H, S
      CALL GETRU (RU)
      CALL GETS (KK,T,S)
      CALL GETH (KK,T,H)
      G = H - T*S
      RETURN
      END

      SUBROUTINE GETS (KK,T,S)
C     entropy at standard state, s=f(P,T)
C     units: erg/(mol*K)      erg=g*cm^2/s^2=1e-7 J
      INTEGER KK, i
      DOUBLE PRECISION LOGT, T, T2, T3, T4, RU
      DOUBLE PRECISION, dimension(KK) :: S, S_OVER_RU, COMMONT
      DOUBLE PRECISION, dimension(14,KK) :: N
      CALL GETRU (RU)
      CALL GETCOMMONT(KK, COMMONT)
      CALL GETTHERMODATA(KK,N)
      T2=T*T
      T3=T2*T
      T4=T3*T
      LOGT=LOG(T)
      do 10 i = 1, KK
        IF (T .GT. COMMONT(i)) THEN
          S_OVER_RU(i) = N(1,i)*logT + N(2,i)*T + N(3,i)*T2/2 +  
     *                   N(4,i)*T3/3 + N(5,i)*T4/4 + N(7,i)
        ELSE
          S_OVER_RU(i) = N(8,i)*logT + N(9,i)*T + N(10,i)*T2/2 + 
     *                   N(11,i)*T3/3 + N(12,i)*T4/4 + N(14,i)
        END IF
  10  continue
      S=S_OVER_RU*RU
      RETURN
      END

      SUBROUTINE GETU (KK,T,U)
C     units: erg/mol
C     note: h = u + pv    where v is volume/mol
C     note: ideal gas: p*v*M = RU*T      
C       where v is volume/mass  and M is mass/mol
      INTEGER KK
      DOUBLE PRECISION RU, T
      DOUBLE PRECISION, dimension(KK) :: U, H
      CALL GETRU (RU)
      CALL GETH (KK,T,H)
      U = H - RU*T
      RETURN
      END

      SUBROUTINE GETCV (KK,T,CV)
C     units: erg/(mol*K)
C     note: cp=cv + RU
      INTEGER KK
      DOUBLE PRECISION RU, T
      DOUBLE PRECISION, dimension(KK) :: CP, CV
      CALL GETRU (RU)
      CALL GETCP (KK,T,CP)
      CV=CP-RU
      RETURN
      END

      SUBROUTINE GETCP (KK,T,CP)
C     units: erg/(mol*K)
C     note: cp=cv + RU
      INTEGER KK, i
      DOUBLE PRECISION T, T2, T3, T4, RU
      DOUBLE PRECISION, dimension(KK) :: CP, CP_OVER_RU, COMMONT
      DOUBLE PRECISION, dimension(14,KK) :: N
      CALL GETRU (RU)
      CALL GETCOMMONT(KK, COMMONT)
      CALL GETTHERMODATA(KK,N)
      T2=T*T
      T3=T2*T
      T4=T3*T
      do 10 i = 1, KK
        IF (T .GT. COMMONT(i)) THEN
          CP_OVER_RU(i) = N(1,i) + N(2,i)*T + N(3,i)*T2 + 
     *                    N(4,i)*T3 + N(5,i)*T4
        ELSE
          CP_OVER_RU(i) = N(8,i) + N(9,i)*T + N(10,i)*T2 + 
     *                    N(11,i)*T3 + N(12,i)*T4
        END IF
  10  continue
      CP=CP_OVER_RU*RU
      RETURN
      END

      SUBROUTINE GETWC (KK,II,T,Y,WDOT)
C     units: mol/(cm^3*sec)
C     Y is mol concentration [mol/cm^3]
      INTEGER KK II
      DOUBLE PRECISION T
      INTEGER, dimension(KK,II) :: NU
      DOUBLE PRECISION, dimension(II) :: NETK,FWDK,REVK
      DOUBLE PRECISION, dimension(KK) :: WDOT, Y
      CALL GETKFKR (KK,II,T,Y,FWDK,REVK)
      CALL GETNUNET (KK,II,NU)
      NETK=FWDK-REVK
      WDOT=matmul(NU,NETK)
      RETURN
      END

      SUBROUTINE GETKFARRAY (II, T, KFARRAY)
      INTEGER II
      DOUBLE PRECISION RU, T
      DOUBLE PRECISION, dimension(II) ::  KFARRAY
      DOUBLE PRECISION, dimension(3,II) :: ABE
      CALL GETRU (RU)
      CALL GETABE (II, ABE)
C     units: 1 cal = 41 840 000 erg
      KFARRAY=(T**ABE(2,:))*EXP(ABE(1,:)-ABE(3,:)*41840000.0/(RU*T))
      RETURN
      END

      SUBROUTINE FINDELEMENT (EE, ELEMENT,I_ELEMENT)
      INTEGER EE, i, I_ELEMENT
      CHARACTER (LEN=3), dimension(EE) :: ELEMENTS
      CHARACTER (LEN=3) ELEMENT
      CALL GETELEMENTS (EE, ELEMENTS)
      do 10 i = 1, EE
        IF (ELEMENT .EQ. ELEMENTS(i)) THEN
          I_ELEMENT = i
          EXIT
        END IF
  10  continue
      RETURN
      END

      SUBROUTINE FINDSPECIES (KK, SPEC, I_SPECIES)
      INTEGER KK,i, I_SPECIES
      CHARACTER (LEN=16), dimension(KK) :: SPECIES
      CHARACTER (LEN=16) SPEC
      CALL GETSPECIES (KK, SPECIES)
      do 10 i = 1, KK
        IF (SPEC .EQ. SPECIES(i)) THEN
          I_SPECIES = i
          EXIT
        END IF
  10  continue
      RETURN
      END

      SUBROUTINE GETRU (RU)
C     units: erg/(mol*K)
      DOUBLE PRECISION RU
      RU = 83145100
      RETURN
      END

      SUBROUTINE GETPATM (PATM)
C     units: BA = (1/10) PA
      DOUBLE PRECISION PATM
      PATM = 1.01325E6
      RETURN
      END









      SUBROUTINE DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
C-----------------------------------------------------------------------
C NOTE:  Users of this solver, DDASSL, are encouraged to use the
C solver DDASPK instead.  DDASPK has a much improved initial condition
C calculation algorithm.  In addition, DDASPK includes iterative
C (Krylov) methods for the linear systems that arise, in addition to
C the direct (dense/banded) methods in DDASSL.
C-----------------------------------------------------------------------
C***BEGIN PROLOGUE  DDASSL
C***PURPOSE  This code solves a system of differential/algebraic
C            equations of the form G(T,Y,YPRIME) = 0.
C***LIBRARY   SLATEC (DASSL)
C***CATEGORY  I1A2
C***TYPE      DOUBLE PRECISION (SDASSL-S, DDASSL-D)
C***KEYWORDS  BACKWARD DIFFERENTIATION FORMULAS, DASSL,
C             DIFFERENTIAL/ALGEBRAIC, IMPLICIT DIFFERENTIAL SYSTEMS
C***AUTHOR  Petzold, Linda R., (LLNL)
C             Computing and Mathematics Research Division
C             Lawrence Livermore National Laboratory
C             L - 316, P.O. Box 808,
C             Livermore, CA.    94550
C***DESCRIPTION
C
C *Usage:
C
C      EXTERNAL RES, JAC
C      INTEGER NEQ, INFO(N), IDID, LRW, LIW, IWORK(LIW), IPAR
C      DOUBLE PRECISION T, Y(NEQ), YPRIME(NEQ), TOUT, RTOL, ATOL,
C     *   RWORK(LRW), RPAR
C
C      CALL DDASSL (RES, NEQ, T, Y, YPRIME, TOUT, INFO, RTOL, ATOL,
C     *   IDID, RWORK, LRW, IWORK, LIW, RPAR, IPAR, JAC)
C
C
C *Arguments:
C  (In the following, all real arrays should be type DOUBLE PRECISION.)
C
C  RES:EXT     This is a subroutine which you provide to define the
C              differential/algebraic system.
C
C  NEQ:IN      This is the number of equations to be solved.
C
C  T:INOUT     This is the current value of the independent variable.
C
C  Y(*):INOUT  This array contains the solution components at T.
C
C  YPRIME(*):INOUT  This array contains the derivatives of the solution
C              components at T.
C
C  TOUT:IN     This is a point at which a solution is desired.
C
C  INFO(N):IN  The basic task of the code is to solve the system from T
C              to TOUT and return an answer at TOUT.  INFO is an integer
C              array which is used to communicate exactly how you want
C              this task to be carried out.  (See below for details.)
C              N must be greater than or equal to 15.
C
C  RTOL,ATOL:INOUT  These quantities represent relative and absolute
C              error tolerances which you provide to indicate how
C              accurately you wish the solution to be computed.  You
C              may choose them to be both scalars or else both vectors.
C              Caution:  In Fortran 77, a scalar is not the same as an
C                        array of length 1.  Some compilers may object
C                        to using scalars for RTOL,ATOL.
C
C  IDID:OUT    This scalar quantity is an indicator reporting what the
C              code did.  You must monitor this integer variable to
C              decide  what action to take next.
C
C  RWORK:WORK  A real work array of length LRW which provides the
C              code with needed storage space.
C
C  LRW:IN      The length of RWORK.  (See below for required length.)
C
C  IWORK:WORK  An integer work array of length LIW which provides the
C              code with needed storage space.
C
C  LIW:IN      The length of IWORK.  (See below for required length.)
C
C  RPAR,IPAR:IN  These are real and integer parameter arrays which
C              you can use for communication between your calling
C              program and the RES subroutine (and the JAC subroutine)
C
C  JAC:EXT     This is the name of a subroutine which you may choose
C              to provide for defining a matrix of partial derivatives
C              described below.
C
C  Quantities which may be altered by DDASSL are:
C     T, Y(*), YPRIME(*), INFO(1), RTOL, ATOL,
C     IDID, RWORK(*) AND IWORK(*)
C
C *Description
C
C  Subroutine DDASSL uses the backward differentiation formulas of
C  orders one through five to solve a system of the above form for Y and
C  YPRIME.  Values for Y and YPRIME at the initial time must be given as
C  input.  These values must be consistent, (that is, if T,Y,YPRIME are
C  the given initial values, they must satisfy G(T,Y,YPRIME) = 0.).  The
C  subroutine solves the system from T to TOUT.  It is easy to continue
C  the solution to get results at additional TOUT.  This is the interval
C  mode of operation.  Intermediate results can also be obtained easily
C  by using the intermediate-output capability.
C
C  The following detailed description is divided into subsections:
C    1. Input required for the first call to DDASSL.
C    2. Output after any return from DDASSL.
C    3. What to do to continue the integration.
C    4. Error messages.
C
C
C  -------- INPUT -- WHAT TO DO ON THE FIRST CALL TO DDASSL ------------
C
C  The first call of the code is defined to be the start of each new
C  problem. Read through the descriptions of all the following items,
C  provide sufficient storage space for designated arrays, set
C  appropriate variables for the initialization of the problem, and
C  give information about how you want the problem to be solved.
C
C
C  RES -- Provide a subroutine of the form
C             SUBROUTINE RES(T,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
C         to define the system of differential/algebraic
C         equations which is to be solved. For the given values
C         of T,Y and YPRIME, the subroutine should
C         return the residual of the differential/algebraic
C         system
C             DELTA = G(T,Y,YPRIME)
C         (DELTA(*) is a vector of length NEQ which is
C         output for RES.)
C
C         Subroutine RES must not alter T,Y or YPRIME.
C         You must declare the name RES in an external
C         statement in your program that calls DDASSL.
C         You must dimension Y,YPRIME and DELTA in RES.
C
C         IRES is an integer flag which is always equal to
C         zero on input. Subroutine RES should alter IRES
C         only if it encounters an illegal value of Y or
C         a stop condition. Set IRES = -1 if an input value
C         is illegal, and DDASSL will try to solve the problem
C         without getting IRES = -1. If IRES = -2, DDASSL
C         will return control to the calling program
C         with IDID = -11.
C
C         RPAR and IPAR are real and integer parameter arrays which
C         you can use for communication between your calling program
C         and subroutine RES. They are not altered by DDASSL. If you
C         do not need RPAR or IPAR, ignore these parameters by treat-
C         ing them as dummy arguments. If you do choose to use them,
C         dimension them in your calling program and in RES as arrays
C         of appropriate length.
C
C  NEQ -- Set it to the number of differential equations.
C         (NEQ .GE. 1)
C
C  T -- Set it to the initial point of the integration.
C         T must be defined as a variable.
C
C  Y(*) -- Set this vector to the initial values of the NEQ solution
C         components at the initial point. You must dimension Y of
C         length at least NEQ in your calling program.
C
C  YPRIME(*) -- Set this vector to the initial values of the NEQ
C         first derivatives of the solution components at the initial
C         point.  You must dimension YPRIME at least NEQ in your
C         calling program. If you do not know initial values of some
C         of the solution components, see the explanation of INFO(11).
C
C  TOUT -- Set it to the first point at which a solution
C         is desired. You can not take TOUT = T.
C         integration either forward in T (TOUT .GT. T) or
C         backward in T (TOUT .LT. T) is permitted.
C
C         The code advances the solution from T to TOUT using
C         step sizes which are automatically selected so as to
C         achieve the desired accuracy. If you wish, the code will
C         return with the solution and its derivative at
C         intermediate steps (intermediate-output mode) so that
C         you can monitor them, but you still must provide TOUT in
C         accord with the basic aim of the code.
C
C         The first step taken by the code is a critical one
C         because it must reflect how fast the solution changes near
C         the initial point. The code automatically selects an
C         initial step size which is practically always suitable for
C         the problem. By using the fact that the code will not step
C         past TOUT in the first step, you could, if necessary,
C         restrict the length of the initial step size.
C
C         For some problems it may not be permissible to integrate
C         past a point TSTOP because a discontinuity occurs there
C         or the solution or its derivative is not defined beyond
C         TSTOP. When you have declared a TSTOP point (SEE INFO(4)
C         and RWORK(1)), you have told the code not to integrate
C         past TSTOP. In this case any TOUT beyond TSTOP is invalid
C         input.
C
C  INFO(*) -- Use the INFO array to give the code more details about
C         how you want your problem solved.  This array should be
C         dimensioned of length 15, though DDASSL uses only the first
C         eleven entries.  You must respond to all of the following
C         items, which are arranged as questions.  The simplest use
C         of the code corresponds to answering all questions as yes,
C         i.e. setting all entries of INFO to 0.
C
C       INFO(1) - This parameter enables the code to initialize
C              itself. You must set it to indicate the start of every
C              new problem.
C
C          **** Is this the first call for this problem ...
C                Yes - Set INFO(1) = 0
C                 No - Not applicable here.
C                      See below for continuation calls.  ****
C
C       INFO(2) - How much accuracy you want of your solution
C              is specified by the error tolerances RTOL and ATOL.
C              The simplest use is to take them both to be scalars.
C              To obtain more flexibility, they can both be vectors.
C              The code must be told your choice.
C
C          **** Are both error tolerances RTOL, ATOL scalars ...
C                Yes - Set INFO(2) = 0
C                      and input scalars for both RTOL and ATOL
C                 No - Set INFO(2) = 1
C                      and input arrays for both RTOL and ATOL ****
C
C       INFO(3) - The code integrates from T in the direction
C              of TOUT by steps. If you wish, it will return the
C              computed solution and derivative at the next
C              intermediate step (the intermediate-output mode) or
C              TOUT, whichever comes first. This is a good way to
C              proceed if you want to see the behavior of the solution.
C              If you must have solutions at a great many specific
C              TOUT points, this code will compute them efficiently.
C
C          **** Do you want the solution only at
C                TOUT (and not at the next intermediate step) ...
C                 Yes - Set INFO(3) = 0
C                  No - Set INFO(3) = 1 ****
C
C       INFO(4) - To handle solutions at a great many specific
C              values TOUT efficiently, this code may integrate past
C              TOUT and interpolate to obtain the result at TOUT.
C              Sometimes it is not possible to integrate beyond some
C              point TSTOP because the equation changes there or it is
C              not defined past TSTOP. Then you must tell the code
C              not to go past.
C
C           **** Can the integration be carried out without any
C                restrictions on the independent variable T ...
C                 Yes - Set INFO(4)=0
C                  No - Set INFO(4)=1
C                       and define the stopping point TSTOP by
C                       setting RWORK(1)=TSTOP ****
C
C       INFO(5) - To solve differential/algebraic problems it is
C              necessary to use a matrix of partial derivatives of the
C              system of differential equations. If you do not
C              provide a subroutine to evaluate it analytically (see
C              description of the item JAC in the call list), it will
C              be approximated by numerical differencing in this code.
C              although it is less trouble for you to have the code
C              compute partial derivatives by numerical differencing,
C              the solution will be more reliable if you provide the
C              derivatives via JAC. Sometimes numerical differencing
C              is cheaper than evaluating derivatives in JAC and
C              sometimes it is not - this depends on your problem.
C
C           **** Do you want the code to evaluate the partial
C                derivatives automatically by numerical differences ...
C                   Yes - Set INFO(5)=0
C                    No - Set INFO(5)=1
C                  and provide subroutine JAC for evaluating the
C                  matrix of partial derivatives ****
C
C       INFO(6) - DDASSL will perform much better if the matrix of
C              partial derivatives, DG/DY + CJ*DG/DYPRIME,
C              (here CJ is a scalar determined by DDASSL)
C              is banded and the code is told this. In this
C              case, the storage needed will be greatly reduced,
C              numerical differencing will be performed much cheaper,
C              and a number of important algorithms will execute much
C              faster. The differential equation is said to have
C              half-bandwidths ML (lower) and MU (upper) if equation i
C              involves only unknowns Y(J) with
C                             I-ML .LE. J .LE. I+MU
C              for all I=1,2,...,NEQ. Thus, ML and MU are the widths
C              of the lower and upper parts of the band, respectively,
C              with the main diagonal being excluded. If you do not
C              indicate that the equation has a banded matrix of partial
C              derivatives, the code works with a full matrix of NEQ**2
C              elements (stored in the conventional way). Computations
C              with banded matrices cost less time and storage than with
C              full matrices if 2*ML+MU .LT. NEQ. If you tell the
C              code that the matrix of partial derivatives has a banded
C              structure and you want to provide subroutine JAC to
C              compute the partial derivatives, then you must be careful
C              to store the elements of the matrix in the special form
C              indicated in the description of JAC.
C
C          **** Do you want to solve the problem using a full
C               (dense) matrix (and not a special banded
C               structure) ...
C                Yes - Set INFO(6)=0
C                 No - Set INFO(6)=1
C                       and provide the lower (ML) and upper (MU)
C                       bandwidths by setting
C                       IWORK(1)=ML
C                       IWORK(2)=MU ****
C
C
C        INFO(7) -- You can specify a maximum (absolute value of)
C              stepsize, so that the code
C              will avoid passing over very
C              large regions.
C
C          ****  Do you want the code to decide
C                on its own maximum stepsize?
C                Yes - Set INFO(7)=0
C                 No - Set INFO(7)=1
C                      and define HMAX by setting
C                      RWORK(2)=HMAX ****
C
C        INFO(8) -- Differential/algebraic problems
C              may occasionally suffer from
C              severe scaling difficulties on the
C              first step. If you know a great deal
C              about the scaling of your problem, you can
C              help to alleviate this problem by
C              specifying an initial stepsize HO.
C
C          ****  Do you want the code to define
C                its own initial stepsize?
C                Yes - Set INFO(8)=0
C                 No - Set INFO(8)=1
C                      and define HO by setting
C                      RWORK(3)=HO ****
C
C        INFO(9) -- If storage is a severe problem,
C              you can save some locations by
C              restricting the maximum order MAXORD.
C              the default value is 5. for each
C              order decrease below 5, the code
C              requires NEQ fewer locations, however
C              it is likely to be slower. In any
C              case, you must have 1 .LE. MAXORD .LE. 5
C          ****  Do you want the maximum order to
C                default to 5?
C                Yes - Set INFO(9)=0
C                 No - Set INFO(9)=1
C                      and define MAXORD by setting
C                      IWORK(3)=MAXORD ****
C
C        INFO(10) --If you know that the solutions to your equations
C               will always be nonnegative, it may help to set this
C               parameter. However, it is probably best to
C               try the code without using this option first,
C               and only to use this option if that doesn't
C               work very well.
C           ****  Do you want the code to solve the problem without
C                 invoking any special nonnegativity constraints?
C                  Yes - Set INFO(10)=0
C                   No - Set INFO(10)=1
C
C        INFO(11) --DDASSL normally requires the initial T,
C               Y, and YPRIME to be consistent. That is,
C               you must have G(T,Y,YPRIME) = 0 at the initial
C               time. If you do not know the initial
C               derivative precisely, you can let DDASSL try
C               to compute it.
C          ****   Are the initial T, Y, YPRIME consistent?
C                 Yes - Set INFO(11) = 0
C                  No - Set INFO(11) = 1,
C                       and set YPRIME to an initial approximation
C                       to YPRIME.  (If you have no idea what
C                       YPRIME should be, set it to zero. Note
C                       that the initial Y should be such
C                       that there must exist a YPRIME so that
C                       G(T,Y,YPRIME) = 0.)
C
C  RTOL, ATOL -- You must assign relative (RTOL) and absolute (ATOL
C         error tolerances to tell the code how accurately you
C         want the solution to be computed.  They must be defined
C         as variables because the code may change them.  You
C         have two choices --
C               Both RTOL and ATOL are scalars. (INFO(2)=0)
C               Both RTOL and ATOL are vectors. (INFO(2)=1)
C         in either case all components must be non-negative.
C
C         The tolerances are used by the code in a local error
C         test at each step which requires roughly that
C               ABS(LOCAL ERROR) .LE. RTOL*ABS(Y)+ATOL
C         for each vector component.
C         (More specifically, a root-mean-square norm is used to
C         measure the size of vectors, and the error test uses the
C         magnitude of the solution at the beginning of the step.)
C
C         The true (global) error is the difference between the
C         true solution of the initial value problem and the
C         computed approximation.  Practically all present day
C         codes, including this one, control the local error at
C         each step and do not even attempt to control the global
C         error directly.
C         Usually, but not always, the true accuracy of the
C         computed Y is comparable to the error tolerances. This
C         code will usually, but not always, deliver a more
C         accurate solution if you reduce the tolerances and
C         integrate again.  By comparing two such solutions you
C         can get a fairly reliable idea of the true error in the
C         solution at the bigger tolerances.
C
C         Setting ATOL=0. results in a pure relative error test on
C         that component.  Setting RTOL=0. results in a pure
C         absolute error test on that component.  A mixed test
C         with non-zero RTOL and ATOL corresponds roughly to a
C         relative error test when the solution component is much
C         bigger than ATOL and to an absolute error test when the
C         solution component is smaller than the threshhold ATOL.
C
C         The code will not attempt to compute a solution at an
C         accuracy unreasonable for the machine being used.  It will
C         advise you if you ask for too much accuracy and inform
C         you as to the maximum accuracy it believes possible.
C
C  RWORK(*) --  Dimension this real work array of length LRW in your
C         calling program.
C
C  LRW -- Set it to the declared length of the RWORK array.
C               You must have
C                    LRW .GE. 40+(MAXORD+4)*NEQ+NEQ**2
C               for the full (dense) JACOBIAN case (when INFO(6)=0), or
C                    LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
C               for the banded user-defined JACOBIAN case
C               (when INFO(5)=1 and INFO(6)=1), or
C                     LRW .GE. 40+(MAXORD+4)*NEQ+(2*ML+MU+1)*NEQ
C                           +2*(NEQ/(ML+MU+1)+1)
C               for the banded finite-difference-generated JACOBIAN case
C               (when INFO(5)=0 and INFO(6)=1)
C
C  IWORK(*) --  Dimension this integer work array of length LIW in
C         your calling program.
C
C  LIW -- Set it to the declared length of the IWORK array.
C               You must have LIW .GE. 20+NEQ
C
C  RPAR, IPAR -- These are parameter arrays, of real and integer
C         type, respectively.  You can use them for communication
C         between your program that calls DDASSL and the
C         RES subroutine (and the JAC subroutine).  They are not
C         altered by DDASSL.  If you do not need RPAR or IPAR,
C         ignore these parameters by treating them as dummy
C         arguments.  If you do choose to use them, dimension
C         them in your calling program and in RES (and in JAC)
C         as arrays of appropriate length.
C
C  JAC -- If you have set INFO(5)=0, you can ignore this parameter
C         by treating it as a dummy argument.  Otherwise, you must
C         provide a subroutine of the form
C             SUBROUTINE JAC(T,Y,YPRIME,PD,CJ,RPAR,IPAR)
C         to define the matrix of partial derivatives
C             PD=DG/DY+CJ*DG/DYPRIME
C         CJ is a scalar which is input to JAC.
C         For the given values of T,Y,YPRIME, the
C         subroutine must evaluate the non-zero partial
C         derivatives for each equation and each solution
C         component, and store these values in the
C         matrix PD.  The elements of PD are set to zero
C         before each call to JAC so only non-zero elements
C         need to be defined.
C
C         Subroutine JAC must not alter T,Y,(*),YPRIME(*), or CJ.
C         You must declare the name JAC in an EXTERNAL statement in
C         your program that calls DDASSL.  You must dimension Y,
C         YPRIME and PD in JAC.
C
C         The way you must store the elements into the PD matrix
C         depends on the structure of the matrix which you
C         indicated by INFO(6).
C               *** INFO(6)=0 -- Full (dense) matrix ***
C                   Give PD a first dimension of NEQ.
C                   When you evaluate the (non-zero) partial derivative
C                   of equation I with respect to variable J, you must
C                   store it in PD according to
C                   PD(I,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
C               *** INFO(6)=1 -- Banded JACOBIAN with ML lower and MU
C                   upper diagonal bands (refer to INFO(6) description
C                   of ML and MU) ***
C                   Give PD a first dimension of 2*ML+MU+1.
C                   when you evaluate the (non-zero) partial derivative
C                   of equation I with respect to variable J, you must
C                   store it in PD according to
C                   IROW = I - J + ML + MU + 1
C                   PD(IROW,J) = "DG(I)/DY(J)+CJ*DG(I)/DYPRIME(J)"
C
C         RPAR and IPAR are real and integer parameter arrays
C         which you can use for communication between your calling
C         program and your JACOBIAN subroutine JAC. They are not
C         altered by DDASSL. If you do not need RPAR or IPAR,
C         ignore these parameters by treating them as dummy
C         arguments. If you do choose to use them, dimension
C         them in your calling program and in JAC as arrays of
C         appropriate length.
C
C
C  OPTIONALLY REPLACEABLE NORM ROUTINE:
C
C     DDASSL uses a weighted norm DDANRM to measure the size
C     of vectors such as the estimated error in each step.
C     A FUNCTION subprogram
C       DOUBLE PRECISION FUNCTION DDANRM(NEQ,V,WT,RPAR,IPAR)
C       DIMENSION V(NEQ),WT(NEQ)
C     is used to define this norm. Here, V is the vector
C     whose norm is to be computed, and WT is a vector of
C     weights.  A DDANRM routine has been included with DDASSL
C     which computes the weighted root-mean-square norm
C     given by
C       DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
C     this norm is suitable for most problems. In some
C     special cases, it may be more convenient and/or
C     efficient to define your own norm by writing a function
C     subprogram to be called instead of DDANRM. This should,
C     however, be attempted only after careful thought and
C     consideration.
C
C
C  -------- OUTPUT -- AFTER ANY RETURN FROM DDASSL ---------------------
C
C  The principal aim of the code is to return a computed solution at
C  TOUT, although it is also possible to obtain intermediate results
C  along the way. To find out whether the code achieved its goal
C  or if the integration process was interrupted before the task was
C  completed, you must check the IDID parameter.
C
C
C  T -- The solution was successfully advanced to the
C               output value of T.
C
C  Y(*) -- Contains the computed solution approximation at T.
C
C  YPRIME(*) -- Contains the computed derivative
C               approximation at T.
C
C  IDID -- Reports what the code did.
C
C                     *** Task completed ***
C                Reported by positive values of IDID
C
C           IDID = 1 -- A step was successfully taken in the
C                   intermediate-output mode. The code has not
C                   yet reached TOUT.
C
C           IDID = 2 -- The integration to TSTOP was successfully
C                   completed (T=TSTOP) by stepping exactly to TSTOP.
C
C           IDID = 3 -- The integration to TOUT was successfully
C                   completed (T=TOUT) by stepping past TOUT.
C                   Y(*) is obtained by interpolation.
C                   YPRIME(*) is obtained by interpolation.
C
C                    *** Task interrupted ***
C                Reported by negative values of IDID
C
C           IDID = -1 -- A large amount of work has been expended.
C                   (About 500 steps)
C
C           IDID = -2 -- The error tolerances are too stringent.
C
C           IDID = -3 -- The local error test cannot be satisfied
C                   because you specified a zero component in ATOL
C                   and the corresponding computed solution
C                   component is zero. Thus, a pure relative error
C                   test is impossible for this component.
C
C           IDID = -6 -- DDASSL had repeated error test
C                   failures on the last attempted step.
C
C           IDID = -7 -- The corrector could not converge.
C
C           IDID = -8 -- The matrix of partial derivatives
C                   is singular.
C
C           IDID = -9 -- The corrector could not converge.
C                   there were repeated error test failures
C                   in this step.
C
C           IDID =-10 -- The corrector could not converge
C                   because IRES was equal to minus one.
C
C           IDID =-11 -- IRES equal to -2 was encountered
C                   and control is being returned to the
C                   calling program.
C
C           IDID =-12 -- DDASSL failed to compute the initial
C                   YPRIME.
C
C
C
C           IDID = -13,..,-32 -- Not applicable for this code
C
C                    *** Task terminated ***
C                Reported by the value of IDID=-33
C
C           IDID = -33 -- The code has encountered trouble from which
C                   it cannot recover. A message is printed
C                   explaining the trouble and control is returned
C                   to the calling program. For example, this occurs
C                   when invalid input is detected.
C
C  RTOL, ATOL -- These quantities remain unchanged except when
C               IDID = -2. In this case, the error tolerances have been
C               increased by the code to values which are estimated to
C               be appropriate for continuing the integration. However,
C               the reported solution at T was obtained using the input
C               values of RTOL and ATOL.
C
C  RWORK, IWORK -- Contain information which is usually of no
C               interest to the user but necessary for subsequent calls.
C               However, you may find use for
C
C               RWORK(3)--Which contains the step size H to be
C                       attempted on the next step.
C
C               RWORK(4)--Which contains the current value of the
C                       independent variable, i.e., the farthest point
C                       integration has reached. This will be different
C                       from T only when interpolation has been
C                       performed (IDID=3).
C
C               RWORK(7)--Which contains the stepsize used
C                       on the last successful step.
C
C               IWORK(7)--Which contains the order of the method to
C                       be attempted on the next step.
C
C               IWORK(8)--Which contains the order of the method used
C                       on the last step.
C
C               IWORK(11)--Which contains the number of steps taken so
C                        far.
C
C               IWORK(12)--Which contains the number of calls to RES
C                        so far.
C
C               IWORK(13)--Which contains the number of evaluations of
C                        the matrix of partial derivatives needed so
C                        far.
C
C               IWORK(14)--Which contains the total number
C                        of error test failures so far.
C
C               IWORK(15)--Which contains the total number
C                        of convergence test failures so far.
C                        (includes singular iteration matrix
C                        failures.)
C
C
C  -------- INPUT -- WHAT TO DO TO CONTINUE THE INTEGRATION ------------
C                    (CALLS AFTER THE FIRST)
C
C  This code is organized so that subsequent calls to continue the
C  integration involve little (if any) additional effort on your
C  part. You must monitor the IDID parameter in order to determine
C  what to do next.
C
C  Recalling that the principal task of the code is to integrate
C  from T to TOUT (the interval mode), usually all you will need
C  to do is specify a new TOUT upon reaching the current TOUT.
C
C  Do not alter any quantity not specifically permitted below,
C  in particular do not alter NEQ,T,Y(*),YPRIME(*),RWORK(*),IWORK(*)
C  or the differential equation in subroutine RES. Any such
C  alteration constitutes a new problem and must be treated as such,
C  i.e., you must start afresh.
C
C  You cannot change from vector to scalar error control or vice
C  versa (INFO(2)), but you can change the size of the entries of
C  RTOL, ATOL. Increasing a tolerance makes the equation easier
C  to integrate. Decreasing a tolerance will make the equation
C  harder to integrate and should generally be avoided.
C
C  You can switch from the intermediate-output mode to the
C  interval mode (INFO(3)) or vice versa at any time.
C
C  If it has been necessary to prevent the integration from going
C  past a point TSTOP (INFO(4), RWORK(1)), keep in mind that the
C  code will not integrate to any TOUT beyond the currently
C  specified TSTOP. Once TSTOP has been reached you must change
C  the value of TSTOP or set INFO(4)=0. You may change INFO(4)
C  or TSTOP at any time but you must supply the value of TSTOP in
C  RWORK(1) whenever you set INFO(4)=1.
C
C  Do not change INFO(5), INFO(6), IWORK(1), or IWORK(2)
C  unless you are going to restart the code.
C
C                 *** Following a completed task ***
C  If
C     IDID = 1, call the code again to continue the integration
C                  another step in the direction of TOUT.
C
C     IDID = 2 or 3, define a new TOUT and call the code again.
C                  TOUT must be different from T. You cannot change
C                  the direction of integration without restarting.
C
C                 *** Following an interrupted task ***
C               To show the code that you realize the task was
C               interrupted and that you want to continue, you
C               must take appropriate action and set INFO(1) = 1
C  If
C    IDID = -1, The code has taken about 500 steps.
C                  If you want to continue, set INFO(1) = 1 and
C                  call the code again. An additional 500 steps
C                  will be allowed.
C
C    IDID = -2, The error tolerances RTOL, ATOL have been
C                  increased to values the code estimates appropriate
C                  for continuing. You may want to change them
C                  yourself. If you are sure you want to continue
C                  with relaxed error tolerances, set INFO(1)=1 and
C                  call the code again.
C
C    IDID = -3, A solution component is zero and you set the
C                  corresponding component of ATOL to zero. If you
C                  are sure you want to continue, you must first
C                  alter the error criterion to use positive values
C                  for those components of ATOL corresponding to zero
C                  solution components, then set INFO(1)=1 and call
C                  the code again.
C
C    IDID = -4,-5  --- Cannot occur with this code.
C
C    IDID = -6, Repeated error test failures occurred on the
C                  last attempted step in DDASSL. A singularity in the
C                  solution may be present. If you are absolutely
C                  certain you want to continue, you should restart
C                  the integration. (Provide initial values of Y and
C                  YPRIME which are consistent)
C
C    IDID = -7, Repeated convergence test failures occurred
C                  on the last attempted step in DDASSL. An inaccurate
C                  or ill-conditioned JACOBIAN may be the problem. If
C                  you are absolutely certain you want to continue, you
C                  should restart the integration.
C
C    IDID = -8, The matrix of partial derivatives is singular.
C                  Some of your equations may be redundant.
C                  DDASSL cannot solve the problem as stated.
C                  It is possible that the redundant equations
C                  could be removed, and then DDASSL could
C                  solve the problem. It is also possible
C                  that a solution to your problem either
C                  does not exist or is not unique.
C
C    IDID = -9, DDASSL had multiple convergence test
C                  failures, preceded by multiple error
C                  test failures, on the last attempted step.
C                  It is possible that your problem
C                  is ill-posed, and cannot be solved
C                  using this code. Or, there may be a
C                  discontinuity or a singularity in the
C                  solution. If you are absolutely certain
C                  you want to continue, you should restart
C                  the integration.
C
C    IDID =-10, DDASSL had multiple convergence test failures
C                  because IRES was equal to minus one.
C                  If you are absolutely certain you want
C                  to continue, you should restart the
C                  integration.
C
C    IDID =-11, IRES=-2 was encountered, and control is being
C                  returned to the calling program.
C
C    IDID =-12, DDASSL failed to compute the initial YPRIME.
C                  This could happen because the initial
C                  approximation to YPRIME was not very good, or
C                  if a YPRIME consistent with the initial Y
C                  does not exist. The problem could also be caused
C                  by an inaccurate or singular iteration matrix.
C
C    IDID = -13,..,-32  --- Cannot occur with this code.
C
C
C                 *** Following a terminated task ***
C
C  If IDID= -33, you cannot continue the solution of this problem.
C                  An attempt to do so will result in your
C                  run being terminated.
C
C
C  -------- ERROR MESSAGES ---------------------------------------------
C
C      The SLATEC error print routine XERMSG is called in the event of
C   unsuccessful completion of a task.  Most of these are treated as
C   "recoverable errors", which means that (unless the user has directed
C   otherwise) control will be returned to the calling program for
C   possible action after the message has been printed.
C
C   In the event of a negative value of IDID other than -33, an appro-
C   priate message is printed and the "error number" printed by XERMSG
C   is the value of IDID.  There are quite a number of illegal input
C   errors that can lead to a returned value IDID=-33.  The conditions
C   and their printed "error numbers" are as follows:
C
C   Error number       Condition
C
C        1       Some element of INFO vector is not zero or one.
C        2       NEQ .le. 0
C        3       MAXORD not in range.
C        4       LRW is less than the required length for RWORK.
C        5       LIW is less than the required length for IWORK.
C        6       Some element of RTOL is .lt. 0
C        7       Some element of ATOL is .lt. 0
C        8       All elements of RTOL and ATOL are zero.
C        9       INFO(4)=1 and TSTOP is behind TOUT.
C       10       HMAX .lt. 0.0
C       11       TOUT is behind T.
C       12       INFO(8)=1 and H0=0.0
C       13       Some element of WT is .le. 0.0
C       14       TOUT is too close to T to start integration.
C       15       INFO(4)=1 and TSTOP is behind T.
C       16       --( Not used in this version )--
C       17       ML illegal.  Either .lt. 0 or .gt. NEQ
C       18       MU illegal.  Either .lt. 0 or .gt. NEQ
C       19       TOUT = T.
C
C   If DDASSL is called again without any action taken to remove the
C   cause of an unsuccessful return, XERMSG will be called with a fatal
C   error flag, which will cause unconditional termination of the
C   program.  There are two such fatal errors:
C
C   Error number -998:  The last step was terminated with a negative
C       value of IDID other than -33, and no appropriate action was
C       taken.
C
C   Error number -999:  The previous call was terminated because of
C       illegal input (IDID=-33) and there is illegal input in the
C       present call, as well.  (Suspect infinite loop.)
C
C  ---------------------------------------------------------------------
C
C***REFERENCES  A DESCRIPTION OF DASSL: A DIFFERENTIAL/ALGEBRAIC
C                 SYSTEM SOLVER, L. R. PETZOLD, SAND82-8637,
C                 SANDIA NATIONAL LABORATORIES, SEPTEMBER 1982.
C***ROUTINES CALLED  D1MACH, DDAINI, DDANRM, DDASTP, DDATRP, DDAWTS,
C                    XERMSG
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   880387  Code changes made.  All common statements have been
C           replaced by a DATA statement, which defines pointers into
C           RWORK, and PARAMETER statements which define pointers
C           into IWORK.  As well the documentation has gone through
C           grammatical changes.
C   881005  The prologue has been changed to mixed case.
C           The subordinate routines had revision dates changed to
C           this date, although the documentation for these routines
C           is all upper case.  No code changes.
C   890511  Code changes made.  The DATA statement in the declaration
C           section of DDASSL was replaced with a PARAMETER
C           statement.  Also the statement S = 100.D0 was removed
C           from the top of the Newton iteration in DDASTP.
C           The subordinate routines had revision dates changed to
C           this date.
C   890517  The revision date syntax was replaced with the revision
C           history syntax.  Also the "DECK" comment was added to
C           the top of all subroutines.  These changes are consistent
C           with new SLATEC guidelines.
C           The subordinate routines had revision dates changed to
C           this date.  No code changes.
C   891013  Code changes made.
C           Removed all occurrences of FLOAT or DBLE.  All operations
C           are now performed with "mixed-mode" arithmetic.
C           Also, specific function names were replaced with generic
C           function names to be consistent with new SLATEC guidelines.
C           In particular:
C              Replaced DSQRT with SQRT everywhere.
C              Replaced DABS with ABS everywhere.
C              Replaced DMIN1 with MIN everywhere.
C              Replaced MIN0 with MIN everywhere.
C              Replaced DMAX1 with MAX everywhere.
C              Replaced MAX0 with MAX everywhere.
C              Replaced DSIGN with SIGN everywhere.
C           Also replaced REVISION DATE with REVISION HISTORY in all
C           subordinate routines.
C   901004  Miscellaneous changes to prologue to complete conversion
C           to SLATEC 4.0 format.  No code changes.  (F.N.Fritsch)
C   901009  Corrected GAMS classification code and converted subsidiary
C           routines to 4.0 format.  No code changes.  (F.N.Fritsch)
C   901010  Converted XERRWV calls to XERMSG calls.  (R.Clemens, AFWL)
C   901019  Code changes made.
C           Merged SLATEC 4.0 changes with previous changes made
C           by C. Ulrich.  Below is a history of the changes made by
C           C. Ulrich. (Changes in subsidiary routines are implied
C           by this history)
C           891228  Bug was found and repaired inside the DDASSL
C                   and DDAINI routines.  DDAINI was incorrectly
C                   returning the initial T with Y and YPRIME
C                   computed at T+H.  The routine now returns T+H
C                   rather than the initial T.
C                   Cosmetic changes made to DDASTP.
C           900904  Three modifications were made to fix a bug (inside
C                   DDASSL) re interpolation for continuation calls and
C                   cases where TN is very close to TSTOP:
C
C                   1) In testing for whether H is too large, just
C                      compare H to (TSTOP - TN), rather than
C                      (TSTOP - TN) * (1-4*UROUND), and set H to
C                      TSTOP - TN.  This will force DDASTP to step
C                      exactly to TSTOP under certain situations
C                      (i.e. when H returned from DDASTP would otherwise
C                      take TN beyond TSTOP).
C
C                   2) Inside the DDASTP loop, interpolate exactly to
C                      TSTOP if TN is very close to TSTOP (rather than
C                      interpolating to within roundoff of TSTOP).
C
C                   3) Modified IDID description for IDID = 2 to say
C                      that the solution is returned by stepping exactly
C                      to TSTOP, rather than TOUT.  (In some cases the
C                      solution is actually obtained by extrapolating
C                      over a distance near unit roundoff to TSTOP,
C                      but this small distance is deemed acceptable in
C                      these circumstances.)
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue, removed unreferenced labels,
C           and improved XERMSG calls.  (FNF)
C   901030  Added ERROR MESSAGES section and reworked other sections to
C           be of more uniform format.  (FNF)
C   910624  Fixed minor bug related to HMAX (six lines after label
C           525).  (LRP)
C   000711  Fixed tests on (TN - TOUT) at 420 and 440 (ACH)
C***END PROLOGUE  DDASSL
C
C**End
C
C     Declare arguments.
C
      INTEGER  NEQ, INFO(15), IDID, LRW, IWORK(*), LIW, IPAR(*)
      DOUBLE PRECISION
     *   T, Y(*), YPRIME(*), TOUT, RTOL(*), ATOL(*), RWORK(*),
     *   RPAR(*)
      EXTERNAL  RES, JAC
C
C     Declare externals.
C
      EXTERNAL  D1MACH, DDAINI, DDANRM, DDASTP, DDATRP, DDAWTS, XERMSG
      DOUBLE PRECISION  D1MACH, DDANRM
C
C     Declare local variables.
C
      INTEGER  I, ITEMP, LALPHA, LBETA, LCJ, LCJOLD, LCTF, LDELTA,
     *   LENIW, LENPD, LENRW, LE, LETF, LGAMMA, LH, LHMAX, LHOLD, LIPVT,
     *   LJCALC, LK, LKOLD, LIWM, LML, LMTYPE, LMU, LMXORD, LNJE, LNPD,
     *   LNRE, LNS, LNST, LNSTL, LPD, LPHASE, LPHI, LPSI, LROUND, LS,
     *   LSIGMA, LTN, LTSTOP, LWM, LWT, MBAND, MSAVE, MXORD, NPD, NTEMP,
     *   NZFLG
      DOUBLE PRECISION
     *   ATOLI, H, HMAX, HMIN, HO, R, RH, RTOLI, TDIST, TN, TNEXT,
     *   TSTOP, UROUND, YPNORM
      LOGICAL  DONE
C       Auxiliary variables for conversion of values to be included in
C       error messages.
      CHARACTER*8  XERN1, XERN2
      CHARACTER*16 XERN3, XERN4
C
C     SET POINTERS INTO IWORK
      PARAMETER (LML=1, LMU=2, LMXORD=3, LMTYPE=4, LNST=11,
     *  LNRE=12, LNJE=13, LETF=14, LCTF=15, LNPD=16,
     *  LIPVT=21, LJCALC=5, LPHASE=6, LK=7, LKOLD=8,
     *  LNS=9, LNSTL=10, LIWM=1)
C
C     SET RELATIVE OFFSET INTO RWORK
      PARAMETER (NPD=1)
C
C     SET POINTERS INTO RWORK
      PARAMETER (LTSTOP=1, LHMAX=2, LH=3, LTN=4,
     *  LCJ=5, LCJOLD=6, LHOLD=7, LS=8, LROUND=9,
     *  LALPHA=11, LBETA=17, LGAMMA=23,
     *  LPSI=29, LSIGMA=35, LDELTA=41)
C
C***FIRST EXECUTABLE STATEMENT  DDASSL
      IF(INFO(1).NE.0)GO TO 100
C
C-----------------------------------------------------------------------
C     THIS BLOCK IS EXECUTED FOR THE INITIAL CALL ONLY.
C     IT CONTAINS CHECKING OF INPUTS AND INITIALIZATIONS.
C-----------------------------------------------------------------------
C
C     FIRST CHECK INFO ARRAY TO MAKE SURE ALL ELEMENTS OF INFO
C     ARE EITHER ZERO OR ONE.
      DO 10 I=2,11
         IF(INFO(I).NE.0.AND.INFO(I).NE.1)GO TO 701
10       CONTINUE
C
      IF(NEQ.LE.0)GO TO 702
C
C     CHECK AND COMPUTE MAXIMUM ORDER
      MXORD=5
      IF(INFO(9).EQ.0)GO TO 20
         MXORD=IWORK(LMXORD)
         IF(MXORD.LT.1.OR.MXORD.GT.5)GO TO 703
20       IWORK(LMXORD)=MXORD
C
C     COMPUTE MTYPE,LENPD,LENRW.CHECK ML AND MU.
      IF(INFO(6).NE.0)GO TO 40
         LENPD=NEQ**2
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
         IF(INFO(5).NE.0)GO TO 30
            IWORK(LMTYPE)=2
            GO TO 60
30          IWORK(LMTYPE)=1
            GO TO 60
40    IF(IWORK(LML).LT.0.OR.IWORK(LML).GE.NEQ)GO TO 717
      IF(IWORK(LMU).LT.0.OR.IWORK(LMU).GE.NEQ)GO TO 718
      LENPD=(2*IWORK(LML)+IWORK(LMU)+1)*NEQ
      IF(INFO(5).NE.0)GO TO 50
         IWORK(LMTYPE)=5
         MBAND=IWORK(LML)+IWORK(LMU)+1
         MSAVE=(NEQ/MBAND)+1
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD+2*MSAVE
         GO TO 60
50       IWORK(LMTYPE)=4
         LENRW=40+(IWORK(LMXORD)+4)*NEQ+LENPD
C
C     CHECK LENGTHS OF RWORK AND IWORK
60    LENIW=20+NEQ
      IWORK(LNPD)=LENPD
      IF(LRW.LT.LENRW)GO TO 704
      IF(LIW.LT.LENIW)GO TO 705
C
C     CHECK TO SEE THAT TOUT IS DIFFERENT FROM T
      IF(TOUT .EQ. T)GO TO 719
C
C     CHECK HMAX
      IF(INFO(7).EQ.0)GO TO 70
         HMAX=RWORK(LHMAX)
         IF(HMAX.LE.0.0D0)GO TO 710
70    CONTINUE
C
C     INITIALIZE COUNTERS
      IWORK(LNST)=0
      IWORK(LNRE)=0
      IWORK(LNJE)=0
C
      IWORK(LNSTL)=0
      IDID=1
      GO TO 200
C
C-----------------------------------------------------------------------
C     THIS BLOCK IS FOR CONTINUATION CALLS
C     ONLY. HERE WE CHECK INFO(1), AND IF THE
C     LAST STEP WAS INTERRUPTED WE CHECK WHETHER
C     APPROPRIATE ACTION WAS TAKEN.
C-----------------------------------------------------------------------
C
100   CONTINUE
      IF(INFO(1).EQ.1)GO TO 110
      IF(INFO(1).NE.-1)GO TO 701
C
C     IF WE ARE HERE, THE LAST STEP WAS INTERRUPTED
C     BY AN ERROR CONDITION FROM DDASTP, AND
C     APPROPRIATE ACTION WAS NOT TAKEN. THIS
C     IS A FATAL ERROR.
      WRITE (XERN1, '(I8)') IDID
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'THE LAST STEP TERMINATED WITH A NEGATIVE VALUE OF IDID = ' //
     *   XERN1 // ' AND NO APPROPRIATE ACTION WAS TAKEN.  ' //
     *   'RUN TERMINATED', -998, 2)
      RETURN
110   CONTINUE
      IWORK(LNSTL)=IWORK(LNST)
C
C-----------------------------------------------------------------------
C     THIS BLOCK IS EXECUTED ON ALL CALLS.
C     THE ERROR TOLERANCE PARAMETERS ARE
C     CHECKED, AND THE WORK ARRAY POINTERS
C     ARE SET.
C-----------------------------------------------------------------------
C
200   CONTINUE
C     CHECK RTOL,ATOL
      NZFLG=0
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 210 I=1,NEQ
         IF(INFO(2).EQ.1)RTOLI=RTOL(I)
         IF(INFO(2).EQ.1)ATOLI=ATOL(I)
         IF(RTOLI.GT.0.0D0.OR.ATOLI.GT.0.0D0)NZFLG=1
         IF(RTOLI.LT.0.0D0)GO TO 706
         IF(ATOLI.LT.0.0D0)GO TO 707
210      CONTINUE
      IF(NZFLG.EQ.0)GO TO 708
C
C     SET UP RWORK STORAGE.IWORK STORAGE IS FIXED
C     IN DATA STATEMENT.
      LE=LDELTA+NEQ
      LWT=LE+NEQ
      LPHI=LWT+NEQ
      LPD=LPHI+(IWORK(LMXORD)+1)*NEQ
      LWM=LPD
      NTEMP=NPD+IWORK(LNPD)
      IF(INFO(1).EQ.1)GO TO 400
C
C-----------------------------------------------------------------------
C     THIS BLOCK IS EXECUTED ON THE INITIAL CALL
C     ONLY. SET THE INITIAL STEP SIZE, AND
C     THE ERROR WEIGHT VECTOR, AND PHI.
C     COMPUTE INITIAL YPRIME, IF NECESSARY.
C-----------------------------------------------------------------------
C
      TN=T
      IDID=1
C
C     SET ERROR WEIGHT VECTOR WT
      CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,Y,RWORK(LWT),RPAR,IPAR)
      DO 305 I = 1,NEQ
         IF(RWORK(LWT+I-1).LE.0.0D0) GO TO 713
305      CONTINUE
C
C     COMPUTE UNIT ROUNDOFF AND HMIN
      UROUND = D1MACH(4)
      RWORK(LROUND) = UROUND
      HMIN = 4.0D0*UROUND*MAX(ABS(T),ABS(TOUT))
C
C     CHECK INITIAL INTERVAL TO SEE THAT IT IS LONG ENOUGH
      TDIST = ABS(TOUT - T)
      IF(TDIST .LT. HMIN) GO TO 714
C
C     CHECK HO, IF THIS WAS INPUT
      IF (INFO(8) .EQ. 0) GO TO 310
         HO = RWORK(LH)
         IF ((TOUT - T)*HO .LT. 0.0D0) GO TO 711
         IF (HO .EQ. 0.0D0) GO TO 712
         GO TO 320
310    CONTINUE
C
C     COMPUTE INITIAL STEPSIZE, TO BE USED BY EITHER
C     DDASTP OR DDAINI, DEPENDING ON INFO(11)
      HO = 0.001D0*TDIST
      YPNORM = DDANRM(NEQ,YPRIME,RWORK(LWT),RPAR,IPAR)
      IF (YPNORM .GT. 0.5D0/HO) HO = 0.5D0/YPNORM
      HO = SIGN(HO,TOUT-T)
C     ADJUST HO IF NECESSARY TO MEET HMAX BOUND
320   IF (INFO(7) .EQ. 0) GO TO 330
         RH = ABS(HO)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) HO = HO/RH
C     COMPUTE TSTOP, IF APPLICABLE
330   IF (INFO(4) .EQ. 0) GO TO 340
         TSTOP = RWORK(LTSTOP)
         IF ((TSTOP - T)*HO .LT. 0.0D0) GO TO 715
         IF ((T + HO - TSTOP)*HO .GT. 0.0D0) HO = TSTOP - T
         IF ((TSTOP - TOUT)*HO .LT. 0.0D0) GO TO 709
C
C     COMPUTE INITIAL DERIVATIVE, UPDATING TN AND Y, IF APPLICABLE
340   IF (INFO(11) .EQ. 0) GO TO 350
      CALL DDAINI(TN,Y,YPRIME,NEQ,
     *  RES,JAC,HO,RWORK(LWT),IDID,RPAR,IPAR,
     *  RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *  RWORK(LWM),IWORK(LIWM),HMIN,RWORK(LROUND),
     *  INFO(10),NTEMP)
      IF (IDID .LT. 0) GO TO 390
C
C     LOAD H WITH HO.  STORE H IN RWORK(LH)
350   H = HO
      RWORK(LH) = H
C
C     LOAD Y AND H*YPRIME INTO PHI(*,1) AND PHI(*,2)
      ITEMP = LPHI + NEQ
      DO 370 I = 1,NEQ
         RWORK(LPHI + I - 1) = Y(I)
370      RWORK(ITEMP + I - 1) = H*YPRIME(I)
C
390   GO TO 500
C
C-------------------------------------------------------
C     THIS BLOCK IS FOR CONTINUATION CALLS ONLY. ITS
C     PURPOSE IS TO CHECK STOP CONDITIONS BEFORE
C     TAKING A STEP.
C     ADJUST H IF NECESSARY TO MEET HMAX BOUND
C-------------------------------------------------------
C
400   CONTINUE
      UROUND=RWORK(LROUND)
      DONE = .FALSE.
      TN=RWORK(LTN)
      H=RWORK(LH)
      IF(INFO(7) .EQ. 0) GO TO 410
         RH = ABS(H)/RWORK(LHMAX)
         IF(RH .GT. 1.0D0) H = H/RH
410   CONTINUE
      IF(T .EQ. TOUT) GO TO 719
      IF((T - TOUT)*H .GT. 0.0D0) GO TO 711
      IF(INFO(4) .EQ. 1) GO TO 430
      IF(INFO(3) .EQ. 1) GO TO 420
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 490
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
420   IF((TN-T)*H .LE. 0.0D0) GO TO 490
      IF((TN - TOUT)*H .GE. 0.0D0) GO TO 425
      CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
425   CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
430   IF(INFO(3) .EQ. 1) GO TO 440
      TSTOP=RWORK(LTSTOP)
      IF((TN-TSTOP)*H.GT.0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H.LT.0.0D0)GO TO 709
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 450
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *   RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
440   TSTOP = RWORK(LTSTOP)
      IF((TN-TSTOP)*H .GT. 0.0D0) GO TO 715
      IF((TSTOP-TOUT)*H .LT. 0.0D0) GO TO 709
      IF((TN-T)*H .LE. 0.0D0) GO TO 450
      IF((TN - TOUT)*H .GE. 0.0D0) GO TO 445
      CALL DDATRP(TN,TN,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TN
      IDID = 1
      DONE = .TRUE.
      GO TO 490
445   CONTINUE
      CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      T = TOUT
      IDID = 3
      DONE = .TRUE.
      GO TO 490
450   CONTINUE
C     CHECK WHETHER WE ARE WITHIN ROUNDOFF OF TSTOP
      IF(ABS(TN-TSTOP).GT.100.0D0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 460
      CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,IWORK(LKOLD),
     *  RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      DONE = .TRUE.
      GO TO 490
460   TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 490
      H=TSTOP-TN
      RWORK(LH)=H
C
490   IF (DONE) GO TO 580
C
C-------------------------------------------------------
C     THE NEXT BLOCK CONTAINS THE CALL TO THE
C     ONE-STEP INTEGRATOR DDASTP.
C     THIS IS A LOOPING POINT FOR THE INTEGRATION STEPS.
C     CHECK FOR TOO MANY STEPS.
C     UPDATE WT.
C     CHECK FOR TOO MUCH ACCURACY REQUESTED.
C     COMPUTE MINIMUM STEPSIZE.
C-------------------------------------------------------
C
500   CONTINUE
C     CHECK FOR FAILURE TO COMPUTE INITIAL YPRIME
      IF (IDID .EQ. -12) GO TO 527
C
C     CHECK FOR TOO MANY STEPS
      IF((IWORK(LNST)-IWORK(LNSTL)).LT.500)
     *   GO TO 510
           IDID=-1
           GO TO 527
C
C     UPDATE WT
510   CALL DDAWTS(NEQ,INFO(2),RTOL,ATOL,RWORK(LPHI),
     *  RWORK(LWT),RPAR,IPAR)
      DO 520 I=1,NEQ
         IF(RWORK(I+LWT-1).GT.0.0D0)GO TO 520
           IDID=-3
           GO TO 527
520   CONTINUE
C
C     TEST FOR TOO MUCH ACCURACY REQUESTED.
      R=DDANRM(NEQ,RWORK(LPHI),RWORK(LWT),RPAR,IPAR)*
     *   100.0D0*UROUND
      IF(R.LE.1.0D0)GO TO 525
C     MULTIPLY RTOL AND ATOL BY R AND RETURN
      IF(INFO(2).EQ.1)GO TO 523
           RTOL(1)=R*RTOL(1)
           ATOL(1)=R*ATOL(1)
           IDID=-2
           GO TO 527
523   DO 524 I=1,NEQ
           RTOL(I)=R*RTOL(I)
524        ATOL(I)=R*ATOL(I)
      IDID=-2
      GO TO 527
525   CONTINUE
C
C     COMPUTE MINIMUM STEPSIZE
      HMIN=4.0D0*UROUND*MAX(ABS(TN),ABS(TOUT))
C
C     TEST H VS. HMAX
      IF (INFO(7) .NE. 0) THEN
         RH = ABS(H)/RWORK(LHMAX)
         IF (RH .GT. 1.0D0) H = H/RH
      ENDIF
C
      CALL DDASTP(TN,Y,YPRIME,NEQ,
     *   RES,JAC,H,RWORK(LWT),INFO(1),IDID,RPAR,IPAR,
     *   RWORK(LPHI),RWORK(LDELTA),RWORK(LE),
     *   RWORK(LWM),IWORK(LIWM),
     *   RWORK(LALPHA),RWORK(LBETA),RWORK(LGAMMA),
     *   RWORK(LPSI),RWORK(LSIGMA),
     *   RWORK(LCJ),RWORK(LCJOLD),RWORK(LHOLD),
     *   RWORK(LS),HMIN,RWORK(LROUND),
     *   IWORK(LPHASE),IWORK(LJCALC),IWORK(LK),
     *   IWORK(LKOLD),IWORK(LNS),INFO(10),NTEMP)
527   IF(IDID.LT.0)GO TO 600
C
C--------------------------------------------------------
C     THIS BLOCK HANDLES THE CASE OF A SUCCESSFUL RETURN
C     FROM DDASTP (IDID=1).  TEST FOR STOP CONDITIONS.
C--------------------------------------------------------
C
      IF(INFO(4).NE.0)GO TO 540
           IF(INFO(3).NE.0)GO TO 530
             IF((TN-TOUT)*H.LT.0.0D0)GO TO 500
             CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
530          IF((TN-TOUT)*H.GE.0.0D0)GO TO 535
             T=TN
             IDID=1
             GO TO 580
535          CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *         IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
             IDID=3
             T=TOUT
             GO TO 580
540   IF(INFO(3).NE.0)GO TO 550
      IF((TN-TOUT)*H.LT.0.0D0)GO TO 542
         CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *     IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
         T=TOUT
         IDID=3
         GO TO 580
542   IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*
     *   (ABS(TN)+ABS(H)))GO TO 545
      TNEXT=TN+H
      IF((TNEXT-TSTOP)*H.LE.0.0D0)GO TO 500
      H=TSTOP-TN
      GO TO 500
545   CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,
     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
550   IF((TN-TOUT)*H.GE.0.0D0)GO TO 555
      IF(ABS(TN-TSTOP).LE.100.0D0*UROUND*(ABS(TN)+ABS(H)))GO TO 552
      T=TN
      IDID=1
      GO TO 580
552   CALL DDATRP(TN,TSTOP,Y,YPRIME,NEQ,
     *  IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      IDID=2
      T=TSTOP
      GO TO 580
555   CALL DDATRP(TN,TOUT,Y,YPRIME,NEQ,
     *   IWORK(LKOLD),RWORK(LPHI),RWORK(LPSI))
      T=TOUT
      IDID=3
      GO TO 580
C
C--------------------------------------------------------
C     ALL SUCCESSFUL RETURNS FROM DDASSL ARE MADE FROM
C     THIS BLOCK.
C--------------------------------------------------------
C
580   CONTINUE
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C-----------------------------------------------------------------------
C     THIS BLOCK HANDLES ALL UNSUCCESSFUL
C     RETURNS OTHER THAN FOR ILLEGAL INPUT.
C-----------------------------------------------------------------------
C
600   CONTINUE
      ITEMP=-IDID
      GO TO (610,620,630,690,690,640,650,660,670,675,
     *  680,685), ITEMP
C
C     THE MAXIMUM NUMBER OF STEPS WAS TAKEN BEFORE
C     REACHING TOUT
610   WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT CURRENT T = ' // XERN3 // ' 500 STEPS TAKEN ON THIS ' //
     *   'CALL BEFORE REACHING TOUT', IDID, 1)
      GO TO 690
C
C     TOO MUCH ACCURACY FOR MACHINE PRECISION
620   WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' TOO MUCH ACCURACY REQUESTED FOR ' //
     *   'PRECISION OF MACHINE. RTOL AND ATOL WERE INCREASED TO ' //
     *   'APPROPRIATE VALUES', IDID, 1)
      GO TO 690
C
C     WT(I) .LE. 0.0 FOR SOME I (NOT AT START OF PROBLEM)
630   WRITE (XERN3, '(1P,D15.6)') TN
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' SOME ELEMENT OF WT HAS BECOME .LE. ' //
     *   '0.0', IDID, 1)
      GO TO 690
C
C     ERROR TEST FAILED REPEATEDLY OR WITH H=HMIN
640   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE ERROR TEST FAILED REPEATEDLY OR WITH ABS(H)=HMIN',
     *   IDID, 1)
      GO TO 690
C
C     CORRECTOR CONVERGENCE FAILED REPEATEDLY OR WITH H=HMIN
650   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE CORRECTOR FAILED TO CONVERGE REPEATEDLY OR WITH ' //
     *   'ABS(H)=HMIN', IDID, 1)
      GO TO 690
C
C     THE ITERATION MATRIX IS SINGULAR
660   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE ITERATION MATRIX IS SINGULAR', IDID, 1)
      GO TO 690
C
C     CORRECTOR FAILURE PRECEDED BY ERROR TEST FAILURES.
670   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE CORRECTOR COULD NOT CONVERGE.  ALSO, THE ERROR TEST ' //
     *   'FAILED REPEATEDLY.', IDID, 1)
      GO TO 690
C
C     CORRECTOR FAILURE BECAUSE IRES = -1
675   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE CORRECTOR COULD NOT CONVERGE BECAUSE IRES WAS EQUAL ' //
     *   'TO MINUS ONE', IDID, 1)
      GO TO 690
C
C     FAILURE BECAUSE IRES = -2
680   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') H
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' IRES WAS EQUAL TO MINUS TWO', IDID, 1)
      GO TO 690
C
C     FAILED TO COMPUTE INITIAL YPRIME
685   WRITE (XERN3, '(1P,D15.6)') TN
      WRITE (XERN4, '(1P,D15.6)') HO
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'AT T = ' // XERN3 // ' AND STEPSIZE H = ' // XERN4 //
     *   ' THE INITIAL YPRIME COULD NOT BE COMPUTED', IDID, 1)
      GO TO 690
C
690   CONTINUE
      INFO(1)=-1
      T=TN
      RWORK(LTN)=TN
      RWORK(LH)=H
      RETURN
C
C-----------------------------------------------------------------------
C     THIS BLOCK HANDLES ALL ERROR RETURNS DUE
C     TO ILLEGAL INPUT, AS DETECTED BEFORE CALLING
C     DDASTP. FIRST THE ERROR MESSAGE ROUTINE IS
C     CALLED. IF THIS HAPPENS TWICE IN
C     SUCCESSION, EXECUTION IS TERMINATED
C
C-----------------------------------------------------------------------
701   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'SOME ELEMENT OF INFO VECTOR IS NOT ZERO OR ONE', 1, 1)
      GO TO 750
C
702   WRITE (XERN1, '(I8)') NEQ
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'NEQ = ' // XERN1 // ' .LE. 0', 2, 1)
      GO TO 750
C
703   WRITE (XERN1, '(I8)') MXORD
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'MAXORD = ' // XERN1 // ' NOT IN RANGE', 3, 1)
      GO TO 750
C
704   WRITE (XERN1, '(I8)') LENRW
      WRITE (XERN2, '(I8)') LRW
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'RWORK LENGTH NEEDED, LENRW = ' // XERN1 //
     *   ', EXCEEDS LRW = ' // XERN2, 4, 1)
      GO TO 750
C
705   WRITE (XERN1, '(I8)') LENIW
      WRITE (XERN2, '(I8)') LIW
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'IWORK LENGTH NEEDED, LENIW = ' // XERN1 //
     *   ', EXCEEDS LIW = ' // XERN2, 5, 1)
      GO TO 750
C
706   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'SOME ELEMENT OF RTOL IS .LT. 0', 6, 1)
      GO TO 750
C
707   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'SOME ELEMENT OF ATOL IS .LT. 0', 7, 1)
      GO TO 750
C
708   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'ALL ELEMENTS OF RTOL AND ATOL ARE ZERO', 8, 1)
      GO TO 750
C
709   WRITE (XERN3, '(1P,D15.6)') TSTOP
      WRITE (XERN4, '(1P,D15.6)') TOUT
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'INFO(4) = 1 AND TSTOP = ' // XERN3 // ' BEHIND TOUT = ' //
     *   XERN4, 9, 1)
      GO TO 750
C
710   WRITE (XERN3, '(1P,D15.6)') HMAX
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'HMAX = ' // XERN3 // ' .LT. 0.0', 10, 1)
      GO TO 750
C
711   WRITE (XERN3, '(1P,D15.6)') TOUT
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'TOUT = ' // XERN3 // ' BEHIND T = ' // XERN4, 11, 1)
      GO TO 750
C
712   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'INFO(8)=1 AND H0=0.0', 12, 1)
      GO TO 750
C
713   CALL XERMSG ('SLATEC', 'DDASSL',
     *   'SOME ELEMENT OF WT IS .LE. 0.0', 13, 1)
      GO TO 750
C
714   WRITE (XERN3, '(1P,D15.6)') TOUT
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'TOUT = ' // XERN3 // ' TOO CLOSE TO T = ' // XERN4 //
     *   ' TO START INTEGRATION', 14, 1)
      GO TO 750
C
715   WRITE (XERN3, '(1P,D15.6)') TSTOP
      WRITE (XERN4, '(1P,D15.6)') T
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'INFO(4)=1 AND TSTOP = ' // XERN3 // ' BEHIND T = ' // XERN4,
     *   15, 1)
      GO TO 750
C
717   WRITE (XERN1, '(I8)') IWORK(LML)
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'ML = ' // XERN1 // ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',
     *   17, 1)
      GO TO 750
C
718   WRITE (XERN1, '(I8)') IWORK(LMU)
      CALL XERMSG ('SLATEC', 'DDASSL',
     *   'MU = ' // XERN1 // ' ILLEGAL.  EITHER .LT. 0 OR .GT. NEQ',
     *   18, 1)
      GO TO 750
C
719   WRITE (XERN3, '(1P,D15.6)') TOUT
      CALL XERMSG ('SLATEC', 'DDASSL',
     *  'TOUT = T = ' // XERN3, 19, 1)
      GO TO 750
C
750   IDID=-33
      IF(INFO(1).EQ.-1) THEN
         CALL XERMSG ('SLATEC', 'DDASSL',
     *      'REPEATED OCCURRENCES OF ILLEGAL INPUT$$' //
     *      'RUN TERMINATED. APPARENT INFINITE LOOP', -999, 2)
      ENDIF
C
      INFO(1)=-1
      RETURN
C-----------END OF SUBROUTINE DDASSL------------------------------------
      END

      SUBROUTINE DDAINI (X, Y, YPRIME, NEQ, RES, JAC, H, WT, IDID, RPAR,
     *   IPAR, PHI, DELTA, E, WM, IWM, HMIN, UROUND, NONNEG, NTEMP)
C***BEGIN PROLOGUE  DDAINI
C***SUBSIDIARY
C***PURPOSE  Initialization routine for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAINI-S, DDAINI-D)
C***AUTHOR  Petzold, Linda R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------
C     DDAINI TAKES ONE STEP OF SIZE H OR SMALLER
C     WITH THE BACKWARD EULER METHOD, TO
C     FIND YPRIME.  X AND Y ARE UPDATED TO BE CONSISTENT WITH THE
C     NEW STEP.  A MODIFIED DAMPED NEWTON ITERATION IS USED TO
C     SOLVE THE CORRECTOR ITERATION.
C
C     THE INITIAL GUESS FOR YPRIME IS USED IN THE
C     PREDICTION, AND IN FORMING THE ITERATION
C     MATRIX, BUT IS NOT INVOLVED IN THE
C     ERROR TEST. THIS MAY HAVE TROUBLE
C     CONVERGING IF THE INITIAL GUESS IS NO
C     GOOD, OR IF G(X,Y,YPRIME) DEPENDS
C     NONLINEARLY ON YPRIME.
C
C     THE PARAMETERS REPRESENT:
C     X --         INDEPENDENT VARIABLE
C     Y --         SOLUTION VECTOR AT X
C     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
C     NEQ --       NUMBER OF EQUATIONS
C     H --         STEPSIZE. IMDER MAY USE A STEPSIZE
C                  SMALLER THAN H.
C     WT --        VECTOR OF WEIGHTS FOR ERROR
C                  CRITERION
C     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS
C                  IDID= 1 -- YPRIME WAS FOUND SUCCESSFULLY
C                  IDID=-12 -- DDAINI FAILED TO FIND YPRIME
C     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS
C                  THAT ARE NOT ALTERED BY DDAINI
C     PHI --       WORK SPACE FOR DDAINI
C     DELTA,E --   WORK SPACE FOR DDAINI
C     WM,IWM --    REAL AND INTEGER ARRAYS STORING
C                  MATRIX INFORMATION
C
C-----------------------------------------------------------------
C***ROUTINES CALLED  DDAJAC, DDANRM, DDASLV
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C   901030  Minor corrections to declarations.  (FNF)
C***END PROLOGUE  DDAINI
C
      INTEGER  NEQ, IDID, IPAR(*), IWM(*), NONNEG, NTEMP
      DOUBLE PRECISION
     *   X, Y(*), YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),
     *   E(*), WM(*), HMIN, UROUND
      EXTERNAL  RES, JAC
C
      EXTERNAL  DDAJAC, DDANRM, DDASLV
      DOUBLE PRECISION  DDANRM
C
      INTEGER  I, IER, IRES, JCALC, LNJE, LNRE, M, MAXIT, MJAC, NCF,
     *   NEF, NSF
      DOUBLE PRECISION
     *   CJ, DAMP, DELNRM, ERR, OLDNRM, R, RATE, S, XOLD, YNORM
      LOGICAL  CONVGD
C
      PARAMETER (LNRE=12)
      PARAMETER (LNJE=13)
C
      DATA MAXIT/10/,MJAC/5/
      DATA DAMP/0.75D0/
C
C
C---------------------------------------------------
C     BLOCK 1.
C     INITIALIZATIONS.
C---------------------------------------------------
C
C***FIRST EXECUTABLE STATEMENT  DDAINI
      IDID=1
      NEF=0
      NCF=0
      NSF=0
      XOLD=X
      YNORM=DDANRM(NEQ,Y,WT,RPAR,IPAR)
C
C     SAVE Y AND YPRIME IN PHI
      DO 100 I=1,NEQ
         PHI(I,1)=Y(I)
100      PHI(I,2)=YPRIME(I)
C
C
C----------------------------------------------------
C     BLOCK 2.
C     DO ONE BACKWARD EULER STEP.
C----------------------------------------------------
C
C     SET UP FOR START OF CORRECTOR ITERATION
200   CJ=1.0D0/H
      X=X+H
C
C     PREDICT SOLUTION AND DERIVATIVE
      DO 250 I=1,NEQ
250     Y(I)=Y(I)+H*YPRIME(I)
C
      JCALC=-1
      M=0
      CONVGD=.TRUE.
C
C
C     CORRECTOR LOOP.
300   IWM(LNRE)=IWM(LNRE)+1
      IRES=0
C
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES.LT.0) GO TO 430
C
C
C     EVALUATE THE ITERATION MATRIX
      IF (JCALC.NE.-1) GO TO 310
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
     *   IER,WT,E,WM,IWM,RES,IRES,
     *   UROUND,JAC,RPAR,IPAR,NTEMP)
C
      S=1000000.D0
      IF (IRES.LT.0) GO TO 430
      IF (IER.NE.0) GO TO 430
      NSF=0
C
C
C
C     MULTIPLY RESIDUAL BY DAMPING FACTOR
310   CONTINUE
      DO 320 I=1,NEQ
320      DELTA(I)=DELTA(I)*DAMP
C
C     COMPUTE A NEW ITERATE (BACK SUBSTITUTION)
C     STORE THE CORRECTION IN DELTA
C
      CALL DDASLV(NEQ,DELTA,WM,IWM)
C
C     UPDATE Y AND YPRIME
      DO 330 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
330      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     TEST FOR CONVERGENCE OF THE ITERATION.
C
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.LE.100.D0*UROUND*YNORM)
     *   GO TO 400
C
      IF (M.GT.0) GO TO 340
         OLDNRM=DELNRM
         GO TO 350
C
340   RATE=(DELNRM/OLDNRM)**(1.0D0/M)
      IF (RATE.GT.0.90D0) GO TO 430
      S=RATE/(1.0D0-RATE)
C
350   IF (S*DELNRM .LE. 0.33D0) GO TO 400
C
C
C     THE CORRECTOR HAS NOT YET CONVERGED. UPDATE
C     M AND AND TEST WHETHER THE MAXIMUM
C     NUMBER OF ITERATIONS HAVE BEEN TRIED.
C     EVERY MJAC ITERATIONS, GET A NEW
C     ITERATION MATRIX.
C
      M=M+1
      IF (M.GE.MAXIT) GO TO 430
C
      IF ((M/MJAC)*MJAC.EQ.M) JCALC=-1
      GO TO 300
C
C
C     THE ITERATION HAS CONVERGED.
C     CHECK NONNEGATIVITY CONSTRAINTS
400   IF (NONNEG.EQ.0) GO TO 450
      DO 410 I=1,NEQ
410      DELTA(I)=MIN(Y(I),0.0D0)
C
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (DELNRM.GT.0.33D0) GO TO 430
C
      DO 420 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
420      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
      GO TO 450
C
C
C     EXITS FROM CORRECTOR LOOP.
430   CONVGD=.FALSE.
450   IF (.NOT.CONVGD) GO TO 600
C
C
C
C-----------------------------------------------------
C     BLOCK 3.
C     THE CORRECTOR ITERATION CONVERGED.
C     DO ERROR TEST.
C-----------------------------------------------------
C
      DO 510 I=1,NEQ
510      E(I)=Y(I)-PHI(I,1)
      ERR=DDANRM(NEQ,E,WT,RPAR,IPAR)
C
      IF (ERR.LE.1.0D0) RETURN
C
C
C
C--------------------------------------------------------
C     BLOCK 4.
C     THE BACKWARD EULER STEP FAILED. RESTORE X, Y
C     AND YPRIME TO THEIR ORIGINAL VALUES.
C     REDUCE STEPSIZE AND TRY AGAIN, IF
C     POSSIBLE.
C---------------------------------------------------------
C
600   CONTINUE
      X = XOLD
      DO 610 I=1,NEQ
         Y(I)=PHI(I,1)
610      YPRIME(I)=PHI(I,2)
C
      IF (CONVGD) GO TO 640
      IF (IER.EQ.0) GO TO 620
         NSF=NSF+1
         H=H*0.25D0
         IF (NSF.LT.3.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
620   IF (IRES.GT.-2) GO TO 630
         IDID=-12
         RETURN
630   NCF=NCF+1
      H=H*0.25D0
      IF (NCF.LT.10.AND.ABS(H).GE.HMIN) GO TO 690
         IDID=-12
         RETURN
C
640   NEF=NEF+1
      R=0.90D0/(2.0D0*ERR+0.0001D0)
      R=MAX(0.1D0,MIN(0.5D0,R))
      H=H*R
      IF (ABS(H).GE.HMIN.AND.NEF.LT.10) GO TO 690
         IDID=-12
         RETURN
690      GO TO 200
C
C-------------END OF SUBROUTINE DDAINI----------------------
      END

      SUBROUTINE DDAJAC (NEQ, X, Y, YPRIME, DELTA, CJ, H, IER, WT, E,
     *   WM, IWM, RES, IRES, UROUND, JAC, RPAR, IPAR, NTEMP)
C***BEGIN PROLOGUE  DDAJAC
C***SUBSIDIARY
C***PURPOSE  Compute the iteration matrix for DDASSL and form the
C            LU-decomposition.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAJAC-S, DDAJAC-D)
C***AUTHOR  Petzold, Linda R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS ROUTINE COMPUTES THE ITERATION MATRIX
C     PD=DG/DY+CJ*DG/DYPRIME (WHERE G(X,Y,YPRIME)=0).
C     HERE PD IS COMPUTED BY THE USER-SUPPLIED
C     ROUTINE JAC IF IWM(MTYPE) IS 1 OR 4, AND
C     IT IS COMPUTED BY NUMERICAL FINITE DIFFERENCING
C     IF IWM(MTYPE)IS 2 OR 5
C     THE PARAMETERS HAVE THE FOLLOWING MEANINGS.
C     Y        = ARRAY CONTAINING PREDICTED VALUES
C     YPRIME   = ARRAY CONTAINING PREDICTED DERIVATIVES
C     DELTA    = RESIDUAL EVALUATED AT (X,Y,YPRIME)
C                (USED ONLY IF IWM(MTYPE)=2 OR 5)
C     CJ       = SCALAR PARAMETER DEFINING ITERATION MATRIX
C     H        = CURRENT STEPSIZE IN INTEGRATION
C     IER      = VARIABLE WHICH IS .NE. 0
C                IF ITERATION MATRIX IS SINGULAR,
C                AND 0 OTHERWISE.
C     WT       = VECTOR OF WEIGHTS FOR COMPUTING NORMS
C     E        = WORK SPACE (TEMPORARY) OF LENGTH NEQ
C     WM       = REAL WORK SPACE FOR MATRICES. ON
C                OUTPUT IT CONTAINS THE LU DECOMPOSITION
C                OF THE ITERATION MATRIX.
C     IWM      = INTEGER WORK SPACE CONTAINING
C                MATRIX INFORMATION
C     RES      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
C                TO EVALUATE THE RESIDUAL FUNCTION G(X,Y,YPRIME)
C     IRES     = FLAG WHICH IS EQUAL TO ZERO IF NO ILLEGAL VALUES
C                IN RES, AND LESS THAN ZERO OTHERWISE.  (IF IRES
C                IS LESS THAN ZERO, THE MATRIX WAS NOT COMPLETED)
C                IN THIS CASE (IF IRES .LT. 0), THEN IER = 0.
C     UROUND   = THE UNIT ROUNDOFF ERROR OF THE MACHINE BEING USED.
C     JAC      = NAME OF THE EXTERNAL USER-SUPPLIED ROUTINE
C                TO EVALUATE THE ITERATION MATRIX (THIS ROUTINE
C                IS ONLY USED IF IWM(MTYPE) IS 1 OR 4)
C-----------------------------------------------------------------------
C***ROUTINES CALLED  DGBFA, DGEFA
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901010  Modified three MAX calls to be all on one line.  (FNF)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C   901101  Corrected PURPOSE.  (FNF)
C***END PROLOGUE  DDAJAC
C
      INTEGER  NEQ, IER, IWM(*), IRES, IPAR(*), NTEMP
      DOUBLE PRECISION
     *   X, Y(*), YPRIME(*), DELTA(*), CJ, H, WT(*), E(*), WM(*),
     *   UROUND, RPAR(*)
      EXTERNAL  RES, JAC
C
      EXTERNAL  DGBFA, DGEFA
C
      INTEGER  I, I1, I2, II, IPSAVE, ISAVE, J, K, L, LENPD, LIPVT,
     *   LML, LMTYPE, LMU, MBA, MBAND, MEB1, MEBAND, MSAVE, MTYPE, N,
     *   NPD, NPDM1, NROW
      DOUBLE PRECISION  DEL, DELINV, SQUR, YPSAVE, YSAVE
C
      PARAMETER (NPD=1)
      PARAMETER (LML=1)
      PARAMETER (LMU=2)
      PARAMETER (LMTYPE=4)
      PARAMETER (LIPVT=21)
C
C***FIRST EXECUTABLE STATEMENT  DDAJAC
      IER = 0
      NPDM1=NPD-1
      MTYPE=IWM(LMTYPE)
      GO TO (100,200,300,400,500),MTYPE
C
C
C     DENSE USER-SUPPLIED MATRIX
100   LENPD=NEQ*NEQ
      DO 110 I=1,LENPD
110      WM(NPDM1+I)=0.0D0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      GO TO 230
C
C
C     DENSE FINITE-DIFFERENCE-GENERATED MATRIX
200   IRES=0
      NROW=NPDM1
      SQUR = SQRT(UROUND)
      DO 210 I=1,NEQ
         DEL=SQUR*MAX(ABS(Y(I)),ABS(H*YPRIME(I)),ABS(WT(I)))
         DEL=SIGN(DEL,H*YPRIME(I))
         DEL=(Y(I)+DEL)-Y(I)
         YSAVE=Y(I)
         YPSAVE=YPRIME(I)
         Y(I)=Y(I)+DEL
         YPRIME(I)=YPRIME(I)+CJ*DEL
         CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
         IF (IRES .LT. 0) RETURN
         DELINV=1.0D0/DEL
         DO 220 L=1,NEQ
220      WM(NROW+L)=(E(L)-DELTA(L))*DELINV
      NROW=NROW+NEQ
      Y(I)=YSAVE
      YPRIME(I)=YPSAVE
210   CONTINUE
C
C
C     DO DENSE-MATRIX LU DECOMPOSITION ON PD
230      CALL DGEFA(WM(NPD),NEQ,NEQ,IWM(LIPVT),IER)
      RETURN
C
C
C     DUMMY SECTION FOR IWM(MTYPE)=3
300   RETURN
C
C
C     BANDED USER-SUPPLIED MATRIX
400   LENPD=(2*IWM(LML)+IWM(LMU)+1)*NEQ
      DO 410 I=1,LENPD
410      WM(NPDM1+I)=0.0D0
      CALL JAC(X,Y,YPRIME,WM(NPD),CJ,RPAR,IPAR)
      MEBAND=2*IWM(LML)+IWM(LMU)+1
      GO TO 550
C
C
C     BANDED FINITE-DIFFERENCE-GENERATED MATRIX
500   MBAND=IWM(LML)+IWM(LMU)+1
      MBA=MIN(MBAND,NEQ)
      MEBAND=MBAND+IWM(LML)
      MEB1=MEBAND-1
      MSAVE=(NEQ/MBAND)+1
      ISAVE=NTEMP-1
      IPSAVE=ISAVE+MSAVE
      IRES=0
      SQUR=SQRT(UROUND)
      DO 540 J=1,MBA
         DO 510 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          WM(ISAVE+K)=Y(N)
          WM(IPSAVE+K)=YPRIME(N)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          Y(N)=Y(N)+DEL
510       YPRIME(N)=YPRIME(N)+CJ*DEL
      CALL RES(X,Y,YPRIME,E,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) RETURN
      DO 530 N=J,NEQ,MBAND
          K= (N-J)/MBAND + 1
          Y(N)=WM(ISAVE+K)
          YPRIME(N)=WM(IPSAVE+K)
          DEL=SQUR*MAX(ABS(Y(N)),ABS(H*YPRIME(N)),ABS(WT(N)))
          DEL=SIGN(DEL,H*YPRIME(N))
          DEL=(Y(N)+DEL)-Y(N)
          DELINV=1.0D0/DEL
          I1=MAX(1,(N-IWM(LMU)))
          I2=MIN(NEQ,(N+IWM(LML)))
          II=N*MEB1-IWM(LML)+NPDM1
          DO 520 I=I1,I2
520         WM(II+I)=(E(I)-DELTA(I))*DELINV
530      CONTINUE
540   CONTINUE
C
C
C     DO LU DECOMPOSITION OF BANDED PD
550   CALL DGBFA(WM(NPD),MEBAND,NEQ,
     *    IWM(LML),IWM(LMU),IWM(LIPVT),IER)
      RETURN
C------END OF SUBROUTINE DDAJAC------
      END

      DOUBLE PRECISION FUNCTION DDANRM (NEQ, V, WT, RPAR, IPAR)
C***BEGIN PROLOGUE  DDANRM
C***SUBSIDIARY
C***PURPOSE  Compute vector norm for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDANRM-S, DDANRM-D)
C***AUTHOR  Petzold, Linda R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
C     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
C     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
C     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
C        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
C-----------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDANRM
C
      INTEGER  NEQ, IPAR(*)
      DOUBLE PRECISION  V(NEQ), WT(NEQ), RPAR(*)
C
      INTEGER  I
      DOUBLE PRECISION  SUM, VMAX
C
C***FIRST EXECUTABLE STATEMENT  DDANRM
      DDANRM = 0.0D0
      VMAX = 0.0D0
      DO 10 I = 1,NEQ
        IF(ABS(V(I)/WT(I)) .GT. VMAX) VMAX = ABS(V(I)/WT(I))
10      CONTINUE
      IF(VMAX .LE. 0.0D0) GO TO 30
      SUM = 0.0D0
      DO 20 I = 1,NEQ
20      SUM = SUM + ((V(I)/WT(I))/VMAX)**2
      DDANRM = VMAX*SQRT(SUM/NEQ)
30    CONTINUE
      RETURN
C------END OF FUNCTION DDANRM------
      END

      SUBROUTINE DDASLV (NEQ, DELTA, WM, IWM)
C***BEGIN PROLOGUE  DDASLV
C***SUBSIDIARY
C***PURPOSE  Linear system solver for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDASLV-S, DDASLV-D)
C***AUTHOR  Petzold, Linda R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR
C     SYSTEM ARISING IN THE NEWTON ITERATION.
C     MATRICES AND REAL TEMPORARY STORAGE AND
C     REAL INFORMATION ARE STORED IN THE ARRAY WM.
C     INTEGER MATRIX INFORMATION IS STORED IN
C     THE ARRAY IWM.
C     FOR A DENSE MATRIX, THE LINPACK ROUTINE
C     DGESL IS CALLED.
C     FOR A BANDED MATRIX,THE LINPACK ROUTINE
C     DGBSL IS CALLED.
C-----------------------------------------------------------------------
C***ROUTINES CALLED  DGBSL, DGESL
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDASLV
C
      INTEGER  NEQ, IWM(*)
      DOUBLE PRECISION  DELTA(*), WM(*)
C
      EXTERNAL  DGBSL, DGESL
C
      INTEGER  LIPVT, LML, LMU, LMTYPE, MEBAND, MTYPE, NPD
      PARAMETER (NPD=1)
      PARAMETER (LML=1)
      PARAMETER (LMU=2)
      PARAMETER (LMTYPE=4)
      PARAMETER (LIPVT=21)
C
C***FIRST EXECUTABLE STATEMENT  DDASLV
      MTYPE=IWM(LMTYPE)
      GO TO(100,100,300,400,400),MTYPE
C
C     DENSE MATRIX
100   CALL DGESL(WM(NPD),NEQ,NEQ,IWM(LIPVT),DELTA,0)
      RETURN
C
C     DUMMY SECTION FOR MTYPE=3
300   CONTINUE
      RETURN
C
C     BANDED MATRIX
400   MEBAND=2*IWM(LML)+IWM(LMU)+1
      CALL DGBSL(WM(NPD),MEBAND,NEQ,IWM(LML),
     *  IWM(LMU),IWM(LIPVT),DELTA,0)
      RETURN
C------END OF SUBROUTINE DDASLV------
      END

      SUBROUTINE DDASTP (X, Y, YPRIME, NEQ, RES, JAC, H, WT, JSTART,
     *   IDID, RPAR, IPAR, PHI, DELTA, E, WM, IWM, ALPHA, BETA, GAMMA,
     *   PSI, SIGMA, CJ, CJOLD, HOLD, S, HMIN, UROUND, IPHASE, JCALC, K,
     *   KOLD, NS, NONNEG, NTEMP)
C***BEGIN PROLOGUE  DDASTP
C***SUBSIDIARY
C***PURPOSE  Perform one step of the DDASSL integration.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDASTP-S, DDASTP-D)
C***AUTHOR  Petzold, Linda R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     DDASTP SOLVES A SYSTEM OF DIFFERENTIAL/
C     ALGEBRAIC EQUATIONS OF THE FORM
C     G(X,Y,YPRIME) = 0,  FOR ONE STEP (NORMALLY
C     FROM X TO X+H).
C
C     THE METHODS USED ARE MODIFIED DIVIDED
C     DIFFERENCE,FIXED LEADING COEFFICIENT
C     FORMS OF BACKWARD DIFFERENTIATION
C     FORMULAS. THE CODE ADJUSTS THE STEPSIZE
C     AND ORDER TO CONTROL THE LOCAL ERROR PER
C     STEP.
C
C
C     THE PARAMETERS REPRESENT
C     X  --        INDEPENDENT VARIABLE
C     Y  --        SOLUTION VECTOR AT X
C     YPRIME --    DERIVATIVE OF SOLUTION VECTOR
C                  AFTER SUCCESSFUL STEP
C     NEQ --       NUMBER OF EQUATIONS TO BE INTEGRATED
C     RES --       EXTERNAL USER-SUPPLIED SUBROUTINE
C                  TO EVALUATE THE RESIDUAL.  THE CALL IS
C                  CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
C                  X,Y,YPRIME ARE INPUT.  DELTA IS OUTPUT.
C                  ON INPUT, IRES=0.  RES SHOULD ALTER IRES ONLY
C                  IF IT ENCOUNTERS AN ILLEGAL VALUE OF Y OR A
C                  STOP CONDITION.  SET IRES=-1 IF AN INPUT VALUE
C                  OF Y IS ILLEGAL, AND DDASTP WILL TRY TO SOLVE
C                  THE PROBLEM WITHOUT GETTING IRES = -1.  IF
C                  IRES=-2, DDASTP RETURNS CONTROL TO THE CALLING
C                  PROGRAM WITH IDID = -11.
C     JAC --       EXTERNAL USER-SUPPLIED ROUTINE TO EVALUATE
C                  THE ITERATION MATRIX (THIS IS OPTIONAL)
C                  THE CALL IS OF THE FORM
C                  CALL JAC(X,Y,YPRIME,PD,CJ,RPAR,IPAR)
C                  PD IS THE MATRIX OF PARTIAL DERIVATIVES,
C                  PD=DG/DY+CJ*DG/DYPRIME
C     H --         APPROPRIATE STEP SIZE FOR NEXT STEP.
C                  NORMALLY DETERMINED BY THE CODE
C     WT --        VECTOR OF WEIGHTS FOR ERROR CRITERION.
C     JSTART --    INTEGER VARIABLE SET 0 FOR
C                  FIRST STEP, 1 OTHERWISE.
C     IDID --      COMPLETION CODE WITH THE FOLLOWING MEANINGS:
C                  IDID= 1 -- THE STEP WAS COMPLETED SUCCESSFULLY
C                  IDID=-6 -- THE ERROR TEST FAILED REPEATEDLY
C                  IDID=-7 -- THE CORRECTOR COULD NOT CONVERGE
C                  IDID=-8 -- THE ITERATION MATRIX IS SINGULAR
C                  IDID=-9 -- THE CORRECTOR COULD NOT CONVERGE.
C                             THERE WERE REPEATED ERROR TEST
C                             FAILURES ON THIS STEP.
C                  IDID=-10-- THE CORRECTOR COULD NOT CONVERGE
C                             BECAUSE IRES WAS EQUAL TO MINUS ONE
C                  IDID=-11-- IRES EQUAL TO -2 WAS ENCOUNTERED,
C                             AND CONTROL IS BEING RETURNED TO
C                             THE CALLING PROGRAM
C     RPAR,IPAR -- REAL AND INTEGER PARAMETER ARRAYS THAT
C                  ARE USED FOR COMMUNICATION BETWEEN THE
C                  CALLING PROGRAM AND EXTERNAL USER ROUTINES
C                  THEY ARE NOT ALTERED BY DDASTP
C     PHI --       ARRAY OF DIVIDED DIFFERENCES USED BY
C                  DDASTP. THE LENGTH IS NEQ*(K+1),WHERE
C                  K IS THE MAXIMUM ORDER
C     DELTA,E --   WORK VECTORS FOR DDASTP OF LENGTH NEQ
C     WM,IWM --    REAL AND INTEGER ARRAYS STORING
C                  MATRIX INFORMATION SUCH AS THE MATRIX
C                  OF PARTIAL DERIVATIVES,PERMUTATION
C                  VECTOR, AND VARIOUS OTHER INFORMATION.
C
C     THE OTHER PARAMETERS ARE INFORMATION
C     WHICH IS NEEDED INTERNALLY BY DDASTP TO
C     CONTINUE FROM STEP TO STEP.
C
C-----------------------------------------------------------------------
C***ROUTINES CALLED  DDAJAC, DDANRM, DDASLV, DDATRP
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C   951030  Reset PSI(1), PHI(*,2) at 690. (ACH)
C   000711  Fixed Newton convergence test below 360 (ACH)
C***END PROLOGUE  DDASTP
C
      INTEGER  NEQ, JSTART, IDID, IPAR(*), IWM(*), IPHASE, JCALC, K,
     *   KOLD, NS, NONNEG, NTEMP
      DOUBLE PRECISION
     *   X, Y(*), YPRIME(*), H, WT(*), RPAR(*), PHI(NEQ,*), DELTA(*),
     *   E(*), WM(*), ALPHA(*), BETA(*), GAMMA(*), PSI(*), SIGMA(*), CJ,
     *   CJOLD, HOLD, S, HMIN, UROUND
      EXTERNAL  RES, JAC
C
      EXTERNAL  DDAJAC, DDANRM, DDASLV, DDATRP
      DOUBLE PRECISION  DDANRM
C
      INTEGER  I, IER, IRES, J, J1, KDIFF, KM1, KNEW, KP1, KP2, LCTF,
     *   LETF, LMXORD, LNJE, LNRE, LNST, M, MAXIT, NCF, NEF, NSF, NSP1
      DOUBLE PRECISION
     *   ALPHA0, ALPHAS, CJLAST, CK, DELNRM, ENORM, ERK, ERKM1,
     *   ERKM2, ERKP1, ERR, EST, HNEW, OLDNRM, PNORM, R, RATE, TEMP1,
     *   TEMP2, TERK, TERKM1, TERKM2, TERKP1, XOLD, XRATE
      LOGICAL  CONVGD
C
      PARAMETER (LMXORD=3)
      PARAMETER (LNST=11)
      PARAMETER (LNRE=12)
      PARAMETER (LNJE=13)
      PARAMETER (LETF=14)
      PARAMETER (LCTF=15)
C
      DATA MAXIT/4/
      DATA XRATE/0.25D0/
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 1.
C     INITIALIZE. ON THE FIRST CALL,SET
C     THE ORDER TO 1 AND INITIALIZE
C     OTHER VARIABLES.
C-----------------------------------------------------------------------
C
C     INITIALIZATIONS FOR ALL CALLS
C***FIRST EXECUTABLE STATEMENT  DDASTP
      IDID=1
      XOLD=X
      NCF=0
      NSF=0
      NEF=0
      IF(JSTART .NE. 0) GO TO 120
C
C     IF THIS IS THE FIRST STEP,PERFORM
C     OTHER INITIALIZATIONS
      IWM(LETF) = 0
      IWM(LCTF) = 0
      K=1
      KOLD=0
      HOLD=0.0D0
      JSTART=1
      PSI(1)=H
      CJOLD = 1.0D0/H
      CJ = CJOLD
      S = 100.D0
      JCALC = -1
      DELNRM=1.0D0
      IPHASE = 0
      NS=0
120   CONTINUE
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 2
C     COMPUTE COEFFICIENTS OF FORMULAS FOR
C     THIS STEP.
C-----------------------------------------------------------------------
200   CONTINUE
      KP1=K+1
      KP2=K+2
      KM1=K-1
      XOLD=X
      IF(H.NE.HOLD.OR.K .NE. KOLD) NS = 0
      NS=MIN(NS+1,KOLD+2)
      NSP1=NS+1
      IF(KP1 .LT. NS)GO TO 230
C
      BETA(1)=1.0D0
      ALPHA(1)=1.0D0
      TEMP1=H
      GAMMA(1)=0.0D0
      SIGMA(1)=1.0D0
      DO 210 I=2,KP1
         TEMP2=PSI(I-1)
         PSI(I-1)=TEMP1
         BETA(I)=BETA(I-1)*PSI(I-1)/TEMP2
         TEMP1=TEMP2+H
         ALPHA(I)=H/TEMP1
         SIGMA(I)=(I-1)*SIGMA(I-1)*ALPHA(I)
         GAMMA(I)=GAMMA(I-1)+ALPHA(I-1)/H
210      CONTINUE
      PSI(KP1)=TEMP1
230   CONTINUE
C
C     COMPUTE ALPHAS, ALPHA0
      ALPHAS = 0.0D0
      ALPHA0 = 0.0D0
      DO 240 I = 1,K
        ALPHAS = ALPHAS - 1.0D0/I
        ALPHA0 = ALPHA0 - ALPHA(I)
240     CONTINUE
C
C     COMPUTE LEADING COEFFICIENT CJ
      CJLAST = CJ
      CJ = -ALPHAS/H
C
C     COMPUTE VARIABLE STEPSIZE ERROR COEFFICIENT CK
      CK = ABS(ALPHA(KP1) + ALPHAS - ALPHA0)
      CK = MAX(CK,ALPHA(KP1))
C
C     DECIDE WHETHER NEW JACOBIAN IS NEEDED
      TEMP1 = (1.0D0 - XRATE)/(1.0D0 + XRATE)
      TEMP2 = 1.0D0/TEMP1
      IF (CJ/CJOLD .LT. TEMP1 .OR. CJ/CJOLD .GT. TEMP2) JCALC = -1
      IF (CJ .NE. CJLAST) S = 100.D0
C
C     CHANGE PHI TO PHI STAR
      IF(KP1 .LT. NSP1) GO TO 280
      DO 270 J=NSP1,KP1
         DO 260 I=1,NEQ
260         PHI(I,J)=BETA(J)*PHI(I,J)
270      CONTINUE
280   CONTINUE
C
C     UPDATE TIME
      X=X+H
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 3
C     PREDICT THE SOLUTION AND DERIVATIVE,
C     AND SOLVE THE CORRECTOR EQUATION
C-----------------------------------------------------------------------
C
C     FIRST,PREDICT THE SOLUTION AND DERIVATIVE
300   CONTINUE
      DO 310 I=1,NEQ
         Y(I)=PHI(I,1)
310      YPRIME(I)=0.0D0
      DO 330 J=2,KP1
         DO 320 I=1,NEQ
            Y(I)=Y(I)+PHI(I,J)
320         YPRIME(I)=YPRIME(I)+GAMMA(J)*PHI(I,J)
330   CONTINUE
      PNORM = DDANRM (NEQ,Y,WT,RPAR,IPAR)
C
C
C
C     SOLVE THE CORRECTOR EQUATION USING A
C     MODIFIED NEWTON SCHEME.
      CONVGD= .TRUE.
      M=0
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
C
C
C     IF INDICATED,REEVALUATE THE
C     ITERATION MATRIX PD = DG/DY + CJ*DG/DYPRIME
C     (WHERE G(X,Y,YPRIME)=0). SET
C     JCALC TO 0 AS AN INDICATOR THAT
C     THIS HAS BEEN DONE.
      IF(JCALC .NE. -1)GO TO 340
      IWM(LNJE)=IWM(LNJE)+1
      JCALC=0
      CALL DDAJAC(NEQ,X,Y,YPRIME,DELTA,CJ,H,
     * IER,WT,E,WM,IWM,RES,IRES,UROUND,JAC,RPAR,
     * IPAR,NTEMP)
      CJOLD=CJ
      S = 100.D0
      IF (IRES .LT. 0) GO TO 380
      IF(IER .NE. 0)GO TO 380
      NSF=0
C
C
C     INITIALIZE THE ERROR ACCUMULATION VECTOR E.
340   CONTINUE
      DO 345 I=1,NEQ
345      E(I)=0.0D0
C
C
C     CORRECTOR LOOP.
350   CONTINUE
C
C     MULTIPLY RESIDUAL BY TEMP1 TO ACCELERATE CONVERGENCE
      TEMP1 = 2.0D0/(1.0D0 + CJ/CJOLD)
      DO 355 I = 1,NEQ
355     DELTA(I) = DELTA(I) * TEMP1
C
C     COMPUTE A NEW ITERATE (BACK-SUBSTITUTION).
C     STORE THE CORRECTION IN DELTA.
      CALL DDASLV(NEQ,DELTA,WM,IWM)
C
C     UPDATE Y, E, AND YPRIME
      DO 360 I=1,NEQ
         Y(I)=Y(I)-DELTA(I)
         E(I)=E(I)-DELTA(I)
360      YPRIME(I)=YPRIME(I)-CJ*DELTA(I)
C
C     TEST FOR CONVERGENCE OF THE ITERATION
      DELNRM=DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF (M .GT. 0) GO TO 365
         IF (DELNRM .LE. 100.D0*UROUND*PNORM) GO TO 375
         OLDNRM = DELNRM
         GO TO 367
365   RATE = (DELNRM/OLDNRM)**(1.0D0/M)
      IF (RATE .GT. 0.90D0) GO TO 370
      S = RATE/(1.0D0 - RATE)
367   IF (S*DELNRM .LE. 0.33D0) GO TO 375
C
C     THE CORRECTOR HAS NOT YET CONVERGED.
C     UPDATE M AND TEST WHETHER THE
C     MAXIMUM NUMBER OF ITERATIONS HAVE
C     BEEN TRIED.
      M=M+1
      IF(M.GE.MAXIT)GO TO 370
C
C     EVALUATE THE RESIDUAL
C     AND GO BACK TO DO ANOTHER ITERATION
      IWM(LNRE)=IWM(LNRE)+1
      IRES = 0
      CALL RES(X,Y,YPRIME,DELTA,IRES,
     *  RPAR,IPAR)
      IF (IRES .LT. 0) GO TO 380
      GO TO 350
C
C
C     THE CORRECTOR FAILED TO CONVERGE IN MAXIT
C     ITERATIONS. IF THE ITERATION MATRIX
C     IS NOT CURRENT,RE-DO THE STEP WITH
C     A NEW ITERATION MATRIX.
370   CONTINUE
      IF(JCALC.EQ.0)GO TO 380
      JCALC=-1
      GO TO 300
C
C
C     THE ITERATION HAS CONVERGED.  IF NONNEGATIVITY OF SOLUTION IS
C     REQUIRED, SET THE SOLUTION NONNEGATIVE, IF THE PERTURBATION
C     TO DO IT IS SMALL ENOUGH.  IF THE CHANGE IS TOO LARGE, THEN
C     CONSIDER THE CORRECTOR ITERATION TO HAVE FAILED.
375   IF(NONNEG .EQ. 0) GO TO 390
      DO 377 I = 1,NEQ
377      DELTA(I) = MIN(Y(I),0.0D0)
      DELNRM = DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      IF(DELNRM .GT. 0.33D0) GO TO 380
      DO 378 I = 1,NEQ
378      E(I) = E(I) - DELTA(I)
      GO TO 390
C
C
C     EXITS FROM BLOCK 3
C     NO CONVERGENCE WITH CURRENT ITERATION
C     MATRIX,OR SINGULAR ITERATION MATRIX
380   CONVGD= .FALSE.
390   JCALC = 1
      IF(.NOT.CONVGD)GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 4
C     ESTIMATE THE ERRORS AT ORDERS K,K-1,K-2
C     AS IF CONSTANT STEPSIZE WAS USED. ESTIMATE
C     THE LOCAL ERROR AT ORDER K AND TEST
C     WHETHER THE CURRENT STEP IS SUCCESSFUL.
C-----------------------------------------------------------------------
C
C     ESTIMATE ERRORS AT ORDERS K,K-1,K-2
      ENORM = DDANRM(NEQ,E,WT,RPAR,IPAR)
      ERK = SIGMA(K+1)*ENORM
      TERK = (K+1)*ERK
      EST = ERK
      KNEW=K
      IF(K .EQ. 1)GO TO 430
      DO 405 I = 1,NEQ
405     DELTA(I) = PHI(I,KP1) + E(I)
      ERKM1=SIGMA(K)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM1 = K*ERKM1
      IF(K .GT. 2)GO TO 410
      IF(TERKM1 .LE. 0.5D0*TERK)GO TO 420
      GO TO 430
410   CONTINUE
      DO 415 I = 1,NEQ
415     DELTA(I) = PHI(I,K) + DELTA(I)
      ERKM2=SIGMA(K-1)*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKM2 = (K-1)*ERKM2
      IF(MAX(TERKM1,TERKM2).GT.TERK)GO TO 430
C     LOWER THE ORDER
420   CONTINUE
      KNEW=K-1
      EST = ERKM1
C
C
C     CALCULATE THE LOCAL ERROR FOR THE CURRENT STEP
C     TO SEE IF THE STEP WAS SUCCESSFUL
430   CONTINUE
      ERR = CK * ENORM
      IF(ERR .GT. 1.0D0)GO TO 600
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 5
C     THE STEP IS SUCCESSFUL. DETERMINE
C     THE BEST ORDER AND STEPSIZE FOR
C     THE NEXT STEP. UPDATE THE DIFFERENCES
C     FOR THE NEXT STEP.
C-----------------------------------------------------------------------
      IDID=1
      IWM(LNST)=IWM(LNST)+1
      KDIFF=K-KOLD
      KOLD=K
      HOLD=H
C
C
C     ESTIMATE THE ERROR AT ORDER K+1 UNLESS:
C        ALREADY DECIDED TO LOWER ORDER, OR
C        ALREADY USING MAXIMUM ORDER, OR
C        STEPSIZE NOT CONSTANT, OR
C        ORDER RAISED IN PREVIOUS STEP
      IF(KNEW.EQ.KM1.OR.K.EQ.IWM(LMXORD))IPHASE=1
      IF(IPHASE .EQ. 0)GO TO 545
      IF(KNEW.EQ.KM1)GO TO 540
      IF(K.EQ.IWM(LMXORD)) GO TO 550
      IF(KP1.GE.NS.OR.KDIFF.EQ.1)GO TO 550
      DO 510 I=1,NEQ
510      DELTA(I)=E(I)-PHI(I,KP2)
      ERKP1 = (1.0D0/(K+2))*DDANRM(NEQ,DELTA,WT,RPAR,IPAR)
      TERKP1 = (K+2)*ERKP1
      IF(K.GT.1)GO TO 520
      IF(TERKP1.GE.0.5D0*TERK)GO TO 550
      GO TO 530
520   IF(TERKM1.LE.MIN(TERK,TERKP1))GO TO 540
      IF(TERKP1.GE.TERK.OR.K.EQ.IWM(LMXORD))GO TO 550
C
C     RAISE ORDER
530   K=KP1
      EST = ERKP1
      GO TO 550
C
C     LOWER ORDER
540   K=KM1
      EST = ERKM1
      GO TO 550
C
C     IF IPHASE = 0, INCREASE ORDER BY ONE AND MULTIPLY STEPSIZE BY
C     FACTOR TWO
545   K = KP1
      HNEW = H*2.0D0
      H = HNEW
      GO TO 575
C
C
C     DETERMINE THE APPROPRIATE STEPSIZE FOR
C     THE NEXT STEP.
550   HNEW=H
      TEMP2=K+1
      R=(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      IF(R .LT. 2.0D0) GO TO 555
      HNEW = 2.0D0*H
      GO TO 560
555   IF(R .GT. 1.0D0) GO TO 560
      R = MAX(0.5D0,MIN(0.9D0,R))
      HNEW = H*R
560   H=HNEW
C
C
C     UPDATE DIFFERENCES FOR NEXT STEP
575   CONTINUE
      IF(KOLD.EQ.IWM(LMXORD))GO TO 585
      DO 580 I=1,NEQ
580      PHI(I,KP2)=E(I)
585   CONTINUE
      DO 590 I=1,NEQ
590      PHI(I,KP1)=PHI(I,KP1)+E(I)
      DO 595 J1=2,KP1
         J=KP1-J1+1
         DO 595 I=1,NEQ
595      PHI(I,J)=PHI(I,J)+PHI(I,J+1)
      RETURN
C
C
C
C
C
C-----------------------------------------------------------------------
C     BLOCK 6
C     THE STEP IS UNSUCCESSFUL. RESTORE X,PSI,PHI
C     DETERMINE APPROPRIATE STEPSIZE FOR
C     CONTINUING THE INTEGRATION, OR EXIT WITH
C     AN ERROR FLAG IF THERE HAVE BEEN MANY
C     FAILURES.
C-----------------------------------------------------------------------
600   IPHASE = 1
C
C     RESTORE X,PHI,PSI
      X=XOLD
      IF(KP1.LT.NSP1)GO TO 630
      DO 620 J=NSP1,KP1
         TEMP1=1.0D0/BETA(J)
         DO 610 I=1,NEQ
610         PHI(I,J)=TEMP1*PHI(I,J)
620      CONTINUE
630   CONTINUE
      DO 640 I=2,KP1
640      PSI(I-1)=PSI(I)-H
C
C
C     TEST WHETHER FAILURE IS DUE TO CORRECTOR ITERATION
C     OR ERROR TEST
      IF(CONVGD)GO TO 660
      IWM(LCTF)=IWM(LCTF)+1
C
C
C     THE NEWTON ITERATION FAILED TO CONVERGE WITH
C     A CURRENT ITERATION MATRIX.  DETERMINE THE CAUSE
C     OF THE FAILURE AND TAKE APPROPRIATE ACTION.
      IF(IER.EQ.0)GO TO 650
C
C     THE ITERATION MATRIX IS SINGULAR. REDUCE
C     THE STEPSIZE BY A FACTOR OF 4. IF
C     THIS HAPPENS THREE TIMES IN A ROW ON
C     THE SAME STEP, RETURN WITH AN ERROR FLAG
      NSF=NSF+1
      R = 0.25D0
      H=H*R
      IF (NSF .LT. 3 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID=-8
      GO TO 675
C
C
C     THE NEWTON ITERATION FAILED TO CONVERGE FOR A REASON
C     OTHER THAN A SINGULAR ITERATION MATRIX.  IF IRES = -2, THEN
C     RETURN.  OTHERWISE, REDUCE THE STEPSIZE AND TRY AGAIN, UNLESS
C     TOO MANY FAILURES HAVE OCCURRED.
650   CONTINUE
      IF (IRES .GT. -2) GO TO 655
      IDID = -11
      GO TO 675
655   NCF = NCF + 1
      R = 0.25D0
      H = H*R
      IF (NCF .LT. 10 .AND. ABS(H) .GE. HMIN) GO TO 690
      IDID = -7
      IF (IRES .LT. 0) IDID = -10
      IF (NEF .GE. 3) IDID = -9
      GO TO 675
C
C
C     THE NEWTON SCHEME CONVERGED, AND THE CAUSE
C     OF THE FAILURE WAS THE ERROR ESTIMATE
C     EXCEEDING THE TOLERANCE.
660   NEF=NEF+1
      IWM(LETF)=IWM(LETF)+1
      IF (NEF .GT. 1) GO TO 665
C
C     ON FIRST ERROR TEST FAILURE, KEEP CURRENT ORDER OR LOWER
C     ORDER BY ONE.  COMPUTE NEW STEPSIZE BASED ON DIFFERENCES
C     OF THE SOLUTION.
      K = KNEW
      TEMP2 = K + 1
      R = 0.90D0*(2.0D0*EST+0.0001D0)**(-1.0D0/TEMP2)
      R = MAX(0.25D0,MIN(0.9D0,R))
      H = H*R
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     ON SECOND ERROR TEST FAILURE, USE THE CURRENT ORDER OR
C     DECREASE ORDER BY ONE.  REDUCE THE STEPSIZE BY A FACTOR OF
C     FOUR.
665   IF (NEF .GT. 2) GO TO 670
      K = KNEW
      R = 0.25D0
      H = R*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C     ON THIRD AND SUBSEQUENT ERROR TEST FAILURES, SET THE ORDER TO
C     ONE AND REDUCE THE STEPSIZE BY A FACTOR OF FOUR.
670   K = 1
      R = 0.25D0
      H = R*H
      IF (ABS(H) .GE. HMIN) GO TO 690
      IDID = -6
      GO TO 675
C
C
C     FOR ALL CRASHES, RESTORE Y TO ITS LAST VALUE,
C     INTERPOLATE TO FIND YPRIME AT LAST X, AND RETURN
675   CONTINUE
      CALL DDATRP(X,X,Y,YPRIME,NEQ,K,PHI,PSI)
      RETURN
C
C
C     GO BACK AND TRY THIS STEP AGAIN.
C     IF THIS IS THE FIRST STEP, RESET PSI(1) AND RESCALE PHI(*,2).
690   IF (KOLD .EQ. 0) THEN
        PSI(1) = H
        DO 695 I = 1,NEQ
695       PHI(I,2) = R*PHI(I,2)
        ENDIF
      GO TO 200
C
C------END OF SUBROUTINE DDASTP------
      END

      SUBROUTINE DDATRP (X, XOUT, YOUT, YPOUT, NEQ, KOLD, PHI, PSI)
C***BEGIN PROLOGUE  DDATRP
C***SUBSIDIARY
C***PURPOSE  Interpolation routine for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDATRP-S, DDATRP-D)
C***AUTHOR  Petzold, Linda R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THE METHODS IN SUBROUTINE DDASTP USE POLYNOMIALS
C     TO APPROXIMATE THE SOLUTION. DDATRP APPROXIMATES THE
C     SOLUTION AND ITS DERIVATIVE AT TIME XOUT BY EVALUATING
C     ONE OF THESE POLYNOMIALS, AND ITS DERIVATIVE,THERE.
C     INFORMATION DEFINING THIS POLYNOMIAL IS PASSED FROM
C     DDASTP, SO DDATRP CANNOT BE USED ALONE.
C
C     THE PARAMETERS ARE:
C     X     THE CURRENT TIME IN THE INTEGRATION.
C     XOUT  THE TIME AT WHICH THE SOLUTION IS DESIRED
C     YOUT  THE INTERPOLATED APPROXIMATION TO Y AT XOUT
C           (THIS IS OUTPUT)
C     YPOUT THE INTERPOLATED APPROXIMATION TO YPRIME AT XOUT
C           (THIS IS OUTPUT)
C     NEQ   NUMBER OF EQUATIONS
C     KOLD  ORDER USED ON LAST SUCCESSFUL STEP
C     PHI   ARRAY OF SCALED DIVIDED DIFFERENCES OF Y
C     PSI   ARRAY OF PAST STEPSIZE HISTORY
C-----------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDATRP
C
      INTEGER  NEQ, KOLD
      DOUBLE PRECISION  X, XOUT, YOUT(*), YPOUT(*), PHI(NEQ,*), PSI(*)
C
      INTEGER  I, J, KOLDP1
      DOUBLE PRECISION  C, D, GAMMA, TEMP1
C
C***FIRST EXECUTABLE STATEMENT  DDATRP
      KOLDP1=KOLD+1
      TEMP1=XOUT-X
      DO 10 I=1,NEQ
         YOUT(I)=PHI(I,1)
10       YPOUT(I)=0.0D0
      C=1.0D0
      D=0.0D0
      GAMMA=TEMP1/PSI(1)
      DO 30 J=2,KOLDP1
         D=D*GAMMA+C/PSI(J-1)
         C=C*GAMMA
         GAMMA=(TEMP1+PSI(J-1))/PSI(J)
         DO 20 I=1,NEQ
            YOUT(I)=YOUT(I)+C*PHI(I,J)
20          YPOUT(I)=YPOUT(I)+D*PHI(I,J)
30       CONTINUE
      RETURN
C
C------END OF SUBROUTINE DDATRP------
      END

      SUBROUTINE DDAWTS (NEQ, IWT, RTOL, ATOL, Y, WT, RPAR, IPAR)
C***BEGIN PROLOGUE  DDAWTS
C***SUBSIDIARY
C***PURPOSE  Set error weight vector for DDASSL.
C***LIBRARY   SLATEC (DASSL)
C***TYPE      DOUBLE PRECISION (SDAWTS-S, DDAWTS-D)
C***AUTHOR  Petzold, Linda R., (LLNL)
C***DESCRIPTION
C-----------------------------------------------------------------------
C     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
C     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
C     I=1,-,N.
C     RTOL AND ATOL ARE SCALARS IF IWT = 0,
C     AND VECTORS IF IWT = 1.
C-----------------------------------------------------------------------
C***ROUTINES CALLED  (NONE)
C***REVISION HISTORY  (YYMMDD)
C   830315  DATE WRITTEN
C   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
C   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
C   901026  Added explicit declarations for all variables and minor
C           cosmetic changes to prologue.  (FNF)
C***END PROLOGUE  DDAWTS
C
      INTEGER  NEQ, IWT, IPAR(*)
      DOUBLE PRECISION  RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*)
C
      INTEGER  I
      DOUBLE PRECISION  ATOLI, RTOLI
C
C***FIRST EXECUTABLE STATEMENT  DDAWTS
      RTOLI=RTOL(1)
      ATOLI=ATOL(1)
      DO 20 I=1,NEQ
         IF (IWT .EQ.0) GO TO 10
           RTOLI=RTOL(I)
           ATOLI=ATOL(I)
10         WT(I)=RTOLI*ABS(Y(I))+ATOLI
20         CONTINUE
      RETURN
C-----------END OF SUBROUTINE DDAWTS------------------------------------
      END

      DOUBLE PRECISION FUNCTION D1MACH (IDUMMY)
C***BEGIN PROLOGUE  D1MACH
C***PURPOSE  Compute the unit roundoff of the machine.
C***CATEGORY  R1
C***TYPE      DOUBLE PRECISION (R1MACH-S, D1MACH-D)
C***KEYWORDS  MACHINE CONSTANTS
C***AUTHOR  Hindmarsh, Alan C., (LLNL)
C***DESCRIPTION
C *Usage:
C        DOUBLE PRECISION  A, D1MACH
C        A = D1MACH(idummy)  [The argument is ignored.]
C
C *Function Return Values:
C     A : the unit roundoff of the machine.
C
C *Description:
C     The unit roundoff is defined as the smallest positive machine
C     number u such that  1.0 + u .ne. 1.0.  This is computed by D1MACH
C     in a machine-independent manner.
C
C***REFERENCES  (NONE)
C***ROUTINES CALLED  DUMSUM
C***REVISION HISTORY  (YYYYMMDD)
C   19930216  DATE WRITTEN
C   19930818  Added SLATEC-format prologue.  (FNF)
C   20030707  Added DUMSUM to force normal storage of COMP.  (ACH)
C***END PROLOGUE  D1MACH
C
      INTEGER IDUMMY
      DOUBLE PRECISION U, COMP
C***FIRST EXECUTABLE STATEMENT  D1MACH
      U = 1.0D0
 10   U = U*0.5D0
      CALL DUMSUM(1.0D0, U, COMP)
      IF (COMP .NE. 1.0D0) GO TO 10
      D1MACH = U*2.0D0
      RETURN
C----------------------- End of Function D1MACH ------------------------
      END
      SUBROUTINE DUMSUM(A,B,C)
C     Routine to force normal storing of A + B, for D1MACH.
      DOUBLE PRECISION A, B, C
      C = A + B
      RETURN
      END

      SUBROUTINE XERMSG (LIBRAR, SUBROU, MESSG, NERR, LEVEL)
C-----------------------------------------------------------------------
C Subroutines XERMSG, XSETF, XSETUN, and the function routine IXSAV, as
C given here, constitute a simplified version of the SLATEC error 
C handling package.  Written by A. C. Hindmarsh, 18 November 1992.
C
C All arguments are input arguments.
C LIBRAR = Library name (character array).  Prefixed to message.
C SUBROU = Routine name (character array).  Prefixed to message.
C MESSG  = The message (character array).
C NERR   = Integer error number.  Prefixed to message.
C LEVEL  = The error level..
C          0 or 1 means recoverable (control returns to caller).
C          2 means fatal (run is aborted--see note below).
C
C Note..  This routine has been simplified in the following ways..
C 1. A single prefix line is printed with NERR, SUBROU, and LIBRAR.
C 2. The message in MESSG is printed, unaltered, on lines of up to 72
C    characters each using a format of (A).
C 3. If LEVEL = 2, control passes to the statement   STOP
C    to abort the run.  This statement may be machine-dependent.
C
C For a different default logical unit number, change the data
C statement in function routine IXSAV.
C For a different run-abort command, change the statement following
C statement 100 at the end.
C-----------------------------------------------------------------------
C Subroutines called by XERMSG.. None
C Function routines called by XERMSG.. IXSAV
C Intrinsic function used by XERMSG.. LEN
C-----------------------------------------------------------------------
      CHARACTER*(*) LIBRAR, SUBROU, MESSG
      INTEGER NERR, LEVEL
      INTEGER I1, I2, IL, IXSAV, LENMSG, LLEN, LUNIT, MESFLG, NLINES
      PARAMETER (LLEN = 72)
C
C Get message print flag and logical unit number. ----------------------
      LUNIT = IXSAV (1, 0, .FALSE.)
      MESFLG = IXSAV (2, 0, .FALSE.)
      IF (MESFLG .EQ. 0) GO TO 100
C Write NERR, SUBROU, and LIBRAR. --------------------------------------
      I1 = LEN(SUBROU)
      I2 = LEN(LIBRAR)
      WRITE (LUNIT, 10) NERR, SUBROU(1:I1), LIBRAR(1:I2)
  10  FORMAT(/,'***Error number ',I6,' from ',A,' in library ',A,'***')
C Write the message. ---------------------------------------------------
      LENMSG = LEN(MESSG)
      NLINES = ( (LENMSG - 1)/LLEN ) + 1
      DO 20 IL = 1,NLINES
        I1 = 1 + (IL - 1)*LLEN
        I2 = MIN(IL*LLEN,LENMSG)
        WRITE (LUNIT,'(A)') MESSG(I1:I2)
  20    CONTINUE
C Abort the run if LEVEL = 2. ------------------------------------------
 100  IF (LEVEL .NE. 2) RETURN
      STOP
C----------------------- End of Subroutine XERMSG ----------------------
      END

      SUBROUTINE XSETUN (LUN)
C-----------------------------------------------------------------------
C This routine resets the logical unit number for messages.
C
C Subroutines called by XSETUN.. None
C Function routine called by XSETUN.. IXSAV
C-----------------------------------------------------------------------
      INTEGER LUN, JUNK, IXSAV
C
      IF (LUN .GT. 0) JUNK = IXSAV (1,LUN,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETUN ----------------------
      END

      SUBROUTINE XSETF (MFLAG)
C-----------------------------------------------------------------------
C This routine resets the print control flag MFLAG.
C
C Subroutines called by XSETF.. None
C Function routine called by XSETF.. IXSAV
C-----------------------------------------------------------------------
      INTEGER MFLAG, JUNK, IXSAV
C
      IF (MFLAG .EQ. 0 .OR. MFLAG .EQ. 1) JUNK = IXSAV (2,MFLAG,.TRUE.)
      RETURN
C----------------------- End of Subroutine XSETF -----------------------
      END

      INTEGER FUNCTION IXSAV (IPAR, IVALUE, ISET)
      LOGICAL ISET
      INTEGER IPAR, IVALUE
C-----------------------------------------------------------------------
C IXSAV saves and recalls one of two error message parameters:
C   LUNIT, the logical unit number to which messages are printed, and
C   MESFLG, the message print flag.
C This is a modification of the SLATEC library routine J4SAVE.
C
C Saved local variables..
C  LUNIT  = Logical unit number for messages.
C           The default is 6 (machine-dependent).
C  MESFLG = Print control flag..
C           1 means print all messages (the default).
C           0 means no printing.
C
C On input..
C   IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG).
C   IVALUE = The value to be set for the parameter, if ISET = .TRUE.
C   ISET   = Logical flag to indicate whether to read or write.
C            If ISET = .TRUE., the parameter will be given
C            the value IVALUE.  If ISET = .FALSE., the parameter
C            will be unchanged, and IVALUE is a dummy argument.
C
C On return..
C   IXSAV = The (old) value of the parameter.
C
C Subroutines/functions called by IXSAV.. None
C-----------------------------------------------------------------------
      INTEGER LUNIT, MESFLG
C-----------------------------------------------------------------------
C The following Fortran-77 declaration is to cause the values of the
C listed (local) variables to be saved between calls to this routine.
C-----------------------------------------------------------------------
      SAVE LUNIT, MESFLG
      DATA LUNIT/6/, MESFLG/1/
C
      IF (IPAR .EQ. 1) THEN
        IXSAV = LUNIT
        IF (ISET) LUNIT = IVALUE
        ENDIF
C
      IF (IPAR .EQ. 2) THEN
        IXSAV = MESFLG
        IF (ISET) MESFLG = IVALUE
        ENDIF
C
      RETURN
C----------------------- End of Function IXSAV -------------------------
      END




      subroutine daxpy(n,da,dx,incx,dy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),da
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (da .eq. 0.0d0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dy(iy) = dy(iy) + da*dx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dy(i) = dy(i) + da*dx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        dy(i) = dy(i) + da*dx(i)
        dy(i + 1) = dy(i + 1) + da*dx(i + 1)
        dy(i + 2) = dy(i + 2) + da*dx(i + 2)
        dy(i + 3) = dy(i + 3) + da*dx(i + 3)
   50 continue
      return
      end

      subroutine  dscal(n,da,dx,incx)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision da,dx(1)
      integer i,incx,m,mp1,n,nincx
c
      if( n.le.0 .or. incx.le.0 )return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        dx(i) = da*dx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dx(i) = da*dx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dx(i) = da*dx(i)
        dx(i + 1) = da*dx(i + 1)
        dx(i + 2) = da*dx(i + 2)
        dx(i + 3) = da*dx(i + 3)
        dx(i + 4) = da*dx(i + 4)
   50 continue
      return
      end

      integer function idamax(n,dx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c     modified 3/93 to return if incx .le. 0.
c
      double precision dx(1),dmax
      integer i,incx,ix,n
c
      idamax = 0
      if( n.lt.1 .or. incx.le.0 ) return
      idamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      dmax = dabs(dx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(dabs(dx(ix)).le.dmax) go to 5
         idamax = i
         dmax = dabs(dx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 dmax = dabs(dx(1))
      do 30 i = 2,n
         if(dabs(dx(i)).le.dmax) go to 30
         idamax = i
         dmax = dabs(dx(i))
   30 continue
      return
      end

      subroutine dgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      double precision a(lda,1)
c
c     dgefa factors a double precision matrix by gaussian elimination.
c
c     dgefa is usually called by dgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for dgeco) = (1 + 9/n)*(time for dgefa) .
c
c     on entry
c
c        a       double precision(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgesl or dgedi will divide by zero
c                     if called.  use  rcond  in dgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c
c     internal variables
c
      double precision t
      integer idamax,j,k,kp1,l,nm1
c
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = idamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (a(l,k) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -1.0d0/a(k,k)
            call dscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call daxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (a(n,n) .eq. 0.0d0) info = n
      return
      end






      subroutine dgbfa(abd,lda,n,ml,mu,ipvt,info)
      integer lda,n,ml,mu,ipvt(1),info
      double precision abd(lda,1)
c
c     dgbfa factors a double precision band matrix by elimination.
c
c     dgbfa is usually called by dgbco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c
c     on entry
c
c        abd     double precision(lda, n)
c                contains the matrix in band storage.  the columns
c                of the matrix are stored in the columns of  abd  and
c                the diagonals of the matrix are stored in rows
c                ml+1 through 2*ml+mu+1 of  abd .
c                see the comments below for details.
c
c        lda     integer
c                the leading dimension of the array  abd .
c                lda must be .ge. 2*ml + mu + 1 .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c                0 .le. ml .lt. n .
c
c        mu      integer
c                number of diagonals above the main diagonal.
c                0 .le. mu .lt. n .
c                more efficient if  ml .le. mu .
c     on return
c
c        abd     an upper triangular matrix in band storage and
c                the multipliers which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that dgbsl will divide by zero if
c                     called.  use  rcond  in dgbco for a reliable
c                     indication of singularity.
c
c     band storage
c
c           if  a  is a band matrix, the following program segment
c           will set up the input.
c
c                   ml = (band width below the diagonal)
c                   mu = (band width above the diagonal)
c                   m = ml + mu + 1
c                   do 20 j = 1, n
c                      i1 = max0(1, j-mu)
c                      i2 = min0(n, j+ml)
c                      do 10 i = i1, i2
c                         k = i - j + m
c                         abd(k,j) = a(i,j)
c                10    continue
c                20 continue
c
c           this uses rows  ml+1  through  2*ml+mu+1  of  abd .
c           in addition, the first  ml  rows in  abd  are used for
c           elements generated during the triangularization.
c           the total number of rows needed in  abd  is  2*ml+mu+1 .
c           the  ml+mu by ml+mu  upper left triangle and the
c           ml by ml  lower right triangle are not referenced.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,dscal,idamax
c     fortran max0,min0
c
c     internal variables
c
      double precision t
      integer i,idamax,i0,j,ju,jz,j0,j1,k,kp1,l,lm,m,mm,nm1
c
c
      m = ml + mu + 1
      info = 0
c
c     zero initial fill-in columns
c
      j0 = mu + 2
      j1 = min0(n,m) - 1
      if (j1 .lt. j0) go to 30
      do 20 jz = j0, j1
         i0 = m + 1 - jz
         do 10 i = i0, ml
            abd(i,jz) = 0.0d0
   10    continue
   20 continue
   30 continue
      jz = j1
      ju = 0
c
c     gaussian elimination with partial pivoting
c
      nm1 = n - 1
      if (nm1 .lt. 1) go to 130
      do 120 k = 1, nm1
         kp1 = k + 1
c
c        zero next fill-in column
c
         jz = jz + 1
         if (jz .gt. n) go to 50
         if (ml .lt. 1) go to 50
            do 40 i = 1, ml
               abd(i,jz) = 0.0d0
   40       continue
   50    continue
c
c        find l = pivot index
c
         lm = min0(ml,n-k)
         l = idamax(lm+1,abd(m,k),1) + m - 1
         ipvt(k) = l + k - m
c
c        zero pivot implies this column already triangularized
c
         if (abd(l,k) .eq. 0.0d0) go to 100
c
c           interchange if necessary
c
            if (l .eq. m) go to 60
               t = abd(l,k)
               abd(l,k) = abd(m,k)
               abd(m,k) = t
   60       continue
c
c           compute multipliers
c
            t = -1.0d0/abd(m,k)
            call dscal(lm,t,abd(m+1,k),1)
c
c           row elimination with column indexing
c
            ju = min0(max0(ju,mu+ipvt(k)),n)
            mm = m
            if (ju .lt. kp1) go to 90
            do 80 j = kp1, ju
               l = l - 1
               mm = mm - 1
               t = abd(l,j)
               if (l .eq. mm) go to 70
                  abd(l,j) = abd(mm,j)
                  abd(mm,j) = t
   70          continue
               call daxpy(lm,t,abd(m+1,k),1,abd(mm+1,j),1)
   80       continue
   90       continue
         go to 110
  100    continue
            info = k
  110    continue
  120 continue
  130 continue
      ipvt(n) = n
      if (abd(m,n) .eq. 0.0d0) info = n
      return
      end


      double precision function ddot(n,dx,incx,dy,incy)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision dx(1),dy(1),dtemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      ddot = 0.0d0
      dtemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        dtemp = dtemp + dx(ix)*dy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      ddot = dtemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        dtemp = dtemp + dx(i)*dy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        dtemp = dtemp + dx(i)*dy(i) + dx(i + 1)*dy(i + 1) +
     *   dx(i + 2)*dy(i + 2) + dx(i + 3)*dy(i + 3) + dx(i + 4)*dy(i + 4)
   50 continue
   60 ddot = dtemp
      return
      end

      subroutine dgesl(a,lda,n,ipvt,b,job)
      integer lda,n,ipvt(1),job
      double precision a(lda,1),b(1)
c
c     dgesl solves the double precision system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgeco or dgefa.
c
c     on entry
c
c        a       double precision(lda, n)
c                the output from dgeco or dgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from dgeco or dgefa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b  where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgeco has set rcond .gt. 0.0
c        or dgefa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgeco(a,lda,n,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgesl(a,lda,n,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,nm1
c
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve  l*y = b
c
         if (nm1 .lt. 1) go to 30
         do 20 k = 1, nm1
            l = ipvt(k)
            t = b(l)
            if (l .eq. k) go to 10
               b(l) = b(k)
               b(k) = t
   10       continue
            call daxpy(n-k,t,a(k+1,k),1,b(k+1),1)
   20    continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/a(k,k)
            t = -b(k)
            call daxpy(k-1,t,a(1,k),1,b(1),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            t = ddot(k-1,a(1,k),1,b(1),1)
            b(k) = (b(k) - t)/a(k,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .lt. 1) go to 90
         do 80 kb = 1, nm1
            k = n - kb
            b(k) = b(k) + ddot(n-k,a(k+1,k),1,b(k+1),1)
            l = ipvt(k)
            if (l .eq. k) go to 70
               t = b(l)
               b(l) = b(k)
               b(k) = t
   70       continue
   80    continue
   90    continue
  100 continue
      return
      end



      subroutine dgbsl(abd,lda,n,ml,mu,ipvt,b,job)
      integer lda,n,ml,mu,ipvt(1),job
      double precision abd(lda,1),b(1)
c
c     dgbsl solves the double precision band system
c     a * x = b  or  trans(a) * x = b
c     using the factors computed by dgbco or dgbfa.
c
c     on entry
c
c        abd     double precision(lda, n)
c                the output from dgbco or dgbfa.
c
c        lda     integer
c                the leading dimension of the array  abd .
c
c        n       integer
c                the order of the original matrix.
c
c        ml      integer
c                number of diagonals below the main diagonal.
c
c        mu      integer
c                number of diagonals above the main diagonal.
c
c        ipvt    integer(n)
c                the pivot vector from dgbco or dgbfa.
c
c        b       double precision(n)
c                the right hand side vector.
c
c        job     integer
c                = 0         to solve  a*x = b ,
c                = nonzero   to solve  trans(a)*x = b , where
c                            trans(a)  is the transpose.
c
c     on return
c
c        b       the solution vector  x .
c
c     error condition
c
c        a division by zero will occur if the input factor contains a
c        zero on the diagonal.  technically this indicates singularity
c        but it is often caused by improper arguments or improper
c        setting of lda .  it will not occur if the subroutines are
c        called correctly and if dgbco has set rcond .gt. 0.0
c        or dgbfa has set info .eq. 0 .
c
c     to compute  inverse(a) * c  where  c  is a matrix
c     with  p  columns
c           call dgbco(abd,lda,n,ml,mu,ipvt,rcond,z)
c           if (rcond is too small) go to ...
c           do 10 j = 1, p
c              call dgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0)
c        10 continue
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas daxpy,ddot
c     fortran min0
c
c     internal variables
c
      double precision ddot,t
      integer k,kb,l,la,lb,lm,m,nm1
c
      m = mu + ml + 1
      nm1 = n - 1
      if (job .ne. 0) go to 50
c
c        job = 0 , solve  a * x = b
c        first solve l*y = b
c
         if (ml .eq. 0) go to 30
         if (nm1 .lt. 1) go to 30
            do 20 k = 1, nm1
               lm = min0(ml,n-k)
               l = ipvt(k)
               t = b(l)
               if (l .eq. k) go to 10
                  b(l) = b(k)
                  b(k) = t
   10          continue
               call daxpy(lm,t,abd(m+1,k),1,b(k+1),1)
   20       continue
   30    continue
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k)/abd(m,k)
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = -b(k)
            call daxpy(lm,t,abd(la,k),1,b(lb),1)
   40    continue
      go to 100
   50 continue
c
c        job = nonzero, solve  trans(a) * x = b
c        first solve  trans(u)*y = b
c
         do 60 k = 1, n
            lm = min0(k,m) - 1
            la = m - lm
            lb = k - lm
            t = ddot(lm,abd(la,k),1,b(lb),1)
            b(k) = (b(k) - t)/abd(m,k)
   60    continue
c
c        now solve trans(l)*x = y
c
         if (ml .eq. 0) go to 90
         if (nm1 .lt. 1) go to 90
            do 80 kb = 1, nm1
               k = n - kb
               lm = min0(ml,n-k)
               b(k) = b(k) + ddot(lm,abd(m+1,k),1,b(k+1),1)
               l = ipvt(k)
               if (l .eq. k) go to 70
                  t = b(l)
                  b(l) = b(k)
                  b(k) = t
   70          continue
   80       continue
   90    continue
  100 continue
      return
      end

      SUBROUTINE GETII (II)
      INTEGER II
      II=325
      RETURN
      END

      SUBROUTINE GETKK (KK)
      INTEGER KK
      KK=53
      RETURN
      END

      SUBROUTINE GETEE (EE)
      INTEGER EE
      EE=5
      RETURN
      END

      SUBROUTINE GETELEMENTS (EE,ELEMENTS)
      INTEGER EE
      CHARACTER (LEN=3), dimension(EE) :: ELEMENTS
      ELEMENTS(1) = "O"
      ELEMENTS(2) = "H"
      ELEMENTS(3) = "C"
      ELEMENTS(4) = "N"
      ELEMENTS(5) = "AR"
      RETURN
      END

      SUBROUTINE GETSPECIES (KK,SPECIES)
      INTEGER KK
      CHARACTER (LEN=16), dimension(KK) :: SPECIES
      SPECIES(1) = "H2"
      SPECIES(2) = "H"
      SPECIES(3) = "O"
      SPECIES(4) = "O2"
      SPECIES(5) = "OH"
      SPECIES(6) = "H2O"
      SPECIES(7) = "HO2"
      SPECIES(8) = "H2O2"
      SPECIES(9) = "C"
      SPECIES(10) = "CH"
      SPECIES(11) = "CH2"
      SPECIES(12) = "CH2(S)"
      SPECIES(13) = "CH3"
      SPECIES(14) = "CH4"
      SPECIES(15) = "CO"
      SPECIES(16) = "CO2"
      SPECIES(17) = "HCO"
      SPECIES(18) = "CH2O"
      SPECIES(19) = "CH2OH"
      SPECIES(20) = "CH3O"
      SPECIES(21) = "CH3OH"
      SPECIES(22) = "C2H"
      SPECIES(23) = "C2H2"
      SPECIES(24) = "C2H3"
      SPECIES(25) = "C2H4"
      SPECIES(26) = "C2H5"
      SPECIES(27) = "C2H6"
      SPECIES(28) = "HCCO"
      SPECIES(29) = "CH2CO"
      SPECIES(30) = "HCCOH"
      SPECIES(31) = "N"
      SPECIES(32) = "NH"
      SPECIES(33) = "NH2"
      SPECIES(34) = "NH3"
      SPECIES(35) = "NNH"
      SPECIES(36) = "NO"
      SPECIES(37) = "NO2"
      SPECIES(38) = "N2O"
      SPECIES(39) = "HNO"
      SPECIES(40) = "CN"
      SPECIES(41) = "HCN"
      SPECIES(42) = "H2CN"
      SPECIES(43) = "HCNN"
      SPECIES(44) = "HCNO"
      SPECIES(45) = "HOCN"
      SPECIES(46) = "HNCO"
      SPECIES(47) = "NCO"
      SPECIES(48) = "N2"
      SPECIES(49) = "AR"
      SPECIES(50) = "C3H7"
      SPECIES(51) = "C3H8"
      SPECIES(52) = "CH2CHO"
      SPECIES(53) = "CH3CHO"
      RETURN
      END

      SUBROUTINE GETWT(KK, WT)
      INTEGER KK
      DOUBLE PRECISION, dimension(KK) :: WT
      WT(1) = 2.015880
      WT(2) = 1.007940
      WT(3) = 15.999400
      WT(4) = 31.998800
      WT(5) = 17.007340
      WT(6) = 18.015280
      WT(7) = 33.006740
      WT(8) = 34.014680
      WT(9) = 12.010700
      WT(10) = 13.018640
      WT(11) = 14.026580
      WT(12) = 14.026580
      WT(13) = 15.034520
      WT(14) = 16.042460
      WT(15) = 28.010100
      WT(16) = 44.009500
      WT(17) = 29.018040
      WT(18) = 30.025980
      WT(19) = 31.033920
      WT(20) = 31.033920
      WT(21) = 32.041860
      WT(22) = 25.029340
      WT(23) = 26.037280
      WT(24) = 27.045220
      WT(25) = 28.053160
      WT(26) = 29.061100
      WT(27) = 30.069040
      WT(28) = 41.028740
      WT(29) = 42.036680
      WT(30) = 42.036680
      WT(31) = 14.006740
      WT(32) = 15.014680
      WT(33) = 16.022620
      WT(34) = 17.030560
      WT(35) = 29.021420
      WT(36) = 30.006140
      WT(37) = 46.005540
      WT(38) = 44.012880
      WT(39) = 31.014080
      WT(40) = 26.017440
      WT(41) = 27.025380
      WT(42) = 28.033320
      WT(43) = 41.032120
      WT(44) = 43.024780
      WT(45) = 43.024780
      WT(46) = 43.024780
      WT(47) = 42.016840
      WT(48) = 28.013480
      WT(49) = 39.948000
      WT(50) = 43.087680
      WT(51) = 44.095620
      WT(52) = 43.044620
      WT(53) = 44.052560
      RETURN
      END

      SUBROUTINE GETCOMPOSITION(EE,KK,COMP)
      INTEGER EE, KK
      INTEGER, dimension(EE,KK) :: COMP
      COMP(:,1) = (/0,2,0,0,0/)
      COMP(:,2) = (/0,1,0,0,0/)
      COMP(:,3) = (/1,0,0,0,0/)
      COMP(:,4) = (/2,0,0,0,0/)
      COMP(:,5) = (/1,1,0,0,0/)
      COMP(:,6) = (/1,2,0,0,0/)
      COMP(:,7) = (/2,1,0,0,0/)
      COMP(:,8) = (/2,2,0,0,0/)
      COMP(:,9) = (/0,0,1,0,0/)
      COMP(:,10) = (/0,1,1,0,0/)
      COMP(:,11) = (/0,2,1,0,0/)
      COMP(:,12) = (/0,2,1,0,0/)
      COMP(:,13) = (/0,3,1,0,0/)
      COMP(:,14) = (/0,4,1,0,0/)
      COMP(:,15) = (/1,0,1,0,0/)
      COMP(:,16) = (/2,0,1,0,0/)
      COMP(:,17) = (/1,1,1,0,0/)
      COMP(:,18) = (/1,2,1,0,0/)
      COMP(:,19) = (/1,3,1,0,0/)
      COMP(:,20) = (/1,3,1,0,0/)
      COMP(:,21) = (/1,4,1,0,0/)
      COMP(:,22) = (/0,1,2,0,0/)
      COMP(:,23) = (/0,2,2,0,0/)
      COMP(:,24) = (/0,3,2,0,0/)
      COMP(:,25) = (/0,4,2,0,0/)
      COMP(:,26) = (/0,5,2,0,0/)
      COMP(:,27) = (/0,6,2,0,0/)
      COMP(:,28) = (/1,1,2,0,0/)
      COMP(:,29) = (/1,2,2,0,0/)
      COMP(:,30) = (/1,2,2,0,0/)
      COMP(:,31) = (/0,0,0,1,0/)
      COMP(:,32) = (/0,1,0,1,0/)
      COMP(:,33) = (/0,2,0,1,0/)
      COMP(:,34) = (/0,3,0,1,0/)
      COMP(:,35) = (/0,1,0,2,0/)
      COMP(:,36) = (/1,0,0,1,0/)
      COMP(:,37) = (/2,0,0,1,0/)
      COMP(:,38) = (/1,0,0,2,0/)
      COMP(:,39) = (/1,1,0,1,0/)
      COMP(:,40) = (/0,0,1,1,0/)
      COMP(:,41) = (/0,1,1,1,0/)
      COMP(:,42) = (/0,2,1,1,0/)
      COMP(:,43) = (/0,1,1,2,0/)
      COMP(:,44) = (/1,1,1,1,0/)
      COMP(:,45) = (/1,1,1,1,0/)
      COMP(:,46) = (/1,1,1,1,0/)
      COMP(:,47) = (/1,0,1,1,0/)
      COMP(:,48) = (/0,0,0,2,0/)
      COMP(:,49) = (/0,0,0,0,1/)
      COMP(:,50) = (/0,7,3,0,0/)
      COMP(:,51) = (/0,8,3,0,0/)
      COMP(:,52) = (/1,3,2,0,0/)
      COMP(:,53) = (/1,4,2,0,0/)
      RETURN
      END

      SUBROUTINE GETTHERMODATA(KK,N)
      INTEGER KK
      DOUBLE PRECISION, dimension(14,KK) :: N
      N(:,1) = (/
     * 3.33727920E+00, -4.94024731E-05,  4.99456778E-07,
     *-1.79566394E-10,  2.00255376E-14, -9.50158922E+02,
     *-3.20502331E+00,  2.34433112E+00,  7.98052075E-03,
     *-1.94781510E-05,  2.01572094E-08, -7.37611761E-12,
     *-9.17935173E+02,  6.83010238E-01/)
      N(:,2) = (/
     * 2.50000001E+00, -2.30842973E-11,  1.61561948E-14,
     *-4.73515235E-18,  4.98197357E-22,  2.54736599E+04,
     *-4.46682914E-01,  2.50000000E+00,  7.05332819E-13,
     *-1.99591964E-15,  2.30081632E-18, -9.27732332E-22,
     * 2.54736599E+04, -4.46682853E-01/)
      N(:,3) = (/
     * 2.56942078E+00, -8.59741137E-05,  4.19484589E-08,
     *-1.00177799E-11,  1.22833691E-15,  2.92175791E+04,
     * 4.78433864E+00,  3.16826710E+00, -3.27931884E-03,
     * 6.64306396E-06, -6.12806624E-09,  2.11265971E-12,
     * 2.91222592E+04,  2.05193346E+00/)
      N(:,4) = (/
     * 3.28253784E+00,  1.48308754E-03, -7.57966669E-07,
     * 2.09470555E-10, -2.16717794E-14, -1.08845772E+03,
     * 5.45323129E+00,  3.78245636E+00, -2.99673416E-03,
     * 9.84730201E-06, -9.68129509E-09,  3.24372837E-12,
     *-1.06394356E+03,  3.65767573E+00/)
      N(:,5) = (/
     * 3.09288767E+00,  5.48429716E-04,  1.26505228E-07,
     *-8.79461556E-11,  1.17412376E-14,  3.85865700E+03,
     * 4.47669610E+00,  3.99201543E+00, -2.40131752E-03,
     * 4.61793841E-06, -3.88113333E-09,  1.36411470E-12,
     * 3.61508056E+03, -1.03925458E-01/)
      N(:,6) = (/
     * 3.03399249E+00,  2.17691804E-03, -1.64072518E-07,
     *-9.70419870E-11,  1.68200992E-14, -3.00042971E+04,
     * 4.96677010E+00,  4.19864056E+00, -2.03643410E-03,
     * 6.52040211E-06, -5.48797062E-09,  1.77197817E-12,
     *-3.02937267E+04, -8.49032208E-01/)
      N(:,7) = (/
     * 4.01721090E+00,  2.23982013E-03, -6.33658150E-07,
     * 1.14246370E-10, -1.07908535E-14,  1.11856713E+02,
     * 3.78510215E+00,  4.30179801E+00, -4.74912051E-03,
     * 2.11582891E-05, -2.42763894E-08,  9.29225124E-12,
     * 2.94808040E+02,  3.71666245E+00/)
      N(:,8) = (/
     * 4.16500285E+00,  4.90831694E-03, -1.90139225E-06,
     * 3.71185986E-10, -2.87908305E-14, -1.78617877E+04,
     * 2.91615662E+00,  4.27611269E+00, -5.42822417E-04,
     * 1.67335701E-05, -2.15770813E-08,  8.62454363E-12,
     *-1.77025821E+04,  3.43505074E+00/)
      N(:,9) = (/
     * 2.49266888E+00,  4.79889284E-05, -7.24335020E-08,
     * 3.74291029E-11, -4.87277893E-15,  8.54512953E+04,
     * 4.80150373E+00,  2.55423955E+00, -3.21537724E-04,
     * 7.33792245E-07, -7.32234889E-10,  2.66521446E-13,
     * 8.54438832E+04,  4.53130848E+00/)
      N(:,10) = (/
     * 2.87846473E+00,  9.70913681E-04,  1.44445655E-07,
     *-1.30687849E-10,  1.76079383E-14,  7.10124364E+04,
     * 5.48497999E+00,  3.48981665E+00,  3.23835541E-04,
     *-1.68899065E-06,  3.16217327E-09, -1.40609067E-12,
     * 7.07972934E+04,  2.08401108E+00/)
      N(:,11) = (/
     * 2.87410113E+00,  3.65639292E-03, -1.40894597E-06,
     * 2.60179549E-10, -1.87727567E-14,  4.62636040E+04,
     * 6.17119324E+00,  3.76267867E+00,  9.68872143E-04,
     * 2.79489841E-06, -3.85091153E-09,  1.68741719E-12,
     * 4.60040401E+04,  1.56253185E+00/)
      N(:,12) = (/
     * 2.29203842E+00,  4.65588637E-03, -2.01191947E-06,
     * 4.17906000E-10, -3.39716365E-14,  5.09259997E+04,
     * 8.62650169E+00,  4.19860411E+00, -2.36661419E-03,
     * 8.23296220E-06, -6.68815981E-09,  1.94314737E-12,
     * 5.04968163E+04, -7.69118967E-01/)
      N(:,13) = (/
     * 2.28571772E+00,  7.23990037E-03, -2.98714348E-06,
     * 5.95684644E-10, -4.67154394E-14,  1.67755843E+04,
     * 8.48007179E+00,  3.67359040E+00,  2.01095175E-03,
     * 5.73021856E-06, -6.87117425E-09,  2.54385734E-12,
     * 1.64449988E+04,  1.60456433E+00/)
      N(:,14) = (/
     * 7.48514950E-02,  1.33909467E-02, -5.73285809E-06,
     * 1.22292535E-09, -1.01815230E-13, -9.46834459E+03,
     * 1.84373180E+01,  5.14987613E+00, -1.36709788E-02,
     * 4.91800599E-05, -4.84743026E-08,  1.66693956E-11,
     *-1.02466476E+04, -4.64130376E+00/)
      N(:,15) = (/
     * 2.71518561E+00,  2.06252743E-03, -9.98825771E-07,
     * 2.30053008E-10, -2.03647716E-14, -1.41518724E+04,
     * 7.81868772E+00,  3.57953347E+00, -6.10353680E-04,
     * 1.01681433E-06,  9.07005884E-10, -9.04424499E-13,
     *-1.43440860E+04,  3.50840928E+00/)
      N(:,16) = (/
     * 3.85746029E+00,  4.41437026E-03, -2.21481404E-06,
     * 5.23490188E-10, -4.72084164E-14, -4.87591660E+04,
     * 2.27163806E+00,  2.35677352E+00,  8.98459677E-03,
     *-7.12356269E-06,  2.45919022E-09, -1.43699548E-13,
     *-4.83719697E+04,  9.90105222E+00/)
      N(:,17) = (/
     * 2.77217438E+00,  4.95695526E-03, -2.48445613E-06,
     * 5.89161778E-10, -5.33508711E-14,  4.01191815E+03,
     * 9.79834492E+00,  4.22118584E+00, -3.24392532E-03,
     * 1.37799446E-05, -1.33144093E-08,  4.33768865E-12,
     * 3.83956496E+03,  3.39437243E+00/)
      N(:,18) = (/
     * 1.76069008E+00,  9.20000082E-03, -4.42258813E-06,
     * 1.00641212E-09, -8.83855640E-14, -1.39958323E+04,
     * 1.36563230E+01,  4.79372315E+00, -9.90833369E-03,
     * 3.73220008E-05, -3.79285261E-08,  1.31772652E-11,
     *-1.43089567E+04,  6.02812900E-01/)
      N(:,19) = (/
     * 3.69266569E+00,  8.64576797E-03, -3.75101120E-06,
     * 7.87234636E-10, -6.48554201E-14, -3.24250627E+03,
     * 5.81043215E+00,  3.86388918E+00,  5.59672304E-03,
     * 5.93271791E-06, -1.04532012E-08,  4.36967278E-12,
     *-3.19391367E+03,  5.47302243E+00/)
      N(:,20) = (/
     * 0.03770799E+02,  0.07871497E-01, -0.02656384E-04,
     * 0.03944431E-08, -0.02112616E-12,  0.12783252E+03,
     * 0.02929575E+02,  0.02106204E+02,  0.07216595E-01,
     * 0.05338472E-04, -0.07377636E-07,  0.02075610E-10,
     * 0.09786011E+04,  0.13152177E+02/)
      N(:,21) = (/
     * 1.78970791E+00,  1.40938292E-02, -6.36500835E-06,
     * 1.38171085E-09, -1.17060220E-13, -2.53748747E+04,
     * 1.45023623E+01,  5.71539582E+00, -1.52309129E-02,
     * 6.52441155E-05, -7.10806889E-08,  2.61352698E-11,
     *-2.56427656E+04, -1.50409823E+00/)
      N(:,22) = (/
     * 3.16780652E+00,  4.75221902E-03, -1.83787077E-06,
     * 3.04190252E-10, -1.77232770E-14,  6.71210650E+04,
     * 6.63589475E+00,  2.88965733E+00,  1.34099611E-02,
     *-2.84769501E-05,  2.94791045E-08, -1.09331511E-11,
     * 6.68393932E+04,  6.22296438E+00/)
      N(:,23) = (/
     * 4.14756964E+00,  5.96166664E-03, -2.37294852E-06,
     * 4.67412171E-10, -3.61235213E-14,  2.59359992E+04,
     *-1.23028121E+00,  8.08681094E-01,  2.33615629E-02,
     *-3.55171815E-05,  2.80152437E-08, -8.50072974E-12,
     * 2.64289807E+04,  1.39397051E+01/)
      N(:,24) = (/
     * 3.01672400E+00,  1.03302292E-02, -4.68082349E-06,
     * 1.01763288E-09, -8.62607041E-14,  3.46128739E+04,
     * 7.78732378E+00,  3.21246645E+00,  1.51479162E-03,
     * 2.59209412E-05, -3.57657847E-08,  1.47150873E-11,
     * 3.48598468E+04,  8.51054025E+00/)
      N(:,25) = (/
     * 2.03611116E+00,  1.46454151E-02, -6.71077915E-06,
     * 1.47222923E-09, -1.25706061E-13,  4.93988614E+03,
     * 1.03053693E+01,  3.95920148E+00, -7.57052247E-03,
     * 5.70990292E-05, -6.91588753E-08,  2.69884373E-11,
     * 5.08977593E+03,  4.09733096E+00/)
      N(:,26) = (/
     * 1.95465642E+00,  1.73972722E-02, -7.98206668E-06,
     * 1.75217689E-09, -1.49641576E-13,  1.28575200E+04,
     * 1.34624343E+01,  4.30646568E+00, -4.18658892E-03,
     * 4.97142807E-05, -5.99126606E-08,  2.30509004E-11,
     * 1.28416265E+04,  4.70720924E+00/)
      N(:,27) = (/
     * 1.07188150E+00,  2.16852677E-02, -1.00256067E-05,
     * 2.21412001E-09, -1.90002890E-13, -1.14263932E+04,
     * 1.51156107E+01,  4.29142492E+00, -5.50154270E-03,
     * 5.99438288E-05, -7.08466285E-08,  2.68685771E-11,
     *-1.15222055E+04,  2.66682316E+00/)
      N(:,28) = (/
     * 0.56282058E+01,  0.40853401E-02, -0.15934547E-05,
     * 0.28626052E-09, -0.19407832E-13,  0.19327215E+05,
     *-0.39302595E+01,  0.22517214E+01,  0.17655021E-01,
     *-0.23729101E-04,  0.17275759E-07, -0.50664811E-11,
     * 0.20059449E+05,  0.12490417E+02/)
      N(:,29) = (/
     * 4.51129732E+00,  9.00359745E-03, -4.16939635E-06,
     * 9.23345882E-10, -7.94838201E-14, -7.55105311E+03,
     * 6.32247205E-01,  2.13583630E+00,  1.81188721E-02,
     *-1.73947474E-05,  9.34397568E-09, -2.01457615E-12,
     *-7.04291804E+03,  1.22156480E+01/)
      N(:,30) = (/
     * 0.59238291E+01,  0.67923600E-02, -0.25658564E-05,
     * 0.44987841E-09, -0.29940101E-13,  0.72646260E+04,
     *-0.76017742E+01,  0.12423733E+01,  0.31072201E-01,
     *-0.50866864E-04,  0.43137131E-07, -0.14014594E-10,
     * 0.80316143E+04,  0.13874319E+02/)
      N(:,31) = (/
     * 0.24159429E+01,  0.17489065E-03, -0.11902369E-06,
     * 0.30226245E-10, -0.20360982E-14,  0.56133773E+05,
     * 0.46496096E+01,  0.25000000E+01,  0.00000000E+00,
     * 0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
     * 0.56104637E+05,  0.41939087E+01/)
      N(:,32) = (/
     * 0.27836928E+01,  0.13298430E-02, -0.42478047E-06,
     * 0.78348501E-10, -0.55044470E-14,  0.42120848E+05,
     * 0.57407799E+01,  0.34929085E+01,  0.31179198E-03,
     *-0.14890484E-05,  0.24816442E-08, -0.10356967E-11,
     * 0.41880629E+05,  0.18483278E+01/)
      N(:,33) = (/
     * 0.28347421E+01,  0.32073082E-02, -0.93390804E-06,
     * 0.13702953E-09, -0.79206144E-14,  0.22171957E+05,
     * 0.65204163E+01,  0.42040029E+01, -0.21061385E-02,
     * 0.71068348E-05, -0.56115197E-08,  0.16440717E-11,
     * 0.21885910E+05, -0.14184248E+00/)
      N(:,34) = (/
     * 0.26344521E+01,  0.56662560E-02, -0.17278676E-05,
     * 0.23867161E-09, -0.12578786E-13, -0.65446958E+04,
     * 0.65662928E+01,  0.42860274E+01, -0.46605230E-02,
     * 0.21718513E-04, -0.22808887E-07,  0.82638046E-11,
     *-0.67417285E+04, -0.62537277E+00/)
      N(:,35) = (/
     * 0.37667544E+01,  0.28915082E-02, -0.10416620E-05,
     * 0.16842594E-09, -0.10091896E-13,  0.28650697E+05,
     * 0.44705067E+01,  0.43446927E+01, -0.48497072E-02,
     * 0.20059459E-04, -0.21726464E-07,  0.79469539E-11,
     * 0.28791973E+05,  0.29779410E+01/)
      N(:,36) = (/
     * 0.32606056E+01,  0.11911043E-02, -0.42917048E-06,
     * 0.69457669E-10, -0.40336099E-14,  0.99209746E+04,
     * 0.63693027E+01,  0.42184763E+01, -0.46389760E-02,
     * 0.11041022E-04, -0.93361354E-08,  0.28035770E-11,
     * 0.98446230E+04,  0.22808464E+01/)
      N(:,37) = (/
     * 0.48847542E+01,  0.21723956E-02, -0.82806906E-06,
     * 0.15747510E-09, -0.10510895E-13,  0.23164983E+04,
     *-0.11741695E+00,  0.39440312E+01, -0.15854290E-02,
     * 0.16657812E-04, -0.20475426E-07,  0.78350564E-11,
     * 0.28966179E+04,  0.63119917E+01/)
      N(:,38) = (/
     * 0.48230729E+01,  0.26270251E-02, -0.95850874E-06,
     * 0.16000712E-09, -0.97752303E-14,  0.80734048E+04,
     *-0.22017207E+01,  0.22571502E+01,  0.11304728E-01,
     *-0.13671319E-04,  0.96819806E-08, -0.29307182E-11,
     * 0.87417744E+04,  0.10757992E+02/)
      N(:,39) = (/
     * 0.29792509E+01,  0.34944059E-02, -0.78549778E-06,
     * 0.57479594E-10, -0.19335916E-15,  0.11750582E+05,
     * 0.86063728E+01,  0.45334916E+01, -0.56696171E-02,
     * 0.18473207E-04, -0.17137094E-07,  0.55454573E-11,
     * 0.11548297E+05,  0.17498417E+01/)
      N(:,40) = (/
     * 0.37459805E+01,  0.43450775E-04,  0.29705984E-06,
     *-0.68651806E-10,  0.44134173E-14,  0.51536188E+05,
     * 0.27867601E+01,  0.36129351E+01, -0.95551327E-03,
     * 0.21442977E-05, -0.31516323E-09, -0.46430356E-12,
     * 0.51708340E+05,  0.39804995E+01/)
      N(:,41) = (/
     * 0.38022392E+01,  0.31464228E-02, -0.10632185E-05,
     * 0.16619757E-09, -0.97997570E-14,  0.14407292E+05,
     * 0.15754601E+01,  0.22589886E+01,  0.10051170E-01,
     *-0.13351763E-04,  0.10092349E-07, -0.30089028E-11,
     * 0.14712633E+05,  0.89164419E+01/)
      N(:,42) = (/
     * 0.52097030E+01,  0.29692911E-02, -0.28555891E-06,
     *-0.16355500E-09,  0.30432589E-13,  0.27677109E+05,
     *-0.44444780E+01,  0.28516610E+01,  0.56952331E-02,
     * 0.10711400E-05, -0.16226120E-08, -0.23511081E-12,
     * 0.28637820E+05,  0.89927511E+01/)
      N(:,43) = (/
     * 0.58946362E+01,  0.39895959E-02, -0.15982380E-05,
     * 0.29249395E-09, -0.20094686E-13,  0.53452941E+05,
     *-0.51030502E+01,  0.25243194E+01,  0.15960619E-01,
     *-0.18816354E-04,  0.12125540E-07, -0.32357378E-11,
     * 0.54261984E+05,  0.11675870E+02/)
      N(:,44) = (/
     * 6.59860456E+00,  3.02778626E-03, -1.07704346E-06,
     * 1.71666528E-10, -1.01439391E-14,  1.79661339E+04,
     *-1.03306599E+01,  2.64727989E+00,  1.27505342E-02,
     *-1.04794236E-05,  4.41432836E-09, -7.57521466E-13,
     * 1.92990252E+04,  1.07332972E+01/)
      N(:,45) = (/
     * 5.89784885E+00,  3.16789393E-03, -1.11801064E-06,
     * 1.77243144E-10, -1.04339177E-14, -3.70653331E+03,
     *-6.18167825E+00,  3.78604952E+00,  6.88667922E-03,
     *-3.21487864E-06,  5.17195767E-10,  1.19360788E-14,
     *-2.82698400E+03,  5.63292162E+00/)
      N(:,46) = (/
     * 6.22395134E+00,  3.17864004E-03, -1.09378755E-06,
     * 1.70735163E-10, -9.95021955E-15, -1.66599344E+04,
     *-8.38224741E+00,  3.63096317E+00,  7.30282357E-03,
     *-2.28050003E-06, -6.61271298E-10,  3.62235752E-13,
     *-1.55873636E+04,  6.19457727E+00/)
      N(:,47) = (/
     * 0.51521845E+01,  0.23051761E-02, -0.88033153E-06,
     * 0.14789098E-09, -0.90977996E-14,  0.14004123E+05,
     *-0.25442660E+01,  0.28269308E+01,  0.88051688E-02,
     *-0.83866134E-05,  0.48016964E-08, -0.13313595E-11,
     * 0.14682477E+05,  0.95504646E+01/)
      N(:,48) = (/
     * 0.02926640E+02,  0.14879768E-02, -0.05684760E-05,
     * 0.10097038E-09, -0.06753351E-13, -0.09227977E+04,
     * 0.05980528E+02,  0.03298677E+02,  0.14082404E-02,
     *-0.03963222E-04,  0.05641515E-07, -0.02444854E-10,
     *-0.10208999E+04,  0.03950372E+02/)
      N(:,49) = (/
     * 0.02500000E+02,  0.00000000E+00,  0.00000000E+00,
     * 0.00000000E+00,  0.00000000E+00, -0.07453750E+04,
     * 0.04366000E+02,  0.02500000E+02,  0.00000000E+00,
     * 0.00000000E+00,  0.00000000E+00,  0.00000000E+00,
     *-0.07453750E+04,  0.04366000E+02/)
      N(:,50) = (/
     * 0.77026987E+01,  0.16044203E-01, -0.52833220E-05,
     * 0.76298590E-09, -0.39392284E-13,  0.82984336E+04,
     *-0.15480180E+02,  0.10515518E+01,  0.25991980E-01,
     * 0.23800540E-05, -0.19609569E-07,  0.93732470E-11,
     * 0.10631863E+05,  0.21122559E+02/)
      N(:,51) = (/
     * 0.75341368E+01,  0.18872239E-01, -0.62718491E-05,
     * 0.91475649E-09, -0.47838069E-13, -0.16467516E+05,
     *-0.17892349E+02,  0.93355381E+00,  0.26424579E-01,
     * 0.61059727E-05, -0.21977499E-07,  0.95149253E-11,
     *-0.13958520E+05,  0.19201691E+02/)
      N(:,52) = (/
     * 0.05975670E+02,  0.08130591E-01, -0.02743624E-04,
     * 0.04070304E-08, -0.02176017E-12,  0.04903218E+04,
     *-0.05045251E+02,  0.03409062E+02,  0.10738574E-01,
     * 0.01891492E-04, -0.07158583E-07,  0.02867385E-10,
     * 0.15214766E+04,  0.09558290E+02/)
      N(:,53) = (/
     * 0.54041108E+01,  0.11723059E-01, -0.42263137E-05,
     * 0.68372451E-09, -0.40984863E-13, -0.22593122E+05,
     *-0.34807917E+01,  0.47294595E+01, -0.31932858E-02,
     * 0.47534921E-04, -0.57458611E-07,  0.21931112E-10,
     *-0.21572878E+05,  0.41030159E+01/)
      RETURN
      END

      SUBROUTINE GETCOMMONT(KK, COMMONT)
      INTEGER KK
      DOUBLE PRECISION, dimension(KK) :: COMMONT
      COMMONT(1) =   1000.0
      COMMONT(2) =   1000.0
      COMMONT(3) =   1000.0
      COMMONT(4) =   1000.0
      COMMONT(5) =   1000.0
      COMMONT(6) =   1000.0
      COMMONT(7) =   1000.0
      COMMONT(8) =   1000.0
      COMMONT(9) =   1000.0
      COMMONT(10) =   1000.0
      COMMONT(11) =   1000.0
      COMMONT(12) =   1000.0
      COMMONT(13) =   1000.0
      COMMONT(14) =   1000.0
      COMMONT(15) =   1000.0
      COMMONT(16) =   1000.0
      COMMONT(17) =   1000.0
      COMMONT(18) =   1000.0
      COMMONT(19) =   1000.0
      COMMONT(20) =   1000.0
      COMMONT(21) =   1000.0
      COMMONT(22) =   1000.0
      COMMONT(23) =   1000.0
      COMMONT(24) =   1000.0
      COMMONT(25) =   1000.0
      COMMONT(26) =   1000.0
      COMMONT(27) =   1000.0
      COMMONT(28) =   1000.0
      COMMONT(29) =   1000.0
      COMMONT(30) =   1000.0
      COMMONT(31) =   1000.0
      COMMONT(32) =   1000.0
      COMMONT(33) =   1000.0
      COMMONT(34) =   1000.0
      COMMONT(35) =   1000.0
      COMMONT(36) =   1000.0
      COMMONT(37) =   1000.0
      COMMONT(38) =   1000.0
      COMMONT(39) =   1000.0
      COMMONT(40) =   1000.0
      COMMONT(41) =   1000.0
      COMMONT(42) =   1000.0
      COMMONT(43) =   1000.0
      COMMONT(44) =   1382.0
      COMMONT(45) =   1368.0
      COMMONT(46) =   1478.0
      COMMONT(47) =   1000.0
      COMMONT(48) =   1000.0
      COMMONT(49) =   1000.0
      COMMONT(50) =   1000.0
      COMMONT(51) =   1000.0
      COMMONT(52) =   1000.0
      COMMONT(53) =   1000.0
      RETURN
      END

      SUBROUTINE GETABE (II, ABE)
      INTEGER II
      DOUBLE PRECISION, dimension(3,II) :: ABE
      ABE(:,1) = (/3.932627e+01,-1.000,.00/)
      ABE(:,2) = (/4.075338e+01,-1.000,.00/)
      ABE(:,3) = (/1.056359e+01,2.700,6260.00/)
      ABE(:,4) = (/3.062675e+01,.000,.00/)
      ABE(:,5) = (/1.608039e+01,2.000,4000.00/)
      ABE(:,6) = (/3.167407e+01,.000,.00/)
      ABE(:,7) = (/3.201305e+01,.000,.00/)
      ABE(:,8) = (/3.033907e+01,.000,.00/)
      ABE(:,9) = (/3.033907e+01,.000,.00/)
      ABE(:,10) = (/3.155497e+01,.000,.00/)
      ABE(:,11) = (/2.074307e+01,1.500,8600.00/)
      ABE(:,12) = (/2.361364e+01,.000,2385.00/)
      ABE(:,13) = (/3.103222e+01,.000,.00/)
      ABE(:,14) = (/3.103222e+01,.000,.00/)
      ABE(:,15) = (/3.129458e+01,.000,3540.00/)
      ABE(:,16) = (/2.993361e+01,.000,.00/)
      ABE(:,17) = (/2.993361e+01,.000,.00/)
      ABE(:,18) = (/1.286876e+01,2.500,3100.00/)
      ABE(:,19) = (/1.177529e+01,2.500,5000.00/)
      ABE(:,20) = (/3.154304e+01,.000,.00/)
      ABE(:,21) = (/1.641820e+01,2.000,1900.00/)
      ABE(:,22) = (/4.527517e+01,-1.410,28950.00/)
      ABE(:,23) = (/1.575281e+01,2.000,1900.00/)
      ABE(:,24) = (/3.103222e+01,.000,.00/)
      ABE(:,25) = (/1.634124e+01,1.830,220.00/)
      ABE(:,26) = (/3.074008e+01,.000,.00/)
      ABE(:,27) = (/1.831310e+01,1.920,5690.00/)
      ABE(:,28) = (/3.223619e+01,.000,.00/)
      ABE(:,29) = (/2.993361e+01,.000,8000.00/)
      ABE(:,30) = (/2.819064e+01,.000,1350.00/)
      ABE(:,31) = (/2.854731e+01,.000,47800.00/)
      ABE(:,32) = (/3.223619e+01,.000,40000.00/)
      ABE(:,33) = (/4.247615e+01,-.860,.00/)
      ABE(:,34) = (/4.448148e+01,-1.240,.00/)
      ABE(:,35) = (/4.386779e+01,-.760,.00/)
      ABE(:,36) = (/4.470463e+01,-1.240,.00/)
      ABE(:,37) = (/4.108986e+01,-.800,.00/)
      ABE(:,38) = (/3.781592e+01,-.6707,17041.00/)
      ABE(:,39) = (/4.144653e+01,-1.000,.00/)
      ABE(:,40) = (/3.903859e+01,-.600,.00/)
      ABE(:,41) = (/4.554088e+01,-1.250,.00/)
      ABE(:,42) = (/4.775645e+01,-2.000,.00/)
      ABE(:,43) = (/5.144533e+01,-2.000,.00/)
      ABE(:,44) = (/2.900979e+01,.000,671.00/)
      ABE(:,45) = (/3.143323e+01,.000,1068.00/)
      ABE(:,46) = (/3.206184e+01,.000,635.00/)
      ABE(:,47) = (/1.630872e+01,2.000,5200.00/)
      ABE(:,48) = (/2.993361e+01,.000,3600.00/)
      ABE(:,49) = (/3.273697e+01,.000,.00/)
      ABE(:,50) = (/3.402795e+01,.000,.00/)
      ABE(:,51) = (/3.103222e+01,.000,.00/)
      ABE(:,52) = (/3.717067e+01,-.534,536.00/)
      ABE(:,53) = (/2.030775e+01,1.620,10840.00/)
      ABE(:,54) = (/2.771720e+01,.480,-260.00/)
      ABE(:,55) = (/3.192695e+01,.000,.00/)
      ABE(:,56) = (/2.701483e+01,.454,3600.00/)
      ABE(:,57) = (/2.701483e+01,.454,2600.00/)
      ABE(:,58) = (/1.786555e+01,1.900,2742.00/)
      ABE(:,59) = (/2.768456e+01,.500,86.00/)
      ABE(:,60) = (/3.062675e+01,.000,.00/)
      ABE(:,61) = (/2.582921e+01,.650,-284.00/)
      ABE(:,62) = (/3.112145e+01,-.090,610.00/)
      ABE(:,63) = (/2.851891e+01,.515,50.00/)
      ABE(:,64) = (/1.754120e+01,1.630,1924.00/)
      ABE(:,65) = (/3.062675e+01,.000,.00/)
      ABE(:,66) = (/2.803649e+01,.500,-110.00/)
      ABE(:,67) = (/3.319937e+01,-.230,1070.00/)
      ABE(:,68) = (/1.664872e+01,2.100,4870.00/)
      ABE(:,69) = (/1.525060e+01,2.100,4870.00/)
      ABE(:,70) = (/3.914395e+01,-1.000,.00/)
      ABE(:,71) = (/2.935379e+01,.000,2400.00/)
      ABE(:,72) = (/2.943603e+01,.270,280.00/)
      ABE(:,73) = (/3.103222e+01,.000,.00/)
      ABE(:,74) = (/2.701483e+01,.454,1820.00/)
      ABE(:,75) = (/1.409692e+01,2.530,12240.00/)
      ABE(:,76) = (/4.079453e+01,-.990,1580.00/)
      ABE(:,77) = (/2.832417e+01,.000,.00/)
      ABE(:,78) = (/1.856044e+01,1.900,7530.00/)
      ABE(:,79) = (/3.223619e+01,.000,.00/)
      ABE(:,80) = (/3.154304e+01,.000,8000.00/)
      ABE(:,81) = (/3.005582e+01,.000,3428.00/)
      ABE(:,82) = (/2.993361e+01,.000,.00/)
      ABE(:,83) = (/1.757671e+01,1.500,79600.00/)
      ABE(:,84) = (/1.919079e+01,1.510,3430.00/)
      ABE(:,85) = (/3.193509e+01,-.370,.00/)
      ABE(:,86) = (/1.048291e+01,2.400,-2110.00/)
      ABE(:,87) = (/3.030517e+01,.000,-500.00/)
      ABE(:,88) = (/2.832417e+01,.000,427.00/)
      ABE(:,89) = (/4.197716e+01,.000,29410.00/)
      ABE(:,90) = (/3.154304e+01,.000,.00/)
      ABE(:,91) = (/3.103222e+01,.000,.00/)
      ABE(:,92) = (/3.062675e+01,.000,.00/)
      ABE(:,93) = (/1.624031e+01,2.000,3000.00/)
      ABE(:,94) = (/3.103222e+01,.000,.00/)
      ABE(:,95) = (/4.247257e+01,-1.430,1330.00/)
      ABE(:,96) = (/1.784086e+01,1.600,5420.00/)
      ABE(:,97) = (/4.100648e+01,-1.340,1417.00/)
      ABE(:,98) = (/1.842068e+01,1.600,3120.00/)
      ABE(:,99) = (/1.767834e+01,1.228,70.00/)
      ABE(:,100) = (/3.154304e+01,.000,.00/)
      ABE(:,101) = (/2.195583e+01,1.180,-447.00/)
      ABE(:,102) = (/2.924046e+01,.000,.00/)
      ABE(:,103) = (/2.924046e+01,.000,.00/)
      ABE(:,104) = (/1.418015e+01,2.000,-840.00/)
      ABE(:,105) = (/1.565606e+01,2.000,1500.00/)
      ABE(:,106) = (/3.062675e+01,.000,.00/)
      ABE(:,107) = (/-8.431015e+00,4.500,-1000.00/)
      ABE(:,108) = (/1.313033e+01,2.300,13500.00/)
      ABE(:,109) = (/1.733301e+01,2.000,14000.00/)
      ABE(:,110) = (/-7.635494e+00,4.000,-2000.00/)
      ABE(:,111) = (/2.924046e+01,.000,.00/)
      ABE(:,112) = (/1.509644e+01,2.000,2500.00/)
      ABE(:,113) = (/1.507964e+01,2.120,870.00/)
      ABE(:,114) = (/2.964592e+01,.000,2000.00/)
      ABE(:,115) = (/2.559080e+01,.000,-1630.00/)
      ABE(:,116) = (/3.367128e+01,.000,12000.00/)
      ABE(:,117) = (/3.062675e+01,.000,.00/)
      ABE(:,118) = (/2.763102e+01,.000,.00/)
      ABE(:,119) = (/3.126333e+01,.000,.00/)
      ABE(:,120) = (/3.264166e+01,.000,23600.00/)
      ABE(:,121) = (/1.553828e+01,2.000,12000.00/)
      ABE(:,122) = (/3.169146e+01,.000,576.00/)
      ABE(:,123) = (/3.154304e+01,.000,.00/)
      ABE(:,124) = (/3.154304e+01,.000,.00/)
      ABE(:,125) = (/3.183721e+01,.000,.00/)
      ABE(:,126) = (/3.231315e+01,.000,3110.00/)
      ABE(:,127) = (/2.937324e+01,.000,-755.00/)
      ABE(:,128) = (/3.131990e+01,.000,.00/)
      ABE(:,129) = (/3.103222e+01,.000,.00/)
      ABE(:,130) = (/3.172537e+01,.000,.00/)
      ABE(:,131) = (/3.154304e+01,.000,.00/)
      ABE(:,132) = (/3.287805e+01,.000,15792.00/)
      ABE(:,133) = (/3.218068e+01,.000,-515.00/)
      ABE(:,134) = (/3.154304e+01,.000,.00/)
      ABE(:,135) = (/2.924046e+01,.000,1500.00/)
      ABE(:,136) = (/1.312236e+01,2.000,7230.00/)
      ABE(:,137) = (/3.500878e+01,.000,11944.00/)
      ABE(:,138) = (/3.131990e+01,.000,.00/)
      ABE(:,139) = (/1.471567e+01,2.000,8270.00/)
      ABE(:,140) = (/2.742030e+01,.500,4510.00/)
      ABE(:,141) = (/3.103222e+01,.000,.00/)
      ABE(:,142) = (/3.033907e+01,.000,600.00/)
      ABE(:,143) = (/2.982825e+01,.000,600.00/)
      ABE(:,144) = (/3.096323e+01,.000,.00/)
      ABE(:,145) = (/3.011593e+01,.000,.00/)
      ABE(:,146) = (/3.187952e+01,.000,.00/)
      ABE(:,147) = (/4.071672e+01,-1.160,1145.00/)
      ABE(:,148) = (/3.103222e+01,.000,.00/)
      ABE(:,149) = (/3.011593e+01,.000,-570.00/)
      ABE(:,150) = (/3.040361e+01,.000,-570.00/)
      ABE(:,151) = (/2.982825e+01,.000,.00/)
      ABE(:,152) = (/2.957693e+01,.000,.00/)
      ABE(:,153) = (/3.027008e+01,.000,.00/)
      ABE(:,154) = (/3.131990e+01,.000,-550.00/)
      ABE(:,155) = (/3.120337e+01,.000,30480.00/)
      ABE(:,156) = (/2.846827e+01,.000,20315.00/)
      ABE(:,157) = (/1.010643e+01,2.470,5180.00/)
      ABE(:,158) = (/3.875386e+01,-1.180,654.00/)
      ABE(:,159) = (/2.955381e+01,.100,10600.00/)
      ABE(:,160) = (/3.090741e+01,.000,.00/)
      ABE(:,161) = (/8.107720e+00,2.810,5860.00/)
      ABE(:,162) = (/1.721671e+01,1.500,9940.00/)
      ABE(:,163) = (/1.611810e+01,1.500,9940.00/)
      ABE(:,164) = (/1.233271e+01,2.000,9200.00/)
      ABE(:,165) = (/1.563034e+01,1.740,10450.00/)
      ABE(:,166) = (/4.185200e+01,-1.000,17000.00/)
      ABE(:,167) = (/3.976989e+01,-1.000,17000.00/)
      ABE(:,168) = (/3.023000e+01,.000,400.00/)
      ABE(:,169) = (/3.052139e+01,.000,900.00/)
      ABE(:,170) = (/-2.847965e+01,7.600,-3530.00/)
      ABE(:,171) = (/2.993361e+01,.000,-755.00/)
      ABE(:,172) = (/2.476280e+01,0.900,1993.00/)
      ABE(:,173) = (/3.836306e+01,-1.390,1015.00/)
      ABE(:,174) = (/2.971046e+01,.440,86770.00/)
      ABE(:,175) = (/2.745667e+01,.000,3875.00/)
      ABE(:,176) = (/2.879417e+01,.000,854.00/)
      ABE(:,177) = (/2.993361e+01,.000,.00/)
      ABE(:,178) = (/3.092686e+01,.000,355.00/)
      ABE(:,179) = (/2.292049e+01,1.000,6500.00/)
      ABE(:,180) = (/3.114555e+01,.000,385.00/)
      ABE(:,181) = (/2.796749e+01,.000,10810.00/)
      ABE(:,182) = (/3.099832e+01,.000,23150.00/)
      ABE(:,183) = (/3.358945e+01,.000,18880.00/)
      ABE(:,184) = (/2.832417e+01,.000,21060.00/)
      ABE(:,185) = (/2.509398e+01,.000,56020.00/)
      ABE(:,186) = (/2.837771e+01,.000,-480.00/)
      ABE(:,187) = (/4.610997e+01,-1.410,.00/)
      ABE(:,188) = (/2.899200e+01,.000,-240.00/)
      ABE(:,189) = (/3.251382e+01,.000,360.00/)
      ABE(:,190) = (/3.131990e+01,.000,.00/)
      ABE(:,191) = (/3.109676e+01,.000,330.00/)
      ABE(:,192) = (/3.062675e+01,.000,.00/)
      ABE(:,193) = (/2.141641e+01,1.200,.00/)
      ABE(:,194) = (/1.304115e+01,2.000,6500.00/)
      ABE(:,195) = (/1.406237e+01,1.500,100.00/)
      ABE(:,196) = (/3.033907e+01,.000,.00/)
      ABE(:,197) = (/3.062675e+01,.000,13850.00/)
      ABE(:,198) = (/3.070371e+01,-.230,.00/)
      ABE(:,199) = (/3.353092e+01,-.450,.00/)
      ABE(:,200) = (/2.872963e+01,.000,.00/)
      ABE(:,201) = (/3.129458e+01,.000,.00/)
      ABE(:,202) = (/3.131990e+01,.000,3650.00/)
      ABE(:,203) = (/1.831532e+01,1.500,-460.00/)
      ABE(:,204) = (/1.961460e+01,.000,.00/)
      ABE(:,205) = (/3.249856e+01,-.110,4980.00/)
      ABE(:,206) = (/2.924046e+01,.000,.00/)
      ABE(:,207) = (/3.084990e+01,.000,.00/)
      ABE(:,208) = (/3.187952e+01,.000,.00/)
      ABE(:,209) = (/3.154304e+01,.000,.00/)
      ABE(:,210) = (/3.062675e+01,.000,.00/)
      ABE(:,211) = (/3.084990e+01,.000,.00/)
      ABE(:,212) = (/4.524874e+01,-1.320,740.00/)
      ABE(:,213) = (/3.084990e+01,.000,.00/)
      ABE(:,214) = (/2.752566e+01,.720,660.00/)
      ABE(:,215) = (/1.638046e+01,1.900,-950.00/)
      ABE(:,216) = (/2.993361e+01,.000,13000.00/)
      ABE(:,217) = (/3.197483e+01,.000,.00/)
      ABE(:,218) = (/3.131990e+01,.000,.00/)
      ABE(:,219) = (/2.971046e+01,.000,7460.00/)
      ABE(:,220) = (/2.944585e+01,.000,-440.00/)
      ABE(:,221) = (/1.259473e+01,2.450,2240.00/)
      ABE(:,222) = (/3.078802e+01,.000,.00/)
      ABE(:,223) = (/3.162001e+01,.000,.00/)
      ABE(:,224) = (/2.854731e+01,.000,.00/)
      ABE(:,225) = (/3.062675e+01,.000,.00/)
      ABE(:,226) = (/2.832417e+01,.000,20000.00/)
      ABE(:,227) = (/3.336759e+01,.000,54050.00/)
      ABE(:,228) = (/3.978580e+01,-1.520,740.00/)
      ABE(:,229) = (/4.278153e+01,-2.000,800.00/)
      ABE(:,230) = (/6.681419e+01,-3.300,126600.00/)
      ABE(:,231) = (/9.918376e+00,2.640,4980.00/)
      ABE(:,232) = (/8.531096e+00,2.640,4980.00/)
      ABE(:,233) = (/2.208680e+01,1.580,26600.00/)
      ABE(:,234) = (/1.391082e+01,2.030,13370.00/)
      ABE(:,235) = (/8.389360e+00,2.260,6400.00/)
      ABE(:,236) = (/5.075174e+00,2.560,9000.00/)
      ABE(:,237) = (/3.112753e+01,.000,.00/)
      ABE(:,238) = (/3.172537e+01,.000,400.00/)
      ABE(:,239) = (/3.177416e+01,.000,46020.00/)
      ABE(:,240) = (/2.186110e+01,0.880,20130.00/)
      ABE(:,241) = (/2.876242e+01,.150,.00/)
      ABE(:,242) = (/2.993361e+01,.000,74000.00/)
      ABE(:,243) = (/2.532844e+01,.000,65000.00/)
      ABE(:,244) = (/3.057546e+01,.000,.00/)
      ABE(:,245) = (/3.099832e+01,.000,.00/)
      ABE(:,246) = (/3.134459e+01,.000,.00/)
      ABE(:,247) = (/3.041603e+01,.000,.00/)
      ABE(:,248) = (/3.083377e+01,.000,.00/)
      ABE(:,249) = (/4.027535e+01,-1.380,1270.00/)
      ABE(:,250) = (/3.330090e+01,-.690,760.00/)
      ABE(:,251) = (/3.126861e+01,-.360,580.00/)
      ABE(:,252) = (/4.027535e+01,-1.380,1270.00/)
      ABE(:,253) = (/3.330090e+01,-.690,760.00/)
      ABE(:,254) = (/3.126861e+01,-.360,580.00/)
      ABE(:,255) = (/3.219537e+01,.000,28800.00/)
      ABE(:,256) = (/2.763102e+01,.000,21750.00/)
      ABE(:,257) = (/3.072206e+01,.000,.00/)
      ABE(:,258) = (/2.832417e+01,.000,.00/)
      ABE(:,259) = (/3.011593e+01,.000,.00/)
      ABE(:,260) = (/3.011593e+01,.000,.00/)
      ABE(:,261) = (/3.223619e+01,.000,.00/)
      ABE(:,262) = (/1.840048e+01,1.410,8500.00/)
      ABE(:,263) = (/1.882615e+01,1.570,44000.00/)
      ABE(:,264) = (/1.460397e+01,2.110,11400.00/)
      ABE(:,265) = (/1.692903e+01,1.700,3800.00/)
      ABE(:,266) = (/1.156172e+01,2.500,13300.00/)
      ABE(:,267) = (/1.731202e+01,1.500,3600.00/)
      ABE(:,268) = (/1.500943e+01,1.500,3600.00/)
      ABE(:,269) = (/3.700688e+01,.000,84720.00/)
      ABE(:,270) = (/3.528071e+01,-.690,2850.00/)
      ABE(:,271) = (/2.632169e+01,.180,2120.00/)
      ABE(:,272) = (/3.276682e+01,-.750,2890.00/)
      ABE(:,273) = (/1.681124e+01,2.000,2000.00/)
      ABE(:,274) = (/2.982825e+01,.000,.00/)
      ABE(:,275) = (/3.404448e+01,-.310,290.00/)
      ABE(:,276) = (/2.893935e+01,.150,-90.00/)
      ABE(:,277) = (/1.319932e+01,2.400,9915.00/)
      ABE(:,278) = (/1.772753e+01,1.600,955.00/)
      ABE(:,279) = (/1.605622e+01,1.940,6460.00/)
      ABE(:,280) = (/2.993361e+01,.000,14350.00/)
      ABE(:,281) = (/3.635685e+01,-0.752,345.00/)
      ABE(:,282) = (/2.880968e+01,.000,-705.00/)
      ABE(:,283) = (/2.872963e+01,.000,11300.00/)
      ABE(:,284) = (/3.114852e+01,.000,.00/)
      ABE(:,285) = (/1.571762e+01,1.830,220.00/)
      ABE(:,286) = (/3.232786e+01,.000,.00/)
      ABE(:,287) = (/3.614821e+01,.000,17330.00/)
      ABE(:,288) = (/2.280271e+01,.500,-1755.00/)
      ABE(:,289) = (/2.830905e+01,.430,-370.00/)
      ABE(:,290) = (/2.938888e+01,.000,1500.00/)
      ABE(:,291) = (/2.850649e+01,.000,1500.00/)
      ABE(:,292) = (/3.292934e+01,.000,10989.00/)
      ABE(:,293) = (/2.494571e+01,.250,-935.00/)
      ABE(:,294) = (/2.643700e+01,.290,11.00/)
      ABE(:,295) = (/1.410594e+01,1.610,-384.00/)
      ABE(:,296) = (/2.870260e+01,.000,1808.00/)
      ABE(:,297) = (/2.870260e+01,.000,1808.00/)
      ABE(:,298) = (/3.103555e+01,.000,39150.00/)
      ABE(:,299) = (/2.144111e+01,1.160,2405.00/)
      ABE(:,300) = (/2.144111e+01,1.160,2405.00/)
      ABE(:,301) = (/2.387728e+01,0.730,-1113.00/)
      ABE(:,302) = (/2.873296e+01,.000,11923.00/)
      ABE(:,303) = (/1.481614e+01,1.770,5920.00/)
      ABE(:,304) = (/2.691050e+01,0.422,-1755.00/)
      ABE(:,305) = (/3.264166e+01,.000,.00/)
      ABE(:,306) = (/2.361918e+01,.000,.00/)
      ABE(:,307) = (/2.388027e+01,.000,.00/)
      ABE(:,308) = (/3.072206e+01,.000,.00/)
      ABE(:,309) = (/3.002892e+01,.000,.00/)
      ABE(:,310) = (/3.011593e+01,.000,.00/)
      ABE(:,311) = (/3.103555e+01,.000,.00/)
      ABE(:,312) = (/2.987492e+01,.000,.00/)
      ABE(:,313) = (/1.217045e+01,2.680,3716.00/)
      ABE(:,314) = (/1.409314e+01,2.540,6756.00/)
      ABE(:,315) = (/1.726867e+01,1.800,934.00/)
      ABE(:,316) = (/5.934894e+00,2.720,1500.00/)
      ABE(:,317) = (/-1.020327e-01,3.650,7154.00/)
      ABE(:,318) = (/1.475160e+01,1.600,5700.00/)
      ABE(:,319) = (/3.219953e+01,.000,.00/)
      ABE(:,320) = (/3.121814e+01,.000,.00/)
      ABE(:,321) = (/1.521669e+01,2.190,890.00/)
      ABE(:,322) = (/3.081323e+01,.000,.00/)
      ABE(:,323) = (/2.396194e+01,0.255,-943.00/)
      ABE(:,324) = (/3.081323e+01,.000,.00/)
      ABE(:,325) = (/3.058957e+01,-0.320,.00/)
      RETURN
      END

      SUBROUTINE GETNUNET(KK,II,NU)
      INTEGER KK,II
      INTEGER, dimension(KK,II) :: NU
      NU(1:10,1) = (/
     *0,0,-2,1,0,0,0,0,0,0/)
      NU(11:20,1) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,1) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,1) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,1) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,1) = (/
     *0,0,0/)
      NU(1:10,2) = (/
     *0,-1,-1,0,1,0,0,0,0,0/)
      NU(11:20,2) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,2) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,2) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,2) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,2) = (/
     *0,0,0/)
      NU(1:10,3) = (/
     *-1,1,-1,0,1,0,0,0,0,0/)
      NU(11:20,3) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,3) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,3) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,3) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,3) = (/
     *0,0,0/)
      NU(1:10,4) = (/
     *0,0,-1,1,1,0,-1,0,0,0/)
      NU(11:20,4) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,4) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,4) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,4) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,4) = (/
     *0,0,0/)
      NU(1:10,5) = (/
     *0,0,-1,0,1,0,1,-1,0,0/)
      NU(11:20,5) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,5) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,5) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,5) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,5) = (/
     *0,0,0/)
      NU(1:10,6) = (/
     *0,1,-1,0,0,0,0,0,0,-1/)
      NU(11:20,6) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,6) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,6) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,6) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,6) = (/
     *0,0,0/)
      NU(1:10,7) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,7) = (/
     *-1,0,0,0,0,0,1,0,0,0/)
      NU(21:30,7) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,7) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,7) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,7) = (/
     *0,0,0/)
      NU(1:10,8) = (/
     *1,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,8) = (/
     *0,-1,0,0,1,0,0,0,0,0/)
      NU(21:30,8) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,8) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,8) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,8) = (/
     *0,0,0/)
      NU(1:10,9) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,9) = (/
     *0,-1,0,0,0,0,1,0,0,0/)
      NU(21:30,9) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,9) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,9) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,9) = (/
     *0,0,0/)
      NU(1:10,10) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,10) = (/
     *0,0,-1,0,0,0,0,1,0,0/)
      NU(21:30,10) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,10) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,10) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,10) = (/
     *0,0,0/)
      NU(1:10,11) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,11) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(21:30,11) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,11) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,11) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,11) = (/
     *0,0,0/)
      NU(1:10,12) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,12) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(21:30,12) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,12) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,12) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,12) = (/
     *0,0,0/)
      NU(1:10,13) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,13) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(21:30,13) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,13) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,13) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,13) = (/
     *0,0,0/)
      NU(1:10,14) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,14) = (/
     *0,0,0,0,0,1,-1,0,0,0/)
      NU(21:30,14) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,14) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,14) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,14) = (/
     *0,0,0/)
      NU(1:10,15) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,15) = (/
     *0,0,0,0,0,0,1,-1,0,0/)
      NU(21:30,15) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,15) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,15) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,15) = (/
     *0,0,0/)
      NU(1:10,16) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,16) = (/
     *0,0,0,0,0,0,0,1,-1,0/)
      NU(21:30,16) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,16) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,16) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,16) = (/
     *0,0,0/)
      NU(1:10,17) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,17) = (/
     *0,0,0,0,0,0,0,1,0,-1/)
      NU(21:30,17) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,17) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,17) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,17) = (/
     *0,0,0/)
      NU(1:10,18) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,18) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(21:30,18) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,18) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,18) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,18) = (/
     *0,0,0/)
      NU(1:10,19) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,19) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(21:30,19) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,19) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,19) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,19) = (/
     *0,0,0/)
      NU(1:10,20) = (/
     *0,0,-1,0,0,0,0,0,0,1/)
      NU(11:20,20) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,20) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(31:40,20) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,20) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,20) = (/
     *0,0,0/)
      NU(1:10,21) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,21) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,21) = (/
     *0,0,-1,0,0,0,0,1,0,0/)
      NU(31:40,21) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,21) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,21) = (/
     *0,0,0/)
      NU(1:10,22) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,22) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,22) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(31:40,22) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,22) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,22) = (/
     *0,0,0/)
      NU(1:10,23) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,23) = (/
     *1,0,0,0,1,0,0,0,0,0/)
      NU(21:30,23) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(31:40,23) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,23) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,23) = (/
     *0,0,0/)
      NU(1:10,24) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,24) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,24) = (/
     *0,0,0,-1,0,0,0,0,1,0/)
      NU(31:40,24) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,24) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,24) = (/
     *0,0,0/)
      NU(1:10,25) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,25) = (/
     *0,0,1,0,0,0,1,0,0,0/)
      NU(21:30,25) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(31:40,25) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,25) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,25) = (/
     *0,0,0/)
      NU(1:10,26) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,26) = (/
     *0,0,1,0,0,0,0,1,0,0/)
      NU(21:30,26) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(31:40,26) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,26) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,26) = (/
     *0,0,0/)
      NU(1:10,27) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,27) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,27) = (/
     *0,0,0,0,0,1,-1,0,0,0/)
      NU(31:40,27) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,27) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,27) = (/
     *0,0,0/)
      NU(1:10,28) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,28) = (/
     *0,0,0,0,2,0,0,0,0,0/)
      NU(21:30,28) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(31:40,28) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,28) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,28) = (/
     *0,0,0/)
      NU(1:10,29) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,29) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,29) = (/
     *0,0,0,0,0,0,0,1,-1,0/)
      NU(31:40,29) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,29) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,29) = (/
     *0,0,0/)
      NU(1:10,30) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,30) = (/
     *1,0,0,0,0,1,0,0,0,0/)
      NU(21:30,30) = (/
     *0,0,0,0,0,0,0,0,-1,0/)
      NU(31:40,30) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,30) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,30) = (/
     *0,0,0/)
      NU(1:10,31) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(11:20,31) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(21:30,31) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,31) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,31) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,31) = (/
     *0,0,0/)
      NU(1:10,32) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,32) = (/
     *0,0,0,0,0,0,1,-1,0,0/)
      NU(21:30,32) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,32) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,32) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,32) = (/
     *0,0,0/)
      NU(1:10,33) = (/
     *0,-1,0,-1,0,0,1,0,0,0/)
      NU(11:20,33) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,33) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,33) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,33) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,33) = (/
     *0,0,0/)
      NU(1:10,34) = (/
     *0,-1,0,-1,0,0,1,0,0,0/)
      NU(11:20,34) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,34) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,34) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,34) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,34) = (/
     *0,0,0/)
      NU(1:10,35) = (/
     *0,-1,0,-1,0,0,1,0,0,0/)
      NU(11:20,35) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,35) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,35) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,35) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,35) = (/
     *0,0,0/)
      NU(1:10,36) = (/
     *0,-1,0,-1,0,0,1,0,0,0/)
      NU(11:20,36) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,36) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,36) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,36) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,36) = (/
     *0,0,0/)
      NU(1:10,37) = (/
     *0,-1,0,-1,0,0,1,0,0,0/)
      NU(11:20,37) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,37) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,37) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,37) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,37) = (/
     *0,0,0/)
      NU(1:10,38) = (/
     *0,-1,1,-1,1,0,0,0,0,0/)
      NU(11:20,38) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,38) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,38) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,38) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,38) = (/
     *0,0,0/)
      NU(1:10,39) = (/
     *1,-2,0,0,0,0,0,0,0,0/)
      NU(11:20,39) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,39) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,39) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,39) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,39) = (/
     *0,0,0/)
      NU(1:10,40) = (/
     *1,-2,0,0,0,0,0,0,0,0/)
      NU(11:20,40) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,40) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,40) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,40) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,40) = (/
     *0,0,0/)
      NU(1:10,41) = (/
     *1,-2,0,0,0,0,0,0,0,0/)
      NU(11:20,41) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,41) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,41) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,41) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,41) = (/
     *0,0,0/)
      NU(1:10,42) = (/
     *1,-2,0,0,0,0,0,0,0,0/)
      NU(11:20,42) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,42) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,42) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,42) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,42) = (/
     *0,0,0/)
      NU(1:10,43) = (/
     *0,-1,0,0,-1,1,0,0,0,0/)
      NU(11:20,43) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,43) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,43) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,43) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,43) = (/
     *0,0,0/)
      NU(1:10,44) = (/
     *0,-1,1,0,0,1,-1,0,0,0/)
      NU(11:20,44) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,44) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,44) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,44) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,44) = (/
     *0,0,0/)
      NU(1:10,45) = (/
     *1,-1,0,1,0,0,-1,0,0,0/)
      NU(11:20,45) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,45) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,45) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,45) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,45) = (/
     *0,0,0/)
      NU(1:10,46) = (/
     *0,-1,0,0,2,0,-1,0,0,0/)
      NU(11:20,46) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,46) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,46) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,46) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,46) = (/
     *0,0,0/)
      NU(1:10,47) = (/
     *1,-1,0,0,0,0,1,-1,0,0/)
      NU(11:20,47) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,47) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,47) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,47) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,47) = (/
     *0,0,0/)
      NU(1:10,48) = (/
     *0,-1,0,0,1,1,0,-1,0,0/)
      NU(11:20,48) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,48) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,48) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,48) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,48) = (/
     *0,0,0/)
      NU(1:10,49) = (/
     *1,-1,0,0,0,0,0,0,1,-1/)
      NU(11:20,49) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,49) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,49) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,49) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,49) = (/
     *0,0,0/)
      NU(1:10,50) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,50) = (/
     *-1,0,1,0,0,0,0,0,0,0/)
      NU(21:30,50) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,50) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,50) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,50) = (/
     *0,0,0/)
      NU(1:10,51) = (/
     *1,-1,0,0,0,0,0,0,0,1/)
      NU(11:20,51) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,51) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,51) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,51) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,51) = (/
     *0,0,0/)
      NU(1:10,52) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,52) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(21:30,52) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,52) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,52) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,52) = (/
     *0,0,0/)
      NU(1:10,53) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,53) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(21:30,53) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,53) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,53) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,53) = (/
     *0,0,0/)
      NU(1:10,54) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,54) = (/
     *0,0,0,0,0,0,-1,1,0,0/)
      NU(21:30,54) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,54) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,54) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,54) = (/
     *0,0,0/)
      NU(1:10,55) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,55) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(21:30,55) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,55) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,55) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,55) = (/
     *0,0,0/)
      NU(1:10,56) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,56) = (/
     *0,0,0,0,0,0,0,-1,1,0/)
      NU(21:30,56) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,56) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,56) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,56) = (/
     *0,0,0/)
      NU(1:10,57) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,57) = (/
     *0,0,0,0,0,0,0,-1,0,1/)
      NU(21:30,57) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,57) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,57) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,57) = (/
     *0,0,0/)
      NU(1:10,58) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,58) = (/
     *0,0,0,0,0,0,1,-1,0,0/)
      NU(21:30,58) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,58) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,58) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,58) = (/
     *0,0,0/)
      NU(1:10,59) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,59) = (/
     *0,0,0,0,0,0,0,0,-1,0/)
      NU(21:30,59) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,59) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,59) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,59) = (/
     *0,0,0/)
      NU(1:10,60) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,60) = (/
     *0,0,0,0,0,0,0,1,-1,0/)
      NU(21:30,60) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,60) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,60) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,60) = (/
     *0,0,0/)
      NU(1:10,61) = (/
     *0,-1,0,0,1,0,0,0,0,0/)
      NU(11:20,61) = (/
     *0,0,1,0,0,0,0,0,-1,0/)
      NU(21:30,61) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,61) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,61) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,61) = (/
     *0,0,0/)
      NU(1:10,62) = (/
     *0,-1,0,0,0,1,0,0,0,0/)
      NU(11:20,62) = (/
     *0,1,0,0,0,0,0,0,-1,0/)
      NU(21:30,62) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,62) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,62) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,62) = (/
     *0,0,0/)
      NU(1:10,63) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,63) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(21:30,63) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,63) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,63) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,63) = (/
     *0,0,0/)
      NU(1:10,64) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,64) = (/
     *0,0,0,0,0,0,0,0,1,-1/)
      NU(21:30,64) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,64) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,64) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,64) = (/
     *0,0,0/)
      NU(1:10,65) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,65) = (/
     *0,0,0,0,0,0,0,1,0,-1/)
      NU(21:30,65) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,65) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,65) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,65) = (/
     *0,0,0/)
      NU(1:10,66) = (/
     *0,-1,0,0,1,0,0,0,0,0/)
      NU(11:20,66) = (/
     *0,0,1,0,0,0,0,0,0,-1/)
      NU(21:30,66) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,66) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,66) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,66) = (/
     *0,0,0/)
      NU(1:10,67) = (/
     *0,-1,0,0,0,1,0,0,0,0/)
      NU(11:20,67) = (/
     *0,1,0,0,0,0,0,0,0,-1/)
      NU(21:30,67) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,67) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,67) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,67) = (/
     *0,0,0/)
      NU(1:10,68) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,68) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(21:30,68) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,68) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,68) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,68) = (/
     *0,0,0/)
      NU(1:10,69) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,69) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(21:30,69) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,69) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,69) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,69) = (/
     *0,0,0/)
      NU(1:10,70) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,70) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,70) = (/
     *0,-1,1,0,0,0,0,0,0,0/)
      NU(31:40,70) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,70) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,70) = (/
     *0,0,0/)
      NU(1:10,71) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,71) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,71) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(31:40,71) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,71) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,71) = (/
     *0,0,0/)
      NU(1:10,72) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,72) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,72) = (/
     *0,0,0,-1,1,0,0,0,0,0/)
      NU(31:40,72) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,72) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,72) = (/
     *0,0,0/)
      NU(1:10,73) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,73) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,73) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(31:40,73) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,73) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,73) = (/
     *0,0,0/)
      NU(1:10,74) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,74) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,74) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(31:40,74) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,74) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,74) = (/
     *0,0,0/)
      NU(1:10,75) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,75) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,75) = (/
     *0,0,0,1,-1,0,0,0,0,0/)
      NU(31:40,75) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,75) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,75) = (/
     *0,0,0/)
      NU(1:10,76) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,76) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,76) = (/
     *0,0,0,0,0,-1,1,0,0,0/)
      NU(31:40,76) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,76) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,76) = (/
     *0,0,0/)
      NU(1:10,77) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,77) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,77) = (/
     *0,0,0,0,1,-1,0,0,0,0/)
      NU(31:40,77) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,77) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,77) = (/
     *0,0,0/)
      NU(1:10,78) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,78) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,78) = (/
     *0,0,0,0,0,1,-1,0,0,0/)
      NU(31:40,78) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,78) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,78) = (/
     *0,0,0/)
      NU(1:10,79) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,79) = (/
     *0,1,0,0,1,0,0,0,0,0/)
      NU(21:30,79) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(31:40,79) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,79) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,79) = (/
     *0,0,0/)
      NU(1:10,80) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,80) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,80) = (/
     *0,0,0,0,0,0,0,1,-1,0/)
      NU(31:40,80) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,80) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,80) = (/
     *0,0,0/)
      NU(1:10,81) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,81) = (/
     *0,0,1,0,1,0,0,0,0,0/)
      NU(21:30,81) = (/
     *0,0,0,0,0,0,0,0,-1,0/)
      NU(31:40,81) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,81) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,81) = (/
     *0,0,0/)
      NU(1:10,82) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,82) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,82) = (/
     *0,0,0,0,0,0,0,0,1,-1/)
      NU(31:40,82) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,82) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,82) = (/
     *0,0,0/)
      NU(1:10,83) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(11:20,83) = (/
     *0,0,0,0,-1,0,0,1,0,0/)
      NU(21:30,83) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,83) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,83) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,83) = (/
     *0,0,0/)
      NU(1:10,84) = (/
     *-1,1,0,0,-1,1,0,0,0,0/)
      NU(11:20,84) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,84) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,84) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,84) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,84) = (/
     *0,0,0/)
      NU(1:10,85) = (/
     *0,0,0,0,-2,0,0,1,0,0/)
      NU(11:20,85) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,85) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,85) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,85) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,85) = (/
     *0,0,0/)
      NU(1:10,86) = (/
     *0,0,1,0,-2,1,0,0,0,0/)
      NU(11:20,86) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,86) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,86) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,86) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,86) = (/
     *0,0,0/)
      NU(1:10,87) = (/
     *0,0,0,1,-1,1,-1,0,0,0/)
      NU(11:20,87) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,87) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,87) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,87) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,87) = (/
     *0,0,0/)
      NU(1:10,88) = (/
     *0,0,0,0,-1,1,1,-1,0,0/)
      NU(11:20,88) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,88) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,88) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,88) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,88) = (/
     *0,0,0/)
      NU(1:10,89) = (/
     *0,0,0,0,-1,1,1,-1,0,0/)
      NU(11:20,89) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,89) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,89) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,89) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,89) = (/
     *0,0,0/)
      NU(1:10,90) = (/
     *0,1,0,0,-1,0,0,0,-1,0/)
      NU(11:20,90) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,90) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,90) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,90) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,90) = (/
     *0,0,0/)
      NU(1:10,91) = (/
     *0,1,0,0,-1,0,0,0,0,-1/)
      NU(11:20,91) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(21:30,91) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,91) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,91) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,91) = (/
     *0,0,0/)
      NU(1:10,92) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,92) = (/
     *-1,0,0,0,0,0,0,1,0,0/)
      NU(21:30,92) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,92) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,92) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,92) = (/
     *0,0,0/)
      NU(1:10,93) = (/
     *0,0,0,0,-1,1,0,0,0,1/)
      NU(11:20,93) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,93) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,93) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,93) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,93) = (/
     *0,0,0/)
      NU(1:10,94) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,94) = (/
     *0,-1,0,0,0,0,0,1,0,0/)
      NU(21:30,94) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,94) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,94) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,94) = (/
     *0,0,0/)
      NU(1:10,95) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(11:20,95) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,95) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,95) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,95) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,95) = (/
     *0,0,0/)
      NU(1:10,96) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,96) = (/
     *1,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,96) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,96) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,96) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,96) = (/
     *0,0,0/)
      NU(1:10,97) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,97) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(21:30,97) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,97) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,97) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,97) = (/
     *0,0,0/)
      NU(1:10,98) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,98) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(21:30,98) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,98) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,98) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,98) = (/
     *0,0,0/)
      NU(1:10,99) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,99) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(21:30,99) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,99) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,99) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,99) = (/
     *0,0,0/)
      NU(1:10,100) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,100) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(21:30,100) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,100) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,100) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,100) = (/
     *0,0,0/)
      NU(1:10,101) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,101) = (/
     *0,0,0,0,0,0,1,-1,0,0/)
      NU(21:30,101) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,101) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,101) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,101) = (/
     *0,0,0/)
      NU(1:10,102) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,102) = (/
     *0,0,0,0,0,0,0,1,-1,0/)
      NU(21:30,102) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,102) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,102) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,102) = (/
     *0,0,0/)
      NU(1:10,103) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,103) = (/
     *0,0,0,0,0,0,0,1,0,-1/)
      NU(21:30,103) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,103) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,103) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,103) = (/
     *0,0,0/)
      NU(1:10,104) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,104) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(21:30,104) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,104) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,104) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,104) = (/
     *0,0,0/)
      NU(1:10,105) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,105) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(21:30,105) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,105) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,105) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,105) = (/
     *0,0,0/)
      NU(1:10,106) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,106) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,106) = (/
     *0,-1,0,0,0,0,0,1,0,0/)
      NU(31:40,106) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,106) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,106) = (/
     *0,0,0/)
      NU(1:10,107) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,107) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,107) = (/
     *0,0,-1,0,0,0,0,0,1,0/)
      NU(31:40,107) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,107) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,107) = (/
     *0,0,0/)
      NU(1:10,108) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,108) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,108) = (/
     *0,0,-1,0,0,0,0,0,0,1/)
      NU(31:40,108) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,108) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,108) = (/
     *0,0,0/)
      NU(1:10,109) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,109) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,109) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(31:40,109) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,109) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,109) = (/
     *0,0,0/)
      NU(1:10,110) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(11:20,110) = (/
     *0,0,1,0,1,0,0,0,0,0/)
      NU(21:30,110) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(31:40,110) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,110) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,110) = (/
     *0,0,0/)
      NU(1:10,111) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,111) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,111) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(31:40,111) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,111) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,111) = (/
     *0,0,0/)
      NU(1:10,112) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,112) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,112) = (/
     *0,0,0,1,-1,0,0,0,0,0/)
      NU(31:40,112) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,112) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,112) = (/
     *0,0,0/)
      NU(1:10,113) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,113) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,113) = (/
     *0,0,0,0,0,1,-1,0,0,0/)
      NU(31:40,113) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,113) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,113) = (/
     *0,0,0/)
      NU(1:10,114) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,114) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,114) = (/
     *0,0,0,0,0,0,0,1,-1,0/)
      NU(31:40,114) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,114) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,114) = (/
     *0,0,0/)
      NU(1:10,115) = (/
     *0,0,0,1,0,0,-2,1,0,0/)
      NU(11:20,115) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,115) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,115) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,115) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,115) = (/
     *0,0,0/)
      NU(1:10,116) = (/
     *0,0,0,1,0,0,-2,1,0,0/)
      NU(11:20,116) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,116) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,116) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,116) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,116) = (/
     *0,0,0/)
      NU(1:10,117) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(11:20,117) = (/
     *-1,0,0,0,0,0,0,1,0,0/)
      NU(21:30,117) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,117) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,117) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,117) = (/
     *0,0,0/)
      NU(1:10,118) = (/
     *0,0,0,1,0,0,-1,0,0,0/)
      NU(11:20,118) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(21:30,118) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,118) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,118) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,118) = (/
     *0,0,0/)
      NU(1:10,119) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(11:20,119) = (/
     *0,0,-1,0,0,0,0,0,0,1/)
      NU(21:30,119) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,119) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,119) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,119) = (/
     *0,0,0/)
      NU(1:10,120) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(11:20,120) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(21:30,120) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,120) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,120) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,120) = (/
     *0,0,0/)
      NU(1:10,121) = (/
     *0,0,0,0,0,0,-1,1,0,0/)
      NU(11:20,121) = (/
     *0,0,0,0,0,0,1,-1,0,0/)
      NU(21:30,121) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,121) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,121) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,121) = (/
     *0,0,0/)
      NU(1:10,122) = (/
     *0,0,1,-1,0,0,0,0,-1,0/)
      NU(11:20,122) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,122) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,122) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,122) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,122) = (/
     *0,0,0/)
      NU(1:10,123) = (/
     *0,1,0,0,0,0,0,0,-1,0/)
      NU(11:20,123) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,123) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(31:40,123) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,123) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,123) = (/
     *0,0,0/)
      NU(1:10,124) = (/
     *0,1,0,0,0,0,0,0,-1,0/)
      NU(11:20,124) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,124) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(31:40,124) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,124) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,124) = (/
     *0,0,0/)
      NU(1:10,125) = (/
     *0,0,1,-1,0,0,0,0,0,-1/)
      NU(11:20,125) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(21:30,125) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,125) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,125) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,125) = (/
     *0,0,0/)
      NU(1:10,126) = (/
     *-1,1,0,0,0,0,0,0,0,-1/)
      NU(11:20,126) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,126) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,126) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,126) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,126) = (/
     *0,0,0/)
      NU(1:10,127) = (/
     *0,1,0,0,0,-1,0,0,0,-1/)
      NU(11:20,127) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(21:30,127) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,127) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,127) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,127) = (/
     *0,0,0/)
      NU(1:10,128) = (/
     *0,1,0,0,0,0,0,0,0,-1/)
      NU(11:20,128) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,128) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(31:40,128) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,128) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,128) = (/
     *0,0,0/)
      NU(1:10,129) = (/
     *0,1,0,0,0,0,0,0,0,-1/)
      NU(11:20,129) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,129) = (/
     *0,0,0,1,0,0,0,0,0,0/)
      NU(31:40,129) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,129) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,129) = (/
     *0,0,0/)
      NU(1:10,130) = (/
     *0,1,0,0,0,0,0,0,0,-1/)
      NU(11:20,130) = (/
     *0,0,0,-1,0,0,0,0,0,0/)
      NU(21:30,130) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(31:40,130) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,130) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,130) = (/
     *0,0,0/)
      NU(1:10,131) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(11:20,131) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(21:30,131) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(31:40,131) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,131) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,131) = (/
     *0,0,0/)
      NU(1:10,132) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(11:20,132) = (/
     *0,0,0,0,1,-1,1,0,0,0/)
      NU(21:30,132) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,132) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,132) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,132) = (/
     *0,0,0/)
      NU(1:10,133) = (/
     *0,1,0,0,0,0,0,0,0,-1/)
      NU(11:20,133) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(21:30,133) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(31:40,133) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,133) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,133) = (/
     *0,0,0/)
      NU(1:10,134) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(11:20,134) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,134) = (/
     *0,0,1,0,0,0,0,-1,0,0/)
      NU(31:40,134) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,134) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,134) = (/
     *0,0,0/)
      NU(1:10,135) = (/
     *0,1,0,-1,1,0,0,0,0,0/)
      NU(11:20,135) = (/
     *-1,0,0,0,1,0,0,0,0,0/)
      NU(21:30,135) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,135) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,135) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,135) = (/
     *0,0,0/)
      NU(1:10,136) = (/
     *-1,1,0,0,0,0,0,0,0,0/)
      NU(11:20,136) = (/
     *-1,0,1,0,0,0,0,0,0,0/)
      NU(21:30,136) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,136) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,136) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,136) = (/
     *0,0,0/)
      NU(1:10,137) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(11:20,137) = (/
     *-2,0,0,0,0,0,0,0,0,0/)
      NU(21:30,137) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(31:40,137) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,137) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,137) = (/
     *0,0,0/)
      NU(1:10,138) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,138) = (/
     *-1,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,138) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(31:40,138) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,138) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,138) = (/
     *0,0,0/)
      NU(1:10,139) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,139) = (/
     *-1,0,2,-1,0,0,0,0,0,0/)
      NU(21:30,139) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,139) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,139) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,139) = (/
     *0,0,0/)
      NU(1:10,140) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,140) = (/
     *-1,0,0,0,-1,0,0,0,0,0/)
      NU(21:30,140) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(31:40,140) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,140) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,140) = (/
     *0,0,0/)
      NU(1:10,141) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,141) = (/
     *-1,0,0,0,1,0,0,0,0,0/)
      NU(21:30,141) = (/
     *0,0,0,1,0,0,0,-1,0,0/)
      NU(31:40,141) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,141) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,141) = (/
     *0,0,0/)
      NU(1:10,142) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,142) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,142) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,142) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,142) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,142) = (/
     *0,0,0/)
      NU(1:10,143) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,143) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,143) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,143) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,143) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,143) = (/
     *0,0,0/)
      NU(1:10,144) = (/
     *0,1,0,-1,1,0,0,0,0,0/)
      NU(11:20,144) = (/
     *0,-1,0,0,1,0,0,0,0,0/)
      NU(21:30,144) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,144) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,144) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,144) = (/
     *0,0,0/)
      NU(1:10,145) = (/
     *0,0,0,-1,0,1,0,0,0,0/)
      NU(11:20,145) = (/
     *0,-1,0,0,1,0,0,0,0,0/)
      NU(21:30,145) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,145) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,145) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,145) = (/
     *0,0,0/)
      NU(1:10,146) = (/
     *-1,1,0,0,0,0,0,0,0,0/)
      NU(11:20,146) = (/
     *0,-1,1,0,0,0,0,0,0,0/)
      NU(21:30,146) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,146) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,146) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,146) = (/
     *0,0,0/)
      NU(1:10,147) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(11:20,147) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,147) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,147) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,147) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,147) = (/
     *0,0,0/)
      NU(1:10,148) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,148) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,148) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,148) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,148) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,148) = (/
     *0,0,0/)
      NU(1:10,149) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,149) = (/
     *0,-1,-1,0,0,0,0,0,0,0/)
      NU(21:30,149) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(31:40,149) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,149) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,149) = (/
     *0,0,0/)
      NU(1:10,150) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,150) = (/
     *0,-1,2,-1,0,0,0,0,0,0/)
      NU(21:30,150) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,150) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,150) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,150) = (/
     *0,0,0/)
      NU(1:10,151) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,151) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,151) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,151) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,151) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,151) = (/
     *0,0,0/)
      NU(1:10,152) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,152) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,152) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,152) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,152) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,152) = (/
     *0,0,0/)
      NU(1:10,153) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,153) = (/
     *0,-1,0,0,1,-1,0,1,0,0/)
      NU(21:30,153) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,153) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,153) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,153) = (/
     *0,0,0/)
      NU(1:10,154) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,154) = (/
     *0,-1,1,0,0,0,0,0,0,0/)
      NU(21:30,154) = (/
     *0,0,0,0,0,1,-1,0,0,0/)
      NU(31:40,154) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,154) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,154) = (/
     *0,0,0/)
      NU(1:10,155) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(11:20,155) = (/
     *0,0,-1,0,0,0,0,0,0,1/)
      NU(21:30,155) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,155) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,155) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,155) = (/
     *0,0,0/)
      NU(1:10,156) = (/
     *0,0,0,-1,1,0,0,0,0,0/)
      NU(11:20,156) = (/
     *0,0,-1,0,0,0,0,1,0,0/)
      NU(21:30,156) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,156) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,156) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,156) = (/
     *0,0,0/)
      NU(1:10,157) = (/
     *0,0,0,0,0,0,1,-1,0,0/)
      NU(11:20,157) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(21:30,157) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,157) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,157) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,157) = (/
     *0,0,0/)
      NU(1:10,158) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,158) = (/
     *0,0,-2,0,0,0,0,0,0,0/)
      NU(21:30,158) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(31:40,158) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,158) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,158) = (/
     *0,0,0/)
      NU(1:10,159) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,159) = (/
     *0,0,-2,0,0,0,0,0,0,0/)
      NU(21:30,159) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(31:40,159) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,159) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,159) = (/
     *0,0,0/)
      NU(1:10,160) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,160) = (/
     *0,0,-1,1,1,0,-1,0,0,0/)
      NU(21:30,160) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,160) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,160) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,160) = (/
     *0,0,0/)
      NU(1:10,161) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,161) = (/
     *0,0,-1,1,0,0,1,-1,0,0/)
      NU(21:30,161) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,161) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,161) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,161) = (/
     *0,0,0/)
      NU(1:10,162) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,162) = (/
     *0,0,-1,1,0,0,0,0,1,0/)
      NU(21:30,162) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,162) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,162) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,162) = (/
     *0,0,0/)
      NU(1:10,163) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,163) = (/
     *0,0,-1,1,0,0,0,0,0,1/)
      NU(21:30,163) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(31:40,163) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,163) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,163) = (/
     *0,0,0/)
      NU(1:10,164) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,164) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(21:30,164) = (/
     *0,0,0,1,-1,0,0,0,0,0/)
      NU(31:40,164) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,164) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,164) = (/
     *0,0,0/)
      NU(1:10,165) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,165) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(21:30,165) = (/
     *0,0,0,0,0,1,-1,0,0,0/)
      NU(31:40,165) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,165) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,165) = (/
     *0,0,0/)
      NU(1:10,166) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,166) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(21:30,166) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,166) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,166) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,166) = (/
     *0,0,0/)
      NU(1:10,167) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,167) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(21:30,167) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,167) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,167) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,167) = (/
     *0,0,0/)
      NU(1:10,168) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,168) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(21:30,168) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,168) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,168) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,168) = (/
     *0,0,0/)
      NU(1:10,169) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,169) = (/
     *0,0,0,0,0,0,0,1,-1,0/)
      NU(21:30,169) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,169) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,169) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,169) = (/
     *0,0,0/)
      NU(1:10,170) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,170) = (/
     *0,0,0,0,0,0,0,1,0,-1/)
      NU(21:30,170) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,170) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,170) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,170) = (/
     *0,0,0/)
      NU(1:10,171) = (/
     *0,0,0,-1,0,0,0,0,0,0/)
      NU(11:20,171) = (/
     *0,0,0,0,1,0,1,0,0,0/)
      NU(21:30,171) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(31:40,171) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,171) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,171) = (/
     *0,0,0/)
      NU(1:10,172) = (/
     *-1,1,0,0,0,0,0,0,0,0/)
      NU(11:20,172) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,172) = (/
     *0,-1,1,0,0,0,0,0,0,0/)
      NU(31:40,172) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,172) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,172) = (/
     *0,0,0/)
      NU(1:10,173) = (/
     *0,0,0,-1,0,0,0,0,0,0/)
      NU(11:20,173) = (/
     *0,0,0,0,0,0,1,1,0,0/)
      NU(21:30,173) = (/
     *0,0,0,-1,0,0,0,0,0,0/)
      NU(31:40,173) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,173) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,173) = (/
     *0,0,0/)
      NU(1:10,174) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(11:20,174) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,174) = (/
     *0,0,1,0,-1,0,0,0,0,0/)
      NU(31:40,174) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,174) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,174) = (/
     *0,0,0/)
      NU(1:10,175) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,175) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,175) = (/
     *0,0,0,0,1,-1,0,0,0,0/)
      NU(31:40,175) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,175) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,175) = (/
     *0,0,0/)
      NU(1:10,176) = (/
     *0,0,0,-1,1,0,0,0,0,0/)
      NU(11:20,176) = (/
     *0,0,0,0,2,0,0,0,0,0/)
      NU(21:30,176) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(31:40,176) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,176) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,176) = (/
     *0,0,0/)
      NU(1:10,177) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,177) = (/
     *0,0,0,0,2,0,0,0,0,0/)
      NU(21:30,177) = (/
     *0,0,1,0,0,0,0,-2,0,0/)
      NU(31:40,177) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,177) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,177) = (/
     *0,0,0/)
      NU(1:10,178) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(11:20,178) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,178) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,178) = (/
     *-1,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,178) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,178) = (/
     *0,0,0/)
      NU(1:10,179) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(11:20,179) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,179) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,179) = (/
     *-1,0,0,0,0,1,0,0,0,0/)
      NU(41:50,179) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,179) = (/
     *0,0,0/)
      NU(1:10,180) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,180) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,180) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,180) = (/
     *-1,0,0,0,0,1,0,0,0,0/)
      NU(41:50,180) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,180) = (/
     *0,0,0/)
      NU(1:10,181) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(11:20,181) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,181) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,181) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(41:50,181) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,181) = (/
     *0,0,0/)
      NU(1:10,182) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,182) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,182) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,182) = (/
     *0,0,0,0,0,2,0,-1,0,0/)
      NU(41:50,182) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,182) = (/
     *0,0,0/)
      NU(1:10,183) = (/
     *0,-1,0,0,1,0,0,0,0,0/)
      NU(11:20,183) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,183) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,183) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(41:50,183) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,183) = (/
     *0,0,0/)
      NU(1:10,184) = (/
     *0,0,0,0,-1,0,1,0,0,0/)
      NU(11:20,184) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,184) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,184) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(41:50,184) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,184) = (/
     *0,0,0/)
      NU(1:10,185) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(11:20,185) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,185) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,185) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(41:50,185) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,185) = (/
     *0,0,0/)
      NU(1:10,186) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(11:20,186) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,186) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,186) = (/
     *0,0,0,0,0,-1,1,0,0,0/)
      NU(41:50,186) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,186) = (/
     *0,0,0/)
      NU(1:10,187) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,187) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,187) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,187) = (/
     *0,0,0,0,0,-1,1,0,0,0/)
      NU(41:50,187) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,187) = (/
     *0,0,0/)
      NU(1:10,188) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(11:20,188) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,188) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,188) = (/
     *0,0,0,0,0,1,-1,0,0,0/)
      NU(41:50,188) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,188) = (/
     *0,0,0/)
      NU(1:10,189) = (/
     *0,-1,0,0,1,0,0,0,0,0/)
      NU(11:20,189) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,189) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,189) = (/
     *0,0,0,0,0,1,-1,0,0,0/)
      NU(41:50,189) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,189) = (/
     *0,0,0/)
      NU(1:10,190) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,190) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,190) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,190) = (/
     *0,-1,0,0,0,1,0,0,0,0/)
      NU(41:50,190) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,190) = (/
     *0,0,0/)
      NU(1:10,191) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,191) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,191) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,191) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(41:50,191) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,191) = (/
     *0,0,0/)
      NU(1:10,192) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,192) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,192) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,192) = (/
     *0,-1,0,0,0,0,0,0,1,0/)
      NU(41:50,192) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,192) = (/
     *0,0,0/)
      NU(1:10,193) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,193) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,193) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,193) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(41:50,193) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,193) = (/
     *0,0,0/)
      NU(1:10,194) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(11:20,194) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,194) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,194) = (/
     *0,-1,0,0,0,0,0,0,1,0/)
      NU(41:50,194) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,194) = (/
     *0,0,0/)
      NU(1:10,195) = (/
     *0,0,0,-1,1,0,0,0,0,0/)
      NU(11:20,195) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,195) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,195) = (/
     *0,-1,0,0,0,1,0,0,0,0/)
      NU(41:50,195) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,195) = (/
     *0,0,0/)
      NU(1:10,196) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,196) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,196) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,196) = (/
     *-1,-1,0,0,0,0,0,0,0,0/)
      NU(41:50,196) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,196) = (/
     *0,0,0/)
      NU(1:10,197) = (/
     *1,0,0,0,0,-1,0,0,0,0/)
      NU(11:20,197) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,197) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,197) = (/
     *0,-1,0,0,0,0,0,0,1,0/)
      NU(41:50,197) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,197) = (/
     *0,0,0/)
      NU(1:10,198) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(11:20,198) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,198) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,198) = (/
     *0,-1,0,0,0,-1,0,0,0,0/)
      NU(41:50,198) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,198) = (/
     *0,0,0/)
      NU(1:10,199) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,199) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,199) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,199) = (/
     *0,-1,0,0,0,-1,0,1,0,0/)
      NU(41:50,199) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,199) = (/
     *0,0,0/)
      NU(1:10,200) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,200) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,200) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,200) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(41:50,200) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,200) = (/
     *0,0,0/)
      NU(1:10,201) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,201) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,201) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,201) = (/
     *0,0,-1,0,0,0,0,0,1,0/)
      NU(41:50,201) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,201) = (/
     *0,0,0/)
      NU(1:10,202) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,202) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,202) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,202) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(41:50,202) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,202) = (/
     *0,0,0/)
      NU(1:10,203) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,203) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,203) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,203) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(41:50,203) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,203) = (/
     *0,0,0/)
      NU(1:10,204) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,204) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,204) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,204) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(41:50,204) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,204) = (/
     *0,0,0/)
      NU(1:10,205) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,205) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,205) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,205) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(41:50,205) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,205) = (/
     *0,0,0/)
      NU(1:10,206) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,206) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,206) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,206) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(41:50,206) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,206) = (/
     *0,0,0/)
      NU(1:10,207) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,207) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,207) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,207) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(41:50,207) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,207) = (/
     *0,0,0/)
      NU(1:10,208) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,208) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,208) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,208) = (/
     *0,1,0,0,-1,1,0,0,0,0/)
      NU(41:50,208) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,208) = (/
     *0,0,0/)
      NU(1:10,209) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,209) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,209) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,209) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(41:50,209) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,209) = (/
     *0,0,0/)
      NU(1:10,210) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,210) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,210) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,210) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(41:50,210) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,210) = (/
     *0,0,0/)
      NU(1:10,211) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,211) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(21:30,211) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,211) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(41:50,211) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(51:53,211) = (/
     *0,0,0/)
      NU(1:10,212) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,212) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,212) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,212) = (/
     *0,0,0,0,0,-1,0,0,1,0/)
      NU(41:50,212) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,212) = (/
     *0,0,0/)
      NU(1:10,213) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,213) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,213) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,213) = (/
     *0,0,0,0,0,1,0,0,-1,0/)
      NU(41:50,213) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,213) = (/
     *0,0,0/)
      NU(1:10,214) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,214) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,214) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,214) = (/
     *0,0,0,0,0,1,0,0,-1,0/)
      NU(41:50,214) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,214) = (/
     *0,0,0/)
      NU(1:10,215) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,215) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,215) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,215) = (/
     *0,0,0,0,0,1,0,0,-1,0/)
      NU(41:50,215) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,215) = (/
     *0,0,0/)
      NU(1:10,216) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,216) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,216) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,216) = (/
     *0,0,0,0,0,1,0,0,-1,0/)
      NU(41:50,216) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,216) = (/
     *0,0,0/)
      NU(1:10,217) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,217) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,217) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,217) = (/
     *1,0,0,0,0,0,0,0,0,-1/)
      NU(41:50,217) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,217) = (/
     *0,0,0/)
      NU(1:10,218) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,218) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,218) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,218) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(41:50,218) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(51:53,218) = (/
     *0,0,0/)
      NU(1:10,219) = (/
     *0,0,0,0,1,-1,0,0,0,0/)
      NU(11:20,219) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,219) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,219) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(41:50,219) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,219) = (/
     *0,0,0/)
      NU(1:10,220) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(11:20,220) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,220) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,220) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(41:50,220) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(51:53,220) = (/
     *0,0,0/)
      NU(1:10,221) = (/
     *-1,1,0,0,0,0,0,0,0,0/)
      NU(11:20,221) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,221) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,221) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(41:50,221) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,221) = (/
     *0,0,0/)
      NU(1:10,222) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,222) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,222) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,222) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(41:50,222) = (/
     *0,0,0,0,0,0,-1,0,0,0/)
      NU(51:53,222) = (/
     *0,0,0/)
      NU(1:10,223) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,223) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,223) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,223) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(41:50,223) = (/
     *0,0,0,0,0,0,-1,0,0,0/)
      NU(51:53,223) = (/
     *0,0,0/)
      NU(1:10,224) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,224) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,224) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,224) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(41:50,224) = (/
     *0,0,0,0,0,0,-1,0,0,0/)
      NU(51:53,224) = (/
     *0,0,0/)
      NU(1:10,225) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,225) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,225) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,225) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(41:50,225) = (/
     *0,0,0,0,0,0,-1,1,0,0/)
      NU(51:53,225) = (/
     *0,0,0/)
      NU(1:10,226) = (/
     *0,0,0,-1,0,0,0,0,0,0/)
      NU(11:20,226) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(21:30,226) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,226) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(41:50,226) = (/
     *0,0,0,0,0,0,-1,0,0,0/)
      NU(51:53,226) = (/
     *0,0,0/)
      NU(1:10,227) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,227) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,227) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,227) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(41:50,227) = (/
     *0,0,0,0,0,0,-1,0,0,0/)
      NU(51:53,227) = (/
     *0,0,0/)
      NU(1:10,228) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,228) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,228) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,228) = (/
     *0,0,0,0,0,-1,0,1,0,0/)
      NU(41:50,228) = (/
     *0,0,0,0,0,0,-1,0,0,0/)
      NU(51:53,228) = (/
     *0,0,0/)
      NU(1:10,229) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,229) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(21:30,229) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,229) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,229) = (/
     *0,0,0,0,0,0,-1,1,0,0/)
      NU(51:53,229) = (/
     *0,0,0/)
      NU(1:10,230) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,230) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,230) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,230) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(41:50,230) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,230) = (/
     *0,0,0/)
      NU(1:10,231) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,231) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,231) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,231) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,231) = (/
     *-1,0,0,0,0,0,1,0,0,0/)
      NU(51:53,231) = (/
     *0,0,0/)
      NU(1:10,232) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,232) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,232) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,232) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(41:50,232) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,232) = (/
     *0,0,0/)
      NU(1:10,233) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,233) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,233) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,233) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(41:50,233) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,233) = (/
     *0,0,0/)
      NU(1:10,234) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,234) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,234) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,234) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,234) = (/
     *-1,0,0,0,1,0,0,0,0,0/)
      NU(51:53,234) = (/
     *0,0,0/)
      NU(1:10,235) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,235) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,235) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,235) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,235) = (/
     *-1,0,0,0,0,1,0,0,0,0/)
      NU(51:53,235) = (/
     *0,0,0/)
      NU(1:10,236) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(11:20,236) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,236) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,236) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(41:50,236) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,236) = (/
     *0,0,0/)
      NU(1:10,237) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,237) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,237) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,237) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,237) = (/
     *-1,1,0,0,0,0,0,0,0,0/)
      NU(51:53,237) = (/
     *0,0,0/)
      NU(1:10,238) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,238) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,238) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,238) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(41:50,238) = (/
     *0,-1,0,0,0,0,0,1,0,0/)
      NU(51:53,238) = (/
     *0,0,0/)
      NU(1:10,239) = (/
     *0,0,0,0,0,0,0,0,-1,0/)
      NU(11:20,239) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,239) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,239) = (/
     *1,0,0,0,0,0,0,0,0,1/)
      NU(41:50,239) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(51:53,239) = (/
     *0,0,0/)
      NU(1:10,240) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(11:20,240) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,240) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,240) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(41:50,240) = (/
     *1,0,0,0,0,0,0,-1,0,0/)
      NU(51:53,240) = (/
     *0,0,0/)
      NU(1:10,241) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(11:20,241) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,241) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,241) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,241) = (/
     *0,0,1,0,0,0,0,-1,0,0/)
      NU(51:53,241) = (/
     *0,0,0/)
      NU(1:10,242) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,242) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,242) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,242) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(41:50,242) = (/
     *1,0,0,0,0,0,0,-1,0,0/)
      NU(51:53,242) = (/
     *0,0,0/)
      NU(1:10,243) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,243) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,243) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,243) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(41:50,243) = (/
     *1,0,0,0,0,0,0,-1,0,0/)
      NU(51:53,243) = (/
     *0,0,0/)
      NU(1:10,244) = (/
     *0,0,1,0,0,0,0,0,-1,0/)
      NU(11:20,244) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,244) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,244) = (/
     *0,0,0,0,0,-1,0,0,0,1/)
      NU(41:50,244) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,244) = (/
     *0,0,0/)
      NU(1:10,245) = (/
     *0,0,0,0,0,0,0,0,-1,0/)
      NU(11:20,245) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,245) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,245) = (/
     *1,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,245) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,245) = (/
     *0,0,0/)
      NU(1:10,246) = (/
     *0,0,1,0,0,0,0,0,0,-1/)
      NU(11:20,246) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,246) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,246) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,246) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,246) = (/
     *0,0,0/)
      NU(1:10,247) = (/
     *0,1,0,0,0,0,0,0,0,-1/)
      NU(11:20,247) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,247) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,247) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,247) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(51:53,247) = (/
     *0,0,0/)
      NU(1:10,248) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(11:20,248) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(21:30,248) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,248) = (/
     *1,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,248) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,248) = (/
     *0,0,0/)
      NU(1:10,249) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,249) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,249) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,249) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,249) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(51:53,249) = (/
     *0,0,0/)
      NU(1:10,250) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(11:20,250) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,250) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,250) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,250) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,250) = (/
     *0,0,0/)
      NU(1:10,251) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,251) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,251) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,251) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,251) = (/
     *0,0,0,1,0,0,0,0,0,0/)
      NU(51:53,251) = (/
     *0,0,0/)
      NU(1:10,252) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,252) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,252) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,252) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,252) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(51:53,252) = (/
     *0,0,0/)
      NU(1:10,253) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(11:20,253) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,253) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,253) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,253) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,253) = (/
     *0,0,0/)
      NU(1:10,254) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,254) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(21:30,254) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,254) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,254) = (/
     *0,0,0,1,0,0,0,0,0,0/)
      NU(51:53,254) = (/
     *0,0,0/)
      NU(1:10,255) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(11:20,255) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,255) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,255) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,255) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,255) = (/
     *0,0,0/)
      NU(1:10,256) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(11:20,256) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,256) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,256) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,256) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(51:53,256) = (/
     *0,0,0/)
      NU(1:10,257) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,257) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,257) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,257) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,257) = (/
     *0,0,-1,0,0,0,0,1,0,0/)
      NU(51:53,257) = (/
     *0,0,0/)
      NU(1:10,258) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,258) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,258) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,258) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(41:50,258) = (/
     *1,0,-1,0,0,0,0,0,0,0/)
      NU(51:53,258) = (/
     *0,0,0/)
      NU(1:10,259) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(11:20,259) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(21:30,259) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,259) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,259) = (/
     *0,0,-1,0,0,0,0,1,0,0/)
      NU(51:53,259) = (/
     *0,0,0/)
      NU(1:10,260) = (/
     *0,1,0,0,-1,0,0,0,0,0/)
      NU(11:20,260) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(21:30,260) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,260) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,260) = (/
     *0,0,-1,0,0,0,0,1,0,0/)
      NU(51:53,260) = (/
     *0,0,0/)
      NU(1:10,261) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,261) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(21:30,261) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,261) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,261) = (/
     *0,0,-1,0,0,0,0,1,0,0/)
      NU(51:53,261) = (/
     *0,0,0/)
      NU(1:10,262) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,262) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(21:30,262) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,262) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(41:50,262) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(51:53,262) = (/
     *0,0,0/)
      NU(1:10,263) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,263) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,263) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,263) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(41:50,263) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(51:53,263) = (/
     *0,0,0/)
      NU(1:10,264) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,264) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,264) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,264) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,264) = (/
     *0,0,0,0,0,-1,1,0,0,0/)
      NU(51:53,264) = (/
     *0,0,0/)
      NU(1:10,265) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,265) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,265) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,265) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(41:50,265) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(51:53,265) = (/
     *0,0,0/)
      NU(1:10,266) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,266) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,266) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,266) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,266) = (/
     *0,0,0,0,0,-1,1,0,0,0/)
      NU(51:53,266) = (/
     *0,0,0/)
      NU(1:10,267) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,267) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,267) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,267) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,267) = (/
     *0,0,0,0,0,-1,1,0,0,0/)
      NU(51:53,267) = (/
     *0,0,0/)
      NU(1:10,268) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(11:20,268) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(21:30,268) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,268) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(41:50,268) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(51:53,268) = (/
     *0,0,0/)
      NU(1:10,269) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,269) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,269) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,269) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(41:50,269) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(51:53,269) = (/
     *0,0,0/)
      NU(1:10,270) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,270) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,270) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,270) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,270) = (/
     *0,0,0,-1,0,1,0,0,0,0/)
      NU(51:53,270) = (/
     *0,0,0/)
      NU(1:10,271) = (/
     *0,-1,0,0,1,0,0,0,0,0/)
      NU(11:20,271) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,271) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,271) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,271) = (/
     *1,0,0,-1,0,0,0,0,0,0/)
      NU(51:53,271) = (/
     *0,0,0/)
      NU(1:10,272) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,272) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,272) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,272) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(41:50,272) = (/
     *0,0,0,-1,0,0,0,0,0,0/)
      NU(51:53,272) = (/
     *0,0,0/)
      NU(1:10,273) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,273) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,273) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,273) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,273) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(51:53,273) = (/
     *0,0,0/)
      NU(1:10,274) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,274) = (/
     *0,0,0,0,1,0,0,0,0,0/)
      NU(21:30,274) = (/
     *0,0,0,0,0,0,0,-1,0,0/)
      NU(31:40,274) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(41:50,274) = (/
     *0,0,0,1,0,0,0,0,0,0/)
      NU(51:53,274) = (/
     *0,0,0/)
      NU(1:10,275) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(11:20,275) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,275) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,275) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(41:50,275) = (/
     *0,1,0,0,0,0,0,0,0,0/)
      NU(51:53,275) = (/
     *0,0,0/)
      NU(1:10,276) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(11:20,276) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,276) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,276) = (/
     *-1,0,0,0,0,0,0,0,0,0/)
      NU(41:50,276) = (/
     *1,0,0,0,0,0,0,0,0,0/)
      NU(51:53,276) = (/
     *0,0,0/)
      NU(1:10,277) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,277) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,277) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,277) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(41:50,277) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,277) = (/
     *0,0,0/)
      NU(1:10,278) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,278) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,278) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,278) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(41:50,278) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,278) = (/
     *0,0,0/)
      NU(1:10,279) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,279) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,279) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,279) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(41:50,279) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,279) = (/
     *0,0,0/)
      NU(1:10,280) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,280) = (/
     *0,0,0,0,1,-1,0,0,0,0/)
      NU(21:30,280) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,280) = (/
     *0,-1,0,0,0,0,0,0,1,0/)
      NU(41:50,280) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,280) = (/
     *0,0,0/)
      NU(1:10,281) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,281) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,281) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,281) = (/
     *0,0,0,0,0,1,-1,0,0,-1/)
      NU(41:50,281) = (/
     *0,0,0,0,0,0,1,0,0,0/)
      NU(51:53,281) = (/
     *0,0,0/)
      NU(1:10,282) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,282) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(21:30,282) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,282) = (/
     *0,0,0,0,0,0,-1,1,0,0/)
      NU(41:50,282) = (/
     *0,0,0,0,0,0,-1,0,0,0/)
      NU(51:53,282) = (/
     *0,0,0/)
      NU(1:10,283) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,283) = (/
     *0,0,0,0,1,-1,0,0,0,0/)
      NU(21:30,283) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,283) = (/
     *-1,0,0,0,0,1,0,0,0,0/)
      NU(41:50,283) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,283) = (/
     *0,0,0/)
      NU(1:10,284) = (/
     *1,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,284) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(21:30,284) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,284) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,284) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,284) = (/
     *0,0,0/)
      NU(1:10,285) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,285) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,285) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(31:40,285) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,285) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,285) = (/
     *0,1,0/)
      NU(1:10,286) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,286) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,286) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(31:40,286) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,286) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,286) = (/
     *0,0,1/)
      NU(1:10,287) = (/
     *0,0,0,1,-1,1,-1,0,0,0/)
      NU(11:20,287) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,287) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,287) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,287) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,287) = (/
     *0,0,0/)
      NU(1:10,288) = (/
     *1,0,0,0,-1,0,0,0,0,0/)
      NU(11:20,288) = (/
     *0,0,-1,0,0,0,0,1,0,0/)
      NU(21:30,288) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,288) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,288) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,288) = (/
     *0,0,0/)
      NU(1:10,289) = (/
     *-1,0,0,0,0,0,0,0,0,-1/)
      NU(11:20,289) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(21:30,289) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,289) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,289) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,289) = (/
     *0,0,0/)
      NU(1:10,290) = (/
     *0,2,0,-1,0,0,0,0,0,0/)
      NU(11:20,290) = (/
     *-1,0,0,0,0,1,0,0,0,0/)
      NU(21:30,290) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,290) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,290) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,290) = (/
     *0,0,0/)
      NU(1:10,291) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(11:20,291) = (/
     *-1,0,0,0,0,0,0,1,0,0/)
      NU(21:30,291) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,291) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,291) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,291) = (/
     *0,0,0/)
      NU(1:10,292) = (/
     *0,2,0,0,0,0,0,0,0,0/)
      NU(11:20,292) = (/
     *-2,0,0,0,0,0,0,0,0,0/)
      NU(21:30,292) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(31:40,292) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,292) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,292) = (/
     *0,0,0/)
      NU(1:10,293) = (/
     *1,0,0,0,0,-1,0,0,0,0/)
      NU(11:20,293) = (/
     *0,-1,0,0,0,0,0,1,0,0/)
      NU(21:30,293) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,293) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,293) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,293) = (/
     *0,0,0/)
      NU(1:10,294) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(11:20,294) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,294) = (/
     *0,0,0,-1,0,0,0,0,0,0/)
      NU(31:40,294) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,294) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,294) = (/
     *0,1,0/)
      NU(1:10,295) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,295) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,295) = (/
     *0,0,1,-1,0,0,0,0,0,0/)
      NU(31:40,295) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,295) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,295) = (/
     *0,0,0/)
      NU(1:10,296) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,296) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,296) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,296) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,296) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,296) = (/
     *0,1,-1/)
      NU(1:10,297) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,297) = (/
     *0,0,1,0,1,0,0,0,0,0/)
      NU(21:30,297) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,297) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,297) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,297) = (/
     *0,0,-1/)
      NU(1:10,298) = (/
     *0,0,0,-1,0,0,1,0,0,0/)
      NU(11:20,298) = (/
     *0,0,1,0,1,0,0,0,0,0/)
      NU(21:30,298) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,298) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,298) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,298) = (/
     *0,0,-1/)
      NU(1:10,299) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,299) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,299) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,299) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,299) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,299) = (/
     *0,1,-1/)
      NU(1:10,300) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,300) = (/
     *0,0,1,0,1,0,0,0,0,0/)
      NU(21:30,300) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,300) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,300) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,300) = (/
     *0,0,-1/)
      NU(1:10,301) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,301) = (/
     *0,0,1,0,1,0,0,0,0,0/)
      NU(21:30,301) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,301) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,301) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,301) = (/
     *0,0,-1/)
      NU(1:10,302) = (/
     *0,0,0,0,0,0,-1,1,0,0/)
      NU(11:20,302) = (/
     *0,0,1,0,1,0,0,0,0,0/)
      NU(21:30,302) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,302) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,302) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,302) = (/
     *0,0,-1/)
      NU(1:10,303) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,303) = (/
     *0,0,0,1,1,0,0,0,0,0/)
      NU(21:30,303) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,303) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,303) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,303) = (/
     *0,0,-1/)
      NU(1:10,304) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,304) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,304) = (/
     *0,0,0,0,0,0,0,0,-1,0/)
      NU(31:40,304) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,304) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,304) = (/
     *0,1,0/)
      NU(1:10,305) = (/
     *0,1,-1,0,0,0,0,0,0,0/)
      NU(11:20,305) = (/
     *1,0,0,0,0,1,0,0,0,0/)
      NU(21:30,305) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,305) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,305) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,305) = (/
     *0,-1,0/)
      NU(1:10,306) = (/
     *0,0,0,-1,1,0,0,0,0,0/)
      NU(11:20,306) = (/
     *0,0,0,0,1,0,0,1,0,0/)
      NU(21:30,306) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,306) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,306) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,306) = (/
     *0,-1,0/)
      NU(1:10,307) = (/
     *0,0,0,-1,1,0,0,0,0,0/)
      NU(11:20,307) = (/
     *0,0,0,0,0,0,2,0,0,0/)
      NU(21:30,307) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,307) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,307) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,307) = (/
     *0,-1,0/)
      NU(1:10,308) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,308) = (/
     *0,0,1,0,0,0,1,0,0,0/)
      NU(21:30,308) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,308) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,308) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,308) = (/
     *0,-1,0/)
      NU(1:10,309) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,309) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,309) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(31:40,309) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,309) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,309) = (/
     *0,-1,0/)
      NU(1:10,310) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,310) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,310) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(31:40,310) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,310) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,310) = (/
     *0,-1,0/)
      NU(1:10,311) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(11:20,311) = (/
     *0,0,0,0,0,0,1,0,1,0/)
      NU(21:30,311) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,311) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,311) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,311) = (/
     *0,-1,0/)
      NU(1:10,312) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,312) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,312) = (/
     *0,0,0,0,0,-1,0,0,0,0/)
      NU(31:40,312) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,312) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(51:53,312) = (/
     *1,0,0/)
      NU(1:10,313) = (/
     *0,0,-1,0,1,0,0,0,0,0/)
      NU(11:20,313) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,313) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,313) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,313) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(51:53,313) = (/
     *-1,0,0/)
      NU(1:10,314) = (/
     *1,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,314) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,314) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,314) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,314) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(51:53,314) = (/
     *-1,0,0/)
      NU(1:10,315) = (/
     *0,0,0,0,-1,1,0,0,0,0/)
      NU(11:20,315) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,315) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,315) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,315) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(51:53,315) = (/
     *-1,0,0/)
      NU(1:10,316) = (/
     *0,0,0,0,0,0,1,-1,0,0/)
      NU(11:20,316) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,316) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,316) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,316) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(51:53,316) = (/
     *1,0,0/)
      NU(1:10,317) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,317) = (/
     *0,0,-1,1,0,0,0,0,0,0/)
      NU(21:30,317) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,317) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,317) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(51:53,317) = (/
     *-1,0,0/)
      NU(1:10,318) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,318) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,318) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(31:40,318) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,318) = (/
     *0,0,0,0,0,0,0,0,0,1/)
      NU(51:53,318) = (/
     *0,0,0/)
      NU(1:10,319) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(11:20,319) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(21:30,319) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(31:40,319) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,319) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(51:53,319) = (/
     *0,0,0/)
      NU(1:10,320) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,320) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,320) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,320) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,320) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(51:53,320) = (/
     *1,0,0/)
      NU(1:10,321) = (/
     *0,-1,0,0,0,0,0,0,0,0/)
      NU(11:20,321) = (/
     *0,0,1,0,0,0,0,0,0,0/)
      NU(21:30,321) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(31:40,321) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,321) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(51:53,321) = (/
     *0,0,0/)
      NU(1:10,322) = (/
     *0,0,0,0,-1,0,0,0,0,0/)
      NU(11:20,322) = (/
     *0,0,0,0,0,0,0,0,1,0/)
      NU(21:30,322) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(31:40,322) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,322) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(51:53,322) = (/
     *0,0,0/)
      NU(1:10,323) = (/
     *0,0,0,1,0,0,-1,0,0,0/)
      NU(11:20,323) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(21:30,323) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(31:40,323) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,323) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(51:53,323) = (/
     *1,0,0/)
      NU(1:10,324) = (/
     *0,0,0,0,1,0,-1,0,0,0/)
      NU(11:20,324) = (/
     *0,0,0,0,0,0,0,1,0,0/)
      NU(21:30,324) = (/
     *0,0,0,0,0,1,0,0,0,0/)
      NU(31:40,324) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,324) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(51:53,324) = (/
     *0,0,0/)
      NU(1:10,325) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(11:20,325) = (/
     *0,0,-1,0,0,0,0,0,0,0/)
      NU(21:30,325) = (/
     *0,0,0,0,0,2,0,0,0,0/)
      NU(31:40,325) = (/
     *0,0,0,0,0,0,0,0,0,0/)
      NU(41:50,325) = (/
     *0,0,0,0,0,0,0,0,0,-1/)
      NU(51:53,325) = (/
     *0,0,0/)
      RETURN
      END

      SUBROUTINE GETKFKR (KK,II,T,Y,FWDK,REVK)
      INTEGER KK II
      DOUBLE PRECISION RU, PATM, T, d, sumY, kf, Kc, Z, X, k0, kinf
      DOUBLE PRECISION Pr, F, Fcent, a, T1, T2, T3, c, n
      DOUBLE PRECISION, dimension(KK) :: G, Y
      DOUBLE PRECISION, dimension(II) ::  KFARRAY, FWDK, REVK
      CALL GETRU (RU)
      CALL GETPATM (PATM)
      CALL GETG(KK,T,G)
      CALL GETKFARRAY (II, T, KFARRAY)
      PATM=PATM/RU
      G=G/RU
      sumY=SUM(Y)
      d = 0.14

      Z=sumY
      Z=Z+Y(1)*(1.4)
      Z=Z+Y(6)*(14.4)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.75)
      Z=Z+Y(16)*(2.6)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.17)
      kf=KFARRAY(1)*Z
      FWDK(1)=kf*Y(3)*Y(3)
      X=0
      X=X+g(3)*(-2)
      X=X+g(4)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(1)=kf/Kc*Y(4)

      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      kf=KFARRAY(2)*Z
      FWDK(2)=kf*Y(2)*Y(3)
      X=0
      X=X+g(2)*(-1)
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(2)=kf/Kc*Y(5)

      kf=KFARRAY(3)
      FWDK(3)=kf*Y(1)*Y(3)
      X=0
      X=X+g(1)*(-1)
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      Kc=EXP(-X/T)
      REVK(3)=kf/Kc*Y(2)*Y(5)

      kf=KFARRAY(4)
      FWDK(4)=kf*Y(3)*Y(7)
      X=0
      X=X+g(3)*(-1)
      X=X+g(4)*(1)
      X=X+g(5)*(1)
      X=X+g(7)*(-1)
      Kc=EXP(-X/T)
      REVK(4)=kf/Kc*Y(4)*Y(5)

      kf=KFARRAY(5)
      FWDK(5)=kf*Y(3)*Y(8)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(7)*(1)
      X=X+g(8)*(-1)
      Kc=EXP(-X/T)
      REVK(5)=kf/Kc*Y(5)*Y(7)

      kf=KFARRAY(6)
      FWDK(6)=kf*Y(3)*Y(10)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(10)*(-1)
      X=X+g(15)*(1)
      Kc=EXP(-X/T)
      REVK(6)=kf/Kc*Y(2)*Y(15)

      kf=KFARRAY(7)
      FWDK(7)=kf*Y(3)*Y(11)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(11)*(-1)
      X=X+g(17)*(1)
      Kc=EXP(-X/T)
      REVK(7)=kf/Kc*Y(2)*Y(17)

      kf=KFARRAY(8)
      FWDK(8)=kf*Y(3)*Y(12)
      X=0
      X=X+g(1)*(1)
      X=X+g(3)*(-1)
      X=X+g(12)*(-1)
      X=X+g(15)*(1)
      Kc=EXP(-X/T)
      REVK(8)=kf/Kc*Y(1)*Y(15)

      kf=KFARRAY(9)
      FWDK(9)=kf*Y(3)*Y(12)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(12)*(-1)
      X=X+g(17)*(1)
      Kc=EXP(-X/T)
      REVK(9)=kf/Kc*Y(2)*Y(17)

      kf=KFARRAY(10)
      FWDK(10)=kf*Y(3)*Y(13)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(13)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)
      REVK(10)=kf/Kc*Y(2)*Y(18)

      kf=KFARRAY(11)
      FWDK(11)=kf*Y(3)*Y(14)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(13)*(1)
      X=X+g(14)*(-1)
      Kc=EXP(-X/T)
      REVK(11)=kf/Kc*Y(5)*Y(13)

C     lindemann approach
      k0=(T**(.000))*EXP(3.403128e+01 - (1.509650e+03)/T)
      kinf=KFARRAY(12)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(4)*(5)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(2.5)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.5)
      Pr=k0*Z/kinf
      F=1
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(12)=kf*Y(3)*Y(15)
      X=0
      X=X+g(3)*(-1)
      X=X+g(15)*(-1)
      X=X+g(16)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(12)=kf/Kc*Y(16)

      kf=KFARRAY(13)
      FWDK(13)=kf*Y(3)*Y(17)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(15)*(1)
      X=X+g(17)*(-1)
      Kc=EXP(-X/T)
      REVK(13)=kf/Kc*Y(5)*Y(15)

      kf=KFARRAY(14)
      FWDK(14)=kf*Y(3)*Y(17)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(16)*(1)
      X=X+g(17)*(-1)
      Kc=EXP(-X/T)
      REVK(14)=kf/Kc*Y(2)*Y(16)

      kf=KFARRAY(15)
      FWDK(15)=kf*Y(3)*Y(18)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(17)*(1)
      X=X+g(18)*(-1)
      Kc=EXP(-X/T)
      REVK(15)=kf/Kc*Y(5)*Y(17)

      kf=KFARRAY(16)
      FWDK(16)=kf*Y(3)*Y(19)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(18)*(1)
      X=X+g(19)*(-1)
      Kc=EXP(-X/T)
      REVK(16)=kf/Kc*Y(5)*Y(18)

      kf=KFARRAY(17)
      FWDK(17)=kf*Y(3)*Y(20)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(18)*(1)
      X=X+g(20)*(-1)
      Kc=EXP(-X/T)
      REVK(17)=kf/Kc*Y(5)*Y(18)

      kf=KFARRAY(18)
      FWDK(18)=kf*Y(3)*Y(21)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(19)*(1)
      X=X+g(21)*(-1)
      Kc=EXP(-X/T)
      REVK(18)=kf/Kc*Y(5)*Y(19)

      kf=KFARRAY(19)
      FWDK(19)=kf*Y(3)*Y(21)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(20)*(1)
      X=X+g(21)*(-1)
      Kc=EXP(-X/T)
      REVK(19)=kf/Kc*Y(5)*Y(20)

      kf=KFARRAY(20)
      FWDK(20)=kf*Y(3)*Y(22)
      X=0
      X=X+g(3)*(-1)
      X=X+g(10)*(1)
      X=X+g(15)*(1)
      X=X+g(22)*(-1)
      Kc=EXP(-X/T)
      REVK(20)=kf/Kc*Y(10)*Y(15)

      kf=KFARRAY(21)
      FWDK(21)=kf*Y(3)*Y(23)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(23)*(-1)
      X=X+g(28)*(1)
      Kc=EXP(-X/T)
      REVK(21)=kf/Kc*Y(2)*Y(28)

      kf=KFARRAY(22)
      FWDK(22)=kf*Y(3)*Y(23)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(22)*(1)
      X=X+g(23)*(-1)
      Kc=EXP(-X/T)
      REVK(22)=kf/Kc*Y(5)*Y(22)

      kf=KFARRAY(23)
      FWDK(23)=kf*Y(3)*Y(23)
      X=0
      X=X+g(3)*(-1)
      X=X+g(11)*(1)
      X=X+g(15)*(1)
      X=X+g(23)*(-1)
      Kc=EXP(-X/T)
      REVK(23)=kf/Kc*Y(11)*Y(15)

      kf=KFARRAY(24)
      FWDK(24)=kf*Y(3)*Y(24)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(24)*(-1)
      X=X+g(29)*(1)
      Kc=EXP(-X/T)
      REVK(24)=kf/Kc*Y(2)*Y(29)

      kf=KFARRAY(25)
      FWDK(25)=kf*Y(3)*Y(25)
      X=0
      X=X+g(3)*(-1)
      X=X+g(13)*(1)
      X=X+g(17)*(1)
      X=X+g(25)*(-1)
      Kc=EXP(-X/T)
      REVK(25)=kf/Kc*Y(13)*Y(17)

      kf=KFARRAY(26)
      FWDK(26)=kf*Y(3)*Y(26)
      X=0
      X=X+g(3)*(-1)
      X=X+g(13)*(1)
      X=X+g(18)*(1)
      X=X+g(26)*(-1)
      Kc=EXP(-X/T)
      REVK(26)=kf/Kc*Y(13)*Y(18)

      kf=KFARRAY(27)
      FWDK(27)=kf*Y(3)*Y(27)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(26)*(1)
      X=X+g(27)*(-1)
      Kc=EXP(-X/T)
      REVK(27)=kf/Kc*Y(5)*Y(26)

      kf=KFARRAY(28)
      FWDK(28)=kf*Y(3)*Y(28)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(15)*(2)
      X=X+g(28)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(28)=kf/Kc*Y(2)*Y(15)*Y(15)

      kf=KFARRAY(29)
      FWDK(29)=kf*Y(3)*Y(29)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(28)*(1)
      X=X+g(29)*(-1)
      Kc=EXP(-X/T)
      REVK(29)=kf/Kc*Y(5)*Y(28)

      kf=KFARRAY(30)
      FWDK(30)=kf*Y(3)*Y(29)
      X=0
      X=X+g(3)*(-1)
      X=X+g(11)*(1)
      X=X+g(16)*(1)
      X=X+g(29)*(-1)
      Kc=EXP(-X/T)
      REVK(30)=kf/Kc*Y(11)*Y(16)

      kf=KFARRAY(31)
      FWDK(31)=kf*Y(4)*Y(15)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(15)*(-1)
      X=X+g(16)*(1)
      Kc=EXP(-X/T)
      REVK(31)=kf/Kc*Y(3)*Y(16)

      kf=KFARRAY(32)
      FWDK(32)=kf*Y(4)*Y(18)
      X=0
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      X=X+g(17)*(1)
      X=X+g(18)*(-1)
      Kc=EXP(-X/T)
      REVK(32)=kf/Kc*Y(7)*Y(17)

      Z=sumY
      Z=Z+Y(4)*(-1)
      Z=Z+Y(6)*(-1)
      Z=Z+Y(15)*(-0.25)
      Z=Z+Y(16)*(0.5)
      Z=Z+Y(27)*(0.5)
      Z=Z+Y(48)*(-1)
      Z=Z+Y(49)*(-1)
      kf=KFARRAY(33)*Z
      FWDK(33)=kf*Y(2)*Y(4)
      X=0
      X=X+g(2)*(-1)
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(33)=kf/Kc*Y(7)

      kf=KFARRAY(34)
      FWDK(34)=kf*Y(2)*Y(4)*Y(4)
      X=0
      X=X+g(2)*(-1)
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(34)=kf/Kc*Y(4)*Y(7)

      kf=KFARRAY(35)
      FWDK(35)=kf*Y(2)*Y(4)*Y(6)
      X=0
      X=X+g(2)*(-1)
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(35)=kf/Kc*Y(6)*Y(7)

      kf=KFARRAY(36)
      FWDK(36)=kf*Y(2)*Y(4)*Y(48)
      X=0
      X=X+g(2)*(-1)
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(36)=kf/Kc*Y(7)*Y(48)

      kf=KFARRAY(37)
      FWDK(37)=kf*Y(2)*Y(4)*Y(49)
      X=0
      X=X+g(2)*(-1)
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(37)=kf/Kc*Y(7)*Y(49)

      kf=KFARRAY(38)
      FWDK(38)=kf*Y(2)*Y(4)
      X=0
      X=X+g(2)*(-1)
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(5)*(1)
      Kc=EXP(-X/T)
      REVK(38)=kf/Kc*Y(3)*Y(5)

      Z=sumY
      Z=Z+Y(1)*(-1)
      Z=Z+Y(6)*(-1)
      Z=Z+Y(14)*(1)
      Z=Z+Y(16)*(-1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.37)
      kf=KFARRAY(39)*Z
      FWDK(39)=kf*Y(2)*Y(2)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-2)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(39)=kf/Kc*Y(1)

      kf=KFARRAY(40)
      FWDK(40)=kf*Y(1)*Y(2)*Y(2)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-2)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(40)=kf/Kc*Y(1)*Y(1)

      kf=KFARRAY(41)
      FWDK(41)=kf*Y(2)*Y(2)*Y(6)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-2)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(41)=kf/Kc*Y(1)*Y(6)

      kf=KFARRAY(42)
      FWDK(42)=kf*Y(2)*Y(2)*Y(16)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-2)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(42)=kf/Kc*Y(1)*Y(16)

      Z=sumY
      Z=Z+Y(1)*(-0.27)
      Z=Z+Y(6)*(2.65)
      Z=Z+Y(14)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.62)
      kf=KFARRAY(43)*Z
      FWDK(43)=kf*Y(2)*Y(5)
      X=0
      X=X+g(2)*(-1)
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(43)=kf/Kc*Y(6)

      kf=KFARRAY(44)
      FWDK(44)=kf*Y(2)*Y(7)
      X=0
      X=X+g(2)*(-1)
      X=X+g(3)*(1)
      X=X+g(6)*(1)
      X=X+g(7)*(-1)
      Kc=EXP(-X/T)
      REVK(44)=kf/Kc*Y(3)*Y(6)

      kf=KFARRAY(45)
      FWDK(45)=kf*Y(2)*Y(7)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(4)*(1)
      X=X+g(7)*(-1)
      Kc=EXP(-X/T)
      REVK(45)=kf/Kc*Y(1)*Y(4)

      kf=KFARRAY(46)
      FWDK(46)=kf*Y(2)*Y(7)
      X=0
      X=X+g(2)*(-1)
      X=X+g(5)*(2)
      X=X+g(7)*(-1)
      Kc=EXP(-X/T)
      REVK(46)=kf/Kc*Y(5)*Y(5)

      kf=KFARRAY(47)
      FWDK(47)=kf*Y(2)*Y(8)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(7)*(1)
      X=X+g(8)*(-1)
      Kc=EXP(-X/T)
      REVK(47)=kf/Kc*Y(1)*Y(7)

      kf=KFARRAY(48)
      FWDK(48)=kf*Y(2)*Y(8)
      X=0
      X=X+g(2)*(-1)
      X=X+g(5)*(1)
      X=X+g(6)*(1)
      X=X+g(8)*(-1)
      Kc=EXP(-X/T)
      REVK(48)=kf/Kc*Y(5)*Y(6)

      kf=KFARRAY(49)
      FWDK(49)=kf*Y(2)*Y(10)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(9)*(1)
      X=X+g(10)*(-1)
      Kc=EXP(-X/T)
      REVK(49)=kf/Kc*Y(1)*Y(9)

C     lindemann approach
      k0=(T**(-2.760))*EXP(5.990643e+01 - (8.051467e+02)/T)
      kinf=KFARRAY(50)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.5620
      T3=91.00
      T1=5836.00
      T2=8552.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(50)=kf*Y(2)*Y(11)
      X=0
      X=X+g(2)*(-1)
      X=X+g(11)*(-1)
      X=X+g(13)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(50)=kf/Kc*Y(13)

      kf=KFARRAY(51)
      FWDK(51)=kf*Y(2)*Y(12)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(10)*(1)
      X=X+g(12)*(-1)
      Kc=EXP(-X/T)
      REVK(51)=kf/Kc*Y(1)*Y(10)

C     lindemann approach
      k0=(T**(-4.760))*EXP(7.694848e+01 - (1.227849e+03)/T)
      kinf=KFARRAY(52)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(2)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.7830
      T3=74.00
      T1=2941.00
      T2=6964.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(52)=kf*Y(2)*Y(13)
      X=0
      X=X+g(2)*(-1)
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(52)=kf/Kc*Y(14)

      kf=KFARRAY(53)
      FWDK(53)=kf*Y(2)*Y(14)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(13)*(1)
      X=X+g(14)*(-1)
      Kc=EXP(-X/T)
      REVK(53)=kf/Kc*Y(1)*Y(13)

C     lindemann approach
      k0=(T**(-2.570))*EXP(5.616626e+01 - (2.138671e+02)/T)
      kinf=KFARRAY(54)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.7824
      T3=271.00
      T1=2755.00
      T2=6570.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(54)=kf*Y(2)*Y(17)
      X=0
      X=X+g(2)*(-1)
      X=X+g(17)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(54)=kf/Kc*Y(18)

      kf=KFARRAY(55)
      FWDK(55)=kf*Y(2)*Y(17)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(15)*(1)
      X=X+g(17)*(-1)
      Kc=EXP(-X/T)
      REVK(55)=kf/Kc*Y(1)*Y(15)

C     lindemann approach
      k0=(T**(-4.820))*EXP(7.392174e+01 - (3.286005e+03)/T)
      kinf=KFARRAY(56)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Pr=k0*Z/kinf
C     Troe form
      a=.7187
      T3=103.00
      T1=1291.00
      T2=4160.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(56)=kf*Y(2)*Y(18)
      X=0
      X=X+g(2)*(-1)
      X=X+g(18)*(-1)
      X=X+g(19)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(56)=kf/Kc*Y(19)

C     lindemann approach
      k0=(T**(-4.800))*EXP(6.986601e+01 - (2.797885e+03)/T)
      kinf=KFARRAY(57)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Pr=k0*Z/kinf
C     Troe form
      a=.7580
      T3=94.00
      T1=1555.00
      T2=4200.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(57)=kf*Y(2)*Y(18)
      X=0
      X=X+g(2)*(-1)
      X=X+g(18)*(-1)
      X=X+g(20)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(57)=kf/Kc*Y(20)

      kf=KFARRAY(58)
      FWDK(58)=kf*Y(2)*Y(18)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(17)*(1)
      X=X+g(18)*(-1)
      Kc=EXP(-X/T)
      REVK(58)=kf/Kc*Y(1)*Y(17)

C     lindemann approach
      k0=(T**(-4.650))*EXP(7.285261e+01 - (2.556341e+03)/T)
      kinf=KFARRAY(59)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Pr=k0*Z/kinf
C     Troe form
      a=.600
      T3=100.00
      T1=90000.0
      T2=10000.0
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(59)=kf*Y(2)*Y(19)
      X=0
      X=X+g(2)*(-1)
      X=X+g(19)*(-1)
      X=X+g(21)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(59)=kf/Kc*Y(21)

      kf=KFARRAY(60)
      FWDK(60)=kf*Y(2)*Y(19)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(18)*(1)
      X=X+g(19)*(-1)
      Kc=EXP(-X/T)
      REVK(60)=kf/Kc*Y(1)*Y(18)

      kf=KFARRAY(61)
      FWDK(61)=kf*Y(2)*Y(19)
      X=0
      X=X+g(2)*(-1)
      X=X+g(5)*(1)
      X=X+g(13)*(1)
      X=X+g(19)*(-1)
      Kc=EXP(-X/T)
      REVK(61)=kf/Kc*Y(5)*Y(13)

      kf=KFARRAY(62)
      FWDK(62)=kf*Y(2)*Y(19)
      X=0
      X=X+g(2)*(-1)
      X=X+g(6)*(1)
      X=X+g(12)*(1)
      X=X+g(19)*(-1)
      Kc=EXP(-X/T)
      REVK(62)=kf/Kc*Y(6)*Y(12)

C     lindemann approach
      k0=(T**(-7.440))*EXP(9.594500e+01 - (7.085291e+03)/T)
      kinf=KFARRAY(63)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Pr=k0*Z/kinf
C     Troe form
      a=.700
      T3=100.00
      T1=90000.0
      T2=10000.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(63)=kf*Y(2)*Y(20)
      X=0
      X=X+g(2)*(-1)
      X=X+g(20)*(-1)
      X=X+g(21)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(63)=kf/Kc*Y(21)

      kf=KFARRAY(64)
      FWDK(64)=kf*Y(2)*Y(20)
      X=0
      X=X+g(19)*(1)
      X=X+g(20)*(-1)
      Kc=EXP(-X/T)
      REVK(64)=kf/Kc*Y(2)*Y(19)

      kf=KFARRAY(65)
      FWDK(65)=kf*Y(2)*Y(20)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(18)*(1)
      X=X+g(20)*(-1)
      Kc=EXP(-X/T)
      REVK(65)=kf/Kc*Y(1)*Y(18)

      kf=KFARRAY(66)
      FWDK(66)=kf*Y(2)*Y(20)
      X=0
      X=X+g(2)*(-1)
      X=X+g(5)*(1)
      X=X+g(13)*(1)
      X=X+g(20)*(-1)
      Kc=EXP(-X/T)
      REVK(66)=kf/Kc*Y(5)*Y(13)

      kf=KFARRAY(67)
      FWDK(67)=kf*Y(2)*Y(20)
      X=0
      X=X+g(2)*(-1)
      X=X+g(6)*(1)
      X=X+g(12)*(1)
      X=X+g(20)*(-1)
      Kc=EXP(-X/T)
      REVK(67)=kf/Kc*Y(6)*Y(12)

      kf=KFARRAY(68)
      FWDK(68)=kf*Y(2)*Y(21)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(19)*(1)
      X=X+g(21)*(-1)
      Kc=EXP(-X/T)
      REVK(68)=kf/Kc*Y(1)*Y(19)

      kf=KFARRAY(69)
      FWDK(69)=kf*Y(2)*Y(21)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(20)*(1)
      X=X+g(21)*(-1)
      Kc=EXP(-X/T)
      REVK(69)=kf/Kc*Y(1)*Y(20)

C     lindemann approach
      k0=(T**(-4.800))*EXP(7.730706e+01 - (9.561117e+02)/T)
      kinf=KFARRAY(70)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.6464
      T3=132.00
      T1=1315.00
      T2=5566.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(70)=kf*Y(2)*Y(22)
      X=0
      X=X+g(2)*(-1)
      X=X+g(22)*(-1)
      X=X+g(23)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(70)=kf/Kc*Y(23)

C     lindemann approach
      k0=(T**(-7.270))*EXP(9.343840e+01 - (3.633224e+03)/T)
      kinf=KFARRAY(71)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.7507
      T3=98.50
      T1=1302.00
      T2=4167.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(71)=kf*Y(2)*Y(23)
      X=0
      X=X+g(2)*(-1)
      X=X+g(23)*(-1)
      X=X+g(24)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(71)=kf/Kc*Y(24)

C     lindemann approach
      k0=(T**(-3.860))*EXP(6.941403e+01 - (1.670679e+03)/T)
      kinf=KFARRAY(72)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.7820
      T3=207.50
      T1=2663.00
      T2=6095.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(72)=kf*Y(2)*Y(24)
      X=0
      X=X+g(2)*(-1)
      X=X+g(24)*(-1)
      X=X+g(25)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(72)=kf/Kc*Y(25)

      kf=KFARRAY(73)
      FWDK(73)=kf*Y(2)*Y(24)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(23)*(1)
      X=X+g(24)*(-1)
      Kc=EXP(-X/T)
      REVK(73)=kf/Kc*Y(1)*Y(23)

C     lindemann approach
      k0=(T**(-7.620))*EXP(9.619775e+01 - (3.507420e+03)/T)
      kinf=KFARRAY(74)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.9753
      T3=210.00
      T1=984.00
      T2=4374.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(74)=kf*Y(2)*Y(25)
      X=0
      X=X+g(2)*(-1)
      X=X+g(25)*(-1)
      X=X+g(26)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(74)=kf/Kc*Y(26)

      kf=KFARRAY(75)
      FWDK(75)=kf*Y(2)*Y(25)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(24)*(1)
      X=X+g(25)*(-1)
      Kc=EXP(-X/T)
      REVK(75)=kf/Kc*Y(1)*Y(24)

C     lindemann approach
      k0=(T**(-7.080))*EXP(9.509412e+01 - (3.364003e+03)/T)
      kinf=KFARRAY(76)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.8422
      T3=125.00
      T1=2219.00
      T2=6882.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(76)=kf*Y(2)*Y(26)
      X=0
      X=X+g(2)*(-1)
      X=X+g(26)*(-1)
      X=X+g(27)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(76)=kf/Kc*Y(27)

      kf=KFARRAY(77)
      FWDK(77)=kf*Y(2)*Y(26)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(25)*(1)
      X=X+g(26)*(-1)
      Kc=EXP(-X/T)
      REVK(77)=kf/Kc*Y(1)*Y(25)

      kf=KFARRAY(78)
      FWDK(78)=kf*Y(2)*Y(27)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(26)*(1)
      X=X+g(27)*(-1)
      Kc=EXP(-X/T)
      REVK(78)=kf/Kc*Y(1)*Y(26)

      kf=KFARRAY(79)
      FWDK(79)=kf*Y(2)*Y(28)
      X=0
      X=X+g(2)*(-1)
      X=X+g(12)*(1)
      X=X+g(15)*(1)
      X=X+g(28)*(-1)
      Kc=EXP(-X/T)
      REVK(79)=kf/Kc*Y(12)*Y(15)

      kf=KFARRAY(80)
      FWDK(80)=kf*Y(2)*Y(29)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(28)*(1)
      X=X+g(29)*(-1)
      Kc=EXP(-X/T)
      REVK(80)=kf/Kc*Y(1)*Y(28)

      kf=KFARRAY(81)
      FWDK(81)=kf*Y(2)*Y(29)
      X=0
      X=X+g(2)*(-1)
      X=X+g(13)*(1)
      X=X+g(15)*(1)
      X=X+g(29)*(-1)
      Kc=EXP(-X/T)
      REVK(81)=kf/Kc*Y(13)*Y(15)

      kf=KFARRAY(82)
      FWDK(82)=kf*Y(2)*Y(30)
      X=0
      X=X+g(29)*(1)
      X=X+g(30)*(-1)
      Kc=EXP(-X/T)
      REVK(82)=kf/Kc*Y(2)*Y(29)

C     lindemann approach
      k0=(T**(-3.420))*EXP(6.379314e+01 - (4.244633e+04)/T)
      kinf=KFARRAY(83)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.9320
      T3=197.00
      T1=1540.00
      T2=10300.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(83)=kf*Y(1)*Y(15)
      X=0
      X=X+g(1)*(-1)
      X=X+g(15)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(83)=kf/Kc*Y(18)

      kf=KFARRAY(84)
      FWDK(84)=kf*Y(1)*Y(5)
      X=0
      X=X+g(1)*(-1)
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      Kc=EXP(-X/T)
      REVK(84)=kf/Kc*Y(2)*Y(6)

C     lindemann approach
      k0=(T**(-.900))*EXP(4.227944e+01 - (-8.554683e+02)/T)
      kinf=KFARRAY(85)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.7346
      T3=94.00
      T1=1756.00
      T2=5182.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(85)=kf*Y(5)*Y(5)
      X=0
      X=X+g(5)*(-2)
      X=X+g(8)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(85)=kf/Kc*Y(8)

      kf=KFARRAY(86)
      FWDK(86)=kf*Y(5)*Y(5)
      X=0
      X=X+g(3)*(1)
      X=X+g(5)*(-2)
      X=X+g(6)*(1)
      Kc=EXP(-X/T)
      REVK(86)=kf/Kc*Y(3)*Y(6)

      kf=KFARRAY(87)
      FWDK(87)=kf*Y(5)*Y(7)
      X=0
      X=X+g(4)*(1)
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(7)*(-1)
      Kc=EXP(-X/T)
      REVK(87)=kf/Kc*Y(4)*Y(6)

      kf=KFARRAY(88)
      FWDK(88)=kf*Y(5)*Y(8)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(7)*(1)
      X=X+g(8)*(-1)
      Kc=EXP(-X/T)
      REVK(88)=kf/Kc*Y(6)*Y(7)

      kf=KFARRAY(89)
      FWDK(89)=kf*Y(5)*Y(8)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(7)*(1)
      X=X+g(8)*(-1)
      Kc=EXP(-X/T)
      REVK(89)=kf/Kc*Y(6)*Y(7)

      kf=KFARRAY(90)
      FWDK(90)=kf*Y(5)*Y(9)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(9)*(-1)
      X=X+g(15)*(1)
      Kc=EXP(-X/T)
      REVK(90)=kf/Kc*Y(2)*Y(15)

      kf=KFARRAY(91)
      FWDK(91)=kf*Y(5)*Y(10)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(10)*(-1)
      X=X+g(17)*(1)
      Kc=EXP(-X/T)
      REVK(91)=kf/Kc*Y(2)*Y(17)

      kf=KFARRAY(92)
      FWDK(92)=kf*Y(5)*Y(11)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(11)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)
      REVK(92)=kf/Kc*Y(2)*Y(18)

      kf=KFARRAY(93)
      FWDK(93)=kf*Y(5)*Y(11)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(10)*(1)
      X=X+g(11)*(-1)
      Kc=EXP(-X/T)
      REVK(93)=kf/Kc*Y(6)*Y(10)

      kf=KFARRAY(94)
      FWDK(94)=kf*Y(5)*Y(12)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(12)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)
      REVK(94)=kf/Kc*Y(2)*Y(18)

C     lindemann approach
      k0=(T**(-5.920))*EXP(8.427936e+01 - (1.580100e+03)/T)
      kinf=KFARRAY(95)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Pr=k0*Z/kinf
C     Troe form
      a=.4120
      T3=195.0
      T1=5900.00
      T2=6394.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(95)=kf*Y(5)*Y(13)
      X=0
      X=X+g(5)*(-1)
      X=X+g(13)*(-1)
      X=X+g(21)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(95)=kf/Kc*Y(21)

      kf=KFARRAY(96)
      FWDK(96)=kf*Y(5)*Y(13)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(11)*(1)
      X=X+g(13)*(-1)
      Kc=EXP(-X/T)
      REVK(96)=kf/Kc*Y(6)*Y(11)

      kf=KFARRAY(97)
      FWDK(97)=kf*Y(5)*Y(13)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(12)*(1)
      X=X+g(13)*(-1)
      Kc=EXP(-X/T)
      REVK(97)=kf/Kc*Y(6)*Y(12)

      kf=KFARRAY(98)
      FWDK(98)=kf*Y(5)*Y(14)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(13)*(1)
      X=X+g(14)*(-1)
      Kc=EXP(-X/T)
      REVK(98)=kf/Kc*Y(6)*Y(13)

      kf=KFARRAY(99)
      FWDK(99)=kf*Y(5)*Y(15)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(15)*(-1)
      X=X+g(16)*(1)
      Kc=EXP(-X/T)
      REVK(99)=kf/Kc*Y(2)*Y(16)

      kf=KFARRAY(100)
      FWDK(100)=kf*Y(5)*Y(17)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(15)*(1)
      X=X+g(17)*(-1)
      Kc=EXP(-X/T)
      REVK(100)=kf/Kc*Y(6)*Y(15)

      kf=KFARRAY(101)
      FWDK(101)=kf*Y(5)*Y(18)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(17)*(1)
      X=X+g(18)*(-1)
      Kc=EXP(-X/T)
      REVK(101)=kf/Kc*Y(6)*Y(17)

      kf=KFARRAY(102)
      FWDK(102)=kf*Y(5)*Y(19)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(18)*(1)
      X=X+g(19)*(-1)
      Kc=EXP(-X/T)
      REVK(102)=kf/Kc*Y(6)*Y(18)

      kf=KFARRAY(103)
      FWDK(103)=kf*Y(5)*Y(20)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(18)*(1)
      X=X+g(20)*(-1)
      Kc=EXP(-X/T)
      REVK(103)=kf/Kc*Y(6)*Y(18)

      kf=KFARRAY(104)
      FWDK(104)=kf*Y(5)*Y(21)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(19)*(1)
      X=X+g(21)*(-1)
      Kc=EXP(-X/T)
      REVK(104)=kf/Kc*Y(6)*Y(19)

      kf=KFARRAY(105)
      FWDK(105)=kf*Y(5)*Y(21)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(20)*(1)
      X=X+g(21)*(-1)
      Kc=EXP(-X/T)
      REVK(105)=kf/Kc*Y(6)*Y(20)

      kf=KFARRAY(106)
      FWDK(106)=kf*Y(5)*Y(22)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(22)*(-1)
      X=X+g(28)*(1)
      Kc=EXP(-X/T)
      REVK(106)=kf/Kc*Y(2)*Y(28)

      kf=KFARRAY(107)
      FWDK(107)=kf*Y(5)*Y(23)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(23)*(-1)
      X=X+g(29)*(1)
      Kc=EXP(-X/T)
      REVK(107)=kf/Kc*Y(2)*Y(29)

      kf=KFARRAY(108)
      FWDK(108)=kf*Y(5)*Y(23)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(23)*(-1)
      X=X+g(30)*(1)
      Kc=EXP(-X/T)
      REVK(108)=kf/Kc*Y(2)*Y(30)

      kf=KFARRAY(109)
      FWDK(109)=kf*Y(5)*Y(23)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(22)*(1)
      X=X+g(23)*(-1)
      Kc=EXP(-X/T)
      REVK(109)=kf/Kc*Y(6)*Y(22)

      kf=KFARRAY(110)
      FWDK(110)=kf*Y(5)*Y(23)
      X=0
      X=X+g(5)*(-1)
      X=X+g(13)*(1)
      X=X+g(15)*(1)
      X=X+g(23)*(-1)
      Kc=EXP(-X/T)
      REVK(110)=kf/Kc*Y(13)*Y(15)

      kf=KFARRAY(111)
      FWDK(111)=kf*Y(5)*Y(24)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(23)*(1)
      X=X+g(24)*(-1)
      Kc=EXP(-X/T)
      REVK(111)=kf/Kc*Y(6)*Y(23)

      kf=KFARRAY(112)
      FWDK(112)=kf*Y(5)*Y(25)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(24)*(1)
      X=X+g(25)*(-1)
      Kc=EXP(-X/T)
      REVK(112)=kf/Kc*Y(6)*Y(24)

      kf=KFARRAY(113)
      FWDK(113)=kf*Y(5)*Y(27)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(26)*(1)
      X=X+g(27)*(-1)
      Kc=EXP(-X/T)
      REVK(113)=kf/Kc*Y(6)*Y(26)

      kf=KFARRAY(114)
      FWDK(114)=kf*Y(5)*Y(29)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(28)*(1)
      X=X+g(29)*(-1)
      Kc=EXP(-X/T)
      REVK(114)=kf/Kc*Y(6)*Y(28)

      kf=KFARRAY(115)
      FWDK(115)=kf*Y(7)*Y(7)
      X=0
      X=X+g(4)*(1)
      X=X+g(7)*(-2)
      X=X+g(8)*(1)
      Kc=EXP(-X/T)
      REVK(115)=kf/Kc*Y(4)*Y(8)

      kf=KFARRAY(116)
      FWDK(116)=kf*Y(7)*Y(7)
      X=0
      X=X+g(4)*(1)
      X=X+g(7)*(-2)
      X=X+g(8)*(1)
      Kc=EXP(-X/T)
      REVK(116)=kf/Kc*Y(4)*Y(8)

      kf=KFARRAY(117)
      FWDK(117)=kf*Y(7)*Y(11)
      X=0
      X=X+g(5)*(1)
      X=X+g(7)*(-1)
      X=X+g(11)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)
      REVK(117)=kf/Kc*Y(5)*Y(18)

      kf=KFARRAY(118)
      FWDK(118)=kf*Y(7)*Y(13)
      X=0
      X=X+g(4)*(1)
      X=X+g(7)*(-1)
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      Kc=EXP(-X/T)
      REVK(118)=kf/Kc*Y(4)*Y(14)

      kf=KFARRAY(119)
      FWDK(119)=kf*Y(7)*Y(13)
      X=0
      X=X+g(5)*(1)
      X=X+g(7)*(-1)
      X=X+g(13)*(-1)
      X=X+g(20)*(1)
      Kc=EXP(-X/T)
      REVK(119)=kf/Kc*Y(5)*Y(20)

      kf=KFARRAY(120)
      FWDK(120)=kf*Y(7)*Y(15)
      X=0
      X=X+g(5)*(1)
      X=X+g(7)*(-1)
      X=X+g(15)*(-1)
      X=X+g(16)*(1)
      Kc=EXP(-X/T)
      REVK(120)=kf/Kc*Y(5)*Y(16)

      kf=KFARRAY(121)
      FWDK(121)=kf*Y(7)*Y(18)
      X=0
      X=X+g(7)*(-1)
      X=X+g(8)*(1)
      X=X+g(17)*(1)
      X=X+g(18)*(-1)
      Kc=EXP(-X/T)
      REVK(121)=kf/Kc*Y(8)*Y(17)

      kf=KFARRAY(122)
      FWDK(122)=kf*Y(4)*Y(9)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(9)*(-1)
      X=X+g(15)*(1)
      Kc=EXP(-X/T)
      REVK(122)=kf/Kc*Y(3)*Y(15)

      kf=KFARRAY(123)
      FWDK(123)=kf*Y(9)*Y(11)
      X=0
      X=X+g(2)*(1)
      X=X+g(9)*(-1)
      X=X+g(11)*(-1)
      X=X+g(22)*(1)
      Kc=EXP(-X/T)
      REVK(123)=kf/Kc*Y(2)*Y(22)

      kf=KFARRAY(124)
      FWDK(124)=kf*Y(9)*Y(13)
      X=0
      X=X+g(2)*(1)
      X=X+g(9)*(-1)
      X=X+g(13)*(-1)
      X=X+g(23)*(1)
      Kc=EXP(-X/T)
      REVK(124)=kf/Kc*Y(2)*Y(23)

      kf=KFARRAY(125)
      FWDK(125)=kf*Y(4)*Y(10)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(10)*(-1)
      X=X+g(17)*(1)
      Kc=EXP(-X/T)
      REVK(125)=kf/Kc*Y(3)*Y(17)

      kf=KFARRAY(126)
      FWDK(126)=kf*Y(1)*Y(10)
      X=0
      X=X+g(1)*(-1)
      X=X+g(2)*(1)
      X=X+g(10)*(-1)
      X=X+g(11)*(1)
      Kc=EXP(-X/T)
      REVK(126)=kf/Kc*Y(2)*Y(11)

      kf=KFARRAY(127)
      FWDK(127)=kf*Y(6)*Y(10)
      X=0
      X=X+g(2)*(1)
      X=X+g(6)*(-1)
      X=X+g(10)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)
      REVK(127)=kf/Kc*Y(2)*Y(18)

      kf=KFARRAY(128)
      FWDK(128)=kf*Y(10)*Y(11)
      X=0
      X=X+g(2)*(1)
      X=X+g(10)*(-1)
      X=X+g(11)*(-1)
      X=X+g(23)*(1)
      Kc=EXP(-X/T)
      REVK(128)=kf/Kc*Y(2)*Y(23)

      kf=KFARRAY(129)
      FWDK(129)=kf*Y(10)*Y(13)
      X=0
      X=X+g(2)*(1)
      X=X+g(10)*(-1)
      X=X+g(13)*(-1)
      X=X+g(24)*(1)
      Kc=EXP(-X/T)
      REVK(129)=kf/Kc*Y(2)*Y(24)

      kf=KFARRAY(130)
      FWDK(130)=kf*Y(10)*Y(14)
      X=0
      X=X+g(2)*(1)
      X=X+g(10)*(-1)
      X=X+g(14)*(-1)
      X=X+g(25)*(1)
      Kc=EXP(-X/T)
      REVK(130)=kf/Kc*Y(2)*Y(25)

C     lindemann approach
      k0=(T**(-3.740))*EXP(6.546192e+01 - (9.742275e+02)/T)
      kinf=KFARRAY(131)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.5757
      T3=237.00
      T1=1652.00
      T2=5069.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(131)=kf*Y(10)*Y(15)
      X=0
      X=X+g(10)*(-1)
      X=X+g(15)*(-1)
      X=X+g(28)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(131)=kf/Kc*Y(28)

      kf=KFARRAY(132)
      FWDK(132)=kf*Y(10)*Y(16)
      X=0
      X=X+g(10)*(-1)
      X=X+g(15)*(1)
      X=X+g(16)*(-1)
      X=X+g(17)*(1)
      Kc=EXP(-X/T)
      REVK(132)=kf/Kc*Y(15)*Y(17)

      kf=KFARRAY(133)
      FWDK(133)=kf*Y(10)*Y(18)
      X=0
      X=X+g(2)*(1)
      X=X+g(10)*(-1)
      X=X+g(18)*(-1)
      X=X+g(29)*(1)
      Kc=EXP(-X/T)
      REVK(133)=kf/Kc*Y(2)*Y(29)

      kf=KFARRAY(134)
      FWDK(134)=kf*Y(10)*Y(28)
      X=0
      X=X+g(10)*(-1)
      X=X+g(15)*(1)
      X=X+g(23)*(1)
      X=X+g(28)*(-1)
      Kc=EXP(-X/T)
      REVK(134)=kf/Kc*Y(15)*Y(23)

      kf=KFARRAY(135)
      FWDK(135)=kf*Y(4)*Y(11)

      kf=KFARRAY(136)
      FWDK(136)=kf*Y(1)*Y(11)
      X=0
      X=X+g(1)*(-1)
      X=X+g(2)*(1)
      X=X+g(11)*(-1)
      X=X+g(13)*(1)
      Kc=EXP(-X/T)
      REVK(136)=kf/Kc*Y(2)*Y(13)

      kf=KFARRAY(137)
      FWDK(137)=kf*Y(11)*Y(11)
      X=0
      X=X+g(1)*(1)
      X=X+g(11)*(-2)
      X=X+g(23)*(1)
      Kc=EXP(-X/T)
      REVK(137)=kf/Kc*Y(1)*Y(23)

      kf=KFARRAY(138)
      FWDK(138)=kf*Y(11)*Y(13)
      X=0
      X=X+g(2)*(1)
      X=X+g(11)*(-1)
      X=X+g(13)*(-1)
      X=X+g(25)*(1)
      Kc=EXP(-X/T)
      REVK(138)=kf/Kc*Y(2)*Y(25)

      kf=KFARRAY(139)
      FWDK(139)=kf*Y(11)*Y(14)
      X=0
      X=X+g(11)*(-1)
      X=X+g(13)*(2)
      X=X+g(14)*(-1)
      Kc=EXP(-X/T)
      REVK(139)=kf/Kc*Y(13)*Y(13)

C     lindemann approach
      k0=(T**(-5.110))*EXP(7.697485e+01 - (3.570322e+03)/T)
      kinf=KFARRAY(140)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.5907
      T3=275.00
      T1=1226.00
      T2=5185.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(140)=kf*Y(11)*Y(15)
      X=0
      X=X+g(11)*(-1)
      X=X+g(15)*(-1)
      X=X+g(29)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(140)=kf/Kc*Y(29)

      kf=KFARRAY(141)
      FWDK(141)=kf*Y(11)*Y(28)
      X=0
      X=X+g(11)*(-1)
      X=X+g(15)*(1)
      X=X+g(24)*(1)
      X=X+g(28)*(-1)
      Kc=EXP(-X/T)
      REVK(141)=kf/Kc*Y(15)*Y(24)

      kf=KFARRAY(142)
      FWDK(142)=kf*Y(12)*Y(48)
      X=0
      X=X+g(11)*(1)
      X=X+g(12)*(-1)
      Kc=EXP(-X/T)
      REVK(142)=kf/Kc*Y(11)*Y(48)

      kf=KFARRAY(143)
      FWDK(143)=kf*Y(12)*Y(49)
      X=0
      X=X+g(11)*(1)
      X=X+g(12)*(-1)
      Kc=EXP(-X/T)
      REVK(143)=kf/Kc*Y(11)*Y(49)

      kf=KFARRAY(144)
      FWDK(144)=kf*Y(4)*Y(12)
      X=0
      X=X+g(2)*(1)
      X=X+g(4)*(-1)
      X=X+g(5)*(1)
      X=X+g(12)*(-1)
      X=X+g(15)*(1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(144)=kf/Kc*Y(2)*Y(5)*Y(15)

      kf=KFARRAY(145)
      FWDK(145)=kf*Y(4)*Y(12)
      X=0
      X=X+g(4)*(-1)
      X=X+g(6)*(1)
      X=X+g(12)*(-1)
      X=X+g(15)*(1)
      Kc=EXP(-X/T)
      REVK(145)=kf/Kc*Y(6)*Y(15)

      kf=KFARRAY(146)
      FWDK(146)=kf*Y(1)*Y(12)
      X=0
      X=X+g(1)*(-1)
      X=X+g(2)*(1)
      X=X+g(12)*(-1)
      X=X+g(13)*(1)
      Kc=EXP(-X/T)
      REVK(146)=kf/Kc*Y(2)*Y(13)

C     lindemann approach
      k0=(T**(-6.360))*EXP(8.812951e+01 - (2.536212e+03)/T)
      kinf=KFARRAY(147)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Pr=k0*Z/kinf
C     Troe form
      a=.6027
      T3=208.00
      T1=3922.00
      T2=10180.0
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(147)=kf*Y(6)*Y(12)
      X=0
      X=X+g(6)*(-1)
      X=X+g(12)*(-1)
      X=X+g(21)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(147)=kf/Kc*Y(21)

      kf=KFARRAY(148)
      FWDK(148)=kf*Y(6)*Y(12)
      X=0
      X=X+g(11)*(1)
      X=X+g(12)*(-1)
      Kc=EXP(-X/T)
      REVK(148)=kf/Kc*Y(6)*Y(11)

      kf=KFARRAY(149)
      FWDK(149)=kf*Y(12)*Y(13)
      X=0
      X=X+g(2)*(1)
      X=X+g(12)*(-1)
      X=X+g(13)*(-1)
      X=X+g(25)*(1)
      Kc=EXP(-X/T)
      REVK(149)=kf/Kc*Y(2)*Y(25)

      kf=KFARRAY(150)
      FWDK(150)=kf*Y(12)*Y(14)
      X=0
      X=X+g(12)*(-1)
      X=X+g(13)*(2)
      X=X+g(14)*(-1)
      Kc=EXP(-X/T)
      REVK(150)=kf/Kc*Y(13)*Y(13)

      kf=KFARRAY(151)
      FWDK(151)=kf*Y(12)*Y(15)
      X=0
      X=X+g(11)*(1)
      X=X+g(12)*(-1)
      Kc=EXP(-X/T)
      REVK(151)=kf/Kc*Y(11)*Y(15)

      kf=KFARRAY(152)
      FWDK(152)=kf*Y(12)*Y(16)
      X=0
      X=X+g(11)*(1)
      X=X+g(12)*(-1)
      Kc=EXP(-X/T)
      REVK(152)=kf/Kc*Y(11)*Y(16)

      kf=KFARRAY(153)
      FWDK(153)=kf*Y(12)*Y(16)
      X=0
      X=X+g(12)*(-1)
      X=X+g(15)*(1)
      X=X+g(16)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)
      REVK(153)=kf/Kc*Y(15)*Y(18)

      kf=KFARRAY(154)
      FWDK(154)=kf*Y(12)*Y(27)
      X=0
      X=X+g(12)*(-1)
      X=X+g(13)*(1)
      X=X+g(26)*(1)
      X=X+g(27)*(-1)
      Kc=EXP(-X/T)
      REVK(154)=kf/Kc*Y(13)*Y(26)

      kf=KFARRAY(155)
      FWDK(155)=kf*Y(4)*Y(13)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(13)*(-1)
      X=X+g(20)*(1)
      Kc=EXP(-X/T)
      REVK(155)=kf/Kc*Y(3)*Y(20)

      kf=KFARRAY(156)
      FWDK(156)=kf*Y(4)*Y(13)
      X=0
      X=X+g(4)*(-1)
      X=X+g(5)*(1)
      X=X+g(13)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)
      REVK(156)=kf/Kc*Y(5)*Y(18)

      kf=KFARRAY(157)
      FWDK(157)=kf*Y(8)*Y(13)
      X=0
      X=X+g(7)*(1)
      X=X+g(8)*(-1)
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      Kc=EXP(-X/T)
      REVK(157)=kf/Kc*Y(7)*Y(14)

C     lindemann approach
      k0=(T**(-7.030))*EXP(9.562976e+01 - (1.389884e+03)/T)
      kinf=KFARRAY(158)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.6190
      T3=73.20
      T1=1180.00
      T2=9999.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(158)=kf*Y(13)*Y(13)
      X=0
      X=X+g(13)*(-2)
      X=X+g(27)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(158)=kf/Kc*Y(27)

      kf=KFARRAY(159)
      FWDK(159)=kf*Y(13)*Y(13)
      X=0
      X=X+g(2)*(1)
      X=X+g(13)*(-2)
      X=X+g(26)*(1)
      Kc=EXP(-X/T)
      REVK(159)=kf/Kc*Y(2)*Y(26)

      kf=KFARRAY(160)
      FWDK(160)=kf*Y(13)*Y(17)
      X=0
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      X=X+g(15)*(1)
      X=X+g(17)*(-1)
      Kc=EXP(-X/T)
      REVK(160)=kf/Kc*Y(14)*Y(15)

      kf=KFARRAY(161)
      FWDK(161)=kf*Y(13)*Y(18)
      X=0
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      X=X+g(17)*(1)
      X=X+g(18)*(-1)
      Kc=EXP(-X/T)
      REVK(161)=kf/Kc*Y(14)*Y(17)

      kf=KFARRAY(162)
      FWDK(162)=kf*Y(13)*Y(21)
      X=0
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      X=X+g(19)*(1)
      X=X+g(21)*(-1)
      Kc=EXP(-X/T)
      REVK(162)=kf/Kc*Y(14)*Y(19)

      kf=KFARRAY(163)
      FWDK(163)=kf*Y(13)*Y(21)
      X=0
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      X=X+g(20)*(1)
      X=X+g(21)*(-1)
      Kc=EXP(-X/T)
      REVK(163)=kf/Kc*Y(14)*Y(20)

      kf=KFARRAY(164)
      FWDK(164)=kf*Y(13)*Y(25)
      X=0
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      X=X+g(24)*(1)
      X=X+g(25)*(-1)
      Kc=EXP(-X/T)
      REVK(164)=kf/Kc*Y(14)*Y(24)

      kf=KFARRAY(165)
      FWDK(165)=kf*Y(13)*Y(27)
      X=0
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      X=X+g(26)*(1)
      X=X+g(27)*(-1)
      Kc=EXP(-X/T)
      REVK(165)=kf/Kc*Y(14)*Y(26)

      kf=KFARRAY(166)
      FWDK(166)=kf*Y(6)*Y(17)
      X=0
      X=X+g(2)*(1)
      X=X+g(15)*(1)
      X=X+g(17)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(166)=kf/Kc*Y(2)*Y(6)*Y(15)

      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(-1)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      kf=KFARRAY(167)*Z
      FWDK(167)=kf*Y(17)
      X=0
      X=X+g(2)*(1)
      X=X+g(15)*(1)
      X=X+g(17)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(167)=kf/Kc*Y(2)*Y(15)

      kf=KFARRAY(168)
      FWDK(168)=kf*Y(4)*Y(17)
      X=0
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      X=X+g(15)*(1)
      X=X+g(17)*(-1)
      Kc=EXP(-X/T)
      REVK(168)=kf/Kc*Y(7)*Y(15)

      kf=KFARRAY(169)
      FWDK(169)=kf*Y(4)*Y(19)
      X=0
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      X=X+g(18)*(1)
      X=X+g(19)*(-1)
      Kc=EXP(-X/T)
      REVK(169)=kf/Kc*Y(7)*Y(18)

      kf=KFARRAY(170)
      FWDK(170)=kf*Y(4)*Y(20)
      X=0
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      X=X+g(18)*(1)
      X=X+g(20)*(-1)
      Kc=EXP(-X/T)
      REVK(170)=kf/Kc*Y(7)*Y(18)

      kf=KFARRAY(171)
      FWDK(171)=kf*Y(4)*Y(22)
      X=0
      X=X+g(4)*(-1)
      X=X+g(15)*(1)
      X=X+g(17)*(1)
      X=X+g(22)*(-1)
      Kc=EXP(-X/T)
      REVK(171)=kf/Kc*Y(15)*Y(17)

      kf=KFARRAY(172)
      FWDK(172)=kf*Y(1)*Y(22)
      X=0
      X=X+g(1)*(-1)
      X=X+g(2)*(1)
      X=X+g(22)*(-1)
      X=X+g(23)*(1)
      Kc=EXP(-X/T)
      REVK(172)=kf/Kc*Y(2)*Y(23)

      kf=KFARRAY(173)
      FWDK(173)=kf*Y(4)*Y(24)
      X=0
      X=X+g(4)*(-1)
      X=X+g(17)*(1)
      X=X+g(18)*(1)
      X=X+g(24)*(-1)
      Kc=EXP(-X/T)
      REVK(173)=kf/Kc*Y(17)*Y(18)

C     lindemann approach
      k0=(T**(-9.300))*EXP(1.178893e+02 - (4.921459e+04)/T)
      kinf=KFARRAY(174)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.7345
      T3=180.00
      T1=1035.00
      T2=5417.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(174)=kf*Y(25)
      X=0
      X=X+g(1)*(1)
      X=X+g(23)*(1)
      X=X+g(25)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(174)=kf/Kc*Y(1)*Y(23)

      kf=KFARRAY(175)
      FWDK(175)=kf*Y(4)*Y(26)
      X=0
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      X=X+g(25)*(1)
      X=X+g(26)*(-1)
      Kc=EXP(-X/T)
      REVK(175)=kf/Kc*Y(7)*Y(25)

      kf=KFARRAY(176)
      FWDK(176)=kf*Y(4)*Y(28)
      X=0
      X=X+g(4)*(-1)
      X=X+g(5)*(1)
      X=X+g(15)*(2)
      X=X+g(28)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(176)=kf/Kc*Y(5)*Y(15)*Y(15)

      kf=KFARRAY(177)
      FWDK(177)=kf*Y(28)*Y(28)
      X=0
      X=X+g(15)*(2)
      X=X+g(23)*(1)
      X=X+g(28)*(-2)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(177)=kf/Kc*Y(15)*Y(15)*Y(23)

      kf=KFARRAY(178)
      FWDK(178)=kf*Y(31)*Y(36)
      X=0
      X=X+g(3)*(1)
      X=X+g(31)*(-1)
      X=X+g(36)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(178)=kf/Kc*Y(3)*Y(48)

      kf=KFARRAY(179)
      FWDK(179)=kf*Y(4)*Y(31)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(31)*(-1)
      X=X+g(36)*(1)
      Kc=EXP(-X/T)
      REVK(179)=kf/Kc*Y(3)*Y(36)

      kf=KFARRAY(180)
      FWDK(180)=kf*Y(5)*Y(31)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(31)*(-1)
      X=X+g(36)*(1)
      Kc=EXP(-X/T)
      REVK(180)=kf/Kc*Y(2)*Y(36)

      kf=KFARRAY(181)
      FWDK(181)=kf*Y(3)*Y(38)
      X=0
      X=X+g(3)*(-1)
      X=X+g(4)*(1)
      X=X+g(38)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(181)=kf/Kc*Y(4)*Y(48)

      kf=KFARRAY(182)
      FWDK(182)=kf*Y(3)*Y(38)
      X=0
      X=X+g(3)*(-1)
      X=X+g(36)*(2)
      X=X+g(38)*(-1)
      Kc=EXP(-X/T)
      REVK(182)=kf/Kc*Y(36)*Y(36)

      kf=KFARRAY(183)
      FWDK(183)=kf*Y(2)*Y(38)
      X=0
      X=X+g(2)*(-1)
      X=X+g(5)*(1)
      X=X+g(38)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(183)=kf/Kc*Y(5)*Y(48)

      kf=KFARRAY(184)
      FWDK(184)=kf*Y(5)*Y(38)
      X=0
      X=X+g(5)*(-1)
      X=X+g(7)*(1)
      X=X+g(38)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(184)=kf/Kc*Y(7)*Y(48)

C     lindemann approach
      k0=(T**(.000))*EXP(3.408779e+01 - (2.850219e+04)/T)
      kinf=KFARRAY(185)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.375)
      Pr=k0*Z/kinf
      F=1
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(185)=kf*Y(38)
      X=0
      X=X+g(3)*(1)
      X=X+g(38)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(185)=kf/Kc*Y(3)*Y(48)

      kf=KFARRAY(186)
      FWDK(186)=kf*Y(7)*Y(36)
      X=0
      X=X+g(5)*(1)
      X=X+g(7)*(-1)
      X=X+g(36)*(-1)
      X=X+g(37)*(1)
      Kc=EXP(-X/T)
      REVK(186)=kf/Kc*Y(5)*Y(37)

      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      kf=KFARRAY(187)*Z
      FWDK(187)=kf*Y(3)*Y(36)
      X=0
      X=X+g(3)*(-1)
      X=X+g(36)*(-1)
      X=X+g(37)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(187)=kf/Kc*Y(37)

      kf=KFARRAY(188)
      FWDK(188)=kf*Y(3)*Y(37)
      X=0
      X=X+g(3)*(-1)
      X=X+g(4)*(1)
      X=X+g(36)*(1)
      X=X+g(37)*(-1)
      Kc=EXP(-X/T)
      REVK(188)=kf/Kc*Y(4)*Y(36)

      kf=KFARRAY(189)
      FWDK(189)=kf*Y(2)*Y(37)
      X=0
      X=X+g(2)*(-1)
      X=X+g(5)*(1)
      X=X+g(36)*(1)
      X=X+g(37)*(-1)
      Kc=EXP(-X/T)
      REVK(189)=kf/Kc*Y(5)*Y(36)

      kf=KFARRAY(190)
      FWDK(190)=kf*Y(3)*Y(32)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(32)*(-1)
      X=X+g(36)*(1)
      Kc=EXP(-X/T)
      REVK(190)=kf/Kc*Y(2)*Y(36)

      kf=KFARRAY(191)
      FWDK(191)=kf*Y(2)*Y(32)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(31)*(1)
      X=X+g(32)*(-1)
      Kc=EXP(-X/T)
      REVK(191)=kf/Kc*Y(1)*Y(31)

      kf=KFARRAY(192)
      FWDK(192)=kf*Y(5)*Y(32)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(32)*(-1)
      X=X+g(39)*(1)
      Kc=EXP(-X/T)
      REVK(192)=kf/Kc*Y(2)*Y(39)

      kf=KFARRAY(193)
      FWDK(193)=kf*Y(5)*Y(32)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(31)*(1)
      X=X+g(32)*(-1)
      Kc=EXP(-X/T)
      REVK(193)=kf/Kc*Y(6)*Y(31)

      kf=KFARRAY(194)
      FWDK(194)=kf*Y(4)*Y(32)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(32)*(-1)
      X=X+g(39)*(1)
      Kc=EXP(-X/T)
      REVK(194)=kf/Kc*Y(3)*Y(39)

      kf=KFARRAY(195)
      FWDK(195)=kf*Y(4)*Y(32)
      X=0
      X=X+g(4)*(-1)
      X=X+g(5)*(1)
      X=X+g(32)*(-1)
      X=X+g(36)*(1)
      Kc=EXP(-X/T)
      REVK(195)=kf/Kc*Y(5)*Y(36)

      kf=KFARRAY(196)
      FWDK(196)=kf*Y(31)*Y(32)
      X=0
      X=X+g(2)*(1)
      X=X+g(31)*(-1)
      X=X+g(32)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(196)=kf/Kc*Y(2)*Y(48)

      kf=KFARRAY(197)
      FWDK(197)=kf*Y(6)*Y(32)
      X=0
      X=X+g(1)*(1)
      X=X+g(6)*(-1)
      X=X+g(32)*(-1)
      X=X+g(39)*(1)
      Kc=EXP(-X/T)
      REVK(197)=kf/Kc*Y(1)*Y(39)

      kf=KFARRAY(198)
      FWDK(198)=kf*Y(32)*Y(36)
      X=0
      X=X+g(5)*(1)
      X=X+g(32)*(-1)
      X=X+g(36)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(198)=kf/Kc*Y(5)*Y(48)

      kf=KFARRAY(199)
      FWDK(199)=kf*Y(32)*Y(36)
      X=0
      X=X+g(2)*(1)
      X=X+g(32)*(-1)
      X=X+g(36)*(-1)
      X=X+g(38)*(1)
      Kc=EXP(-X/T)
      REVK(199)=kf/Kc*Y(2)*Y(38)

      kf=KFARRAY(200)
      FWDK(200)=kf*Y(3)*Y(33)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(32)*(1)
      X=X+g(33)*(-1)
      Kc=EXP(-X/T)
      REVK(200)=kf/Kc*Y(5)*Y(32)

      kf=KFARRAY(201)
      FWDK(201)=kf*Y(3)*Y(33)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(33)*(-1)
      X=X+g(39)*(1)
      Kc=EXP(-X/T)
      REVK(201)=kf/Kc*Y(2)*Y(39)

      kf=KFARRAY(202)
      FWDK(202)=kf*Y(2)*Y(33)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(32)*(1)
      X=X+g(33)*(-1)
      Kc=EXP(-X/T)
      REVK(202)=kf/Kc*Y(1)*Y(32)

      kf=KFARRAY(203)
      FWDK(203)=kf*Y(5)*Y(33)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(32)*(1)
      X=X+g(33)*(-1)
      Kc=EXP(-X/T)
      REVK(203)=kf/Kc*Y(6)*Y(32)

      kf=KFARRAY(204)
      FWDK(204)=kf*Y(35)
      X=0
      X=X+g(2)*(1)
      X=X+g(35)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(204)=kf/Kc*Y(2)*Y(48)

      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      kf=KFARRAY(205)*Z
      FWDK(205)=kf*Y(35)
      X=0
      X=X+g(2)*(1)
      X=X+g(35)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(205)=kf/Kc*Y(2)*Y(48)

      kf=KFARRAY(206)
      FWDK(206)=kf*Y(4)*Y(35)
      X=0
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      X=X+g(35)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(206)=kf/Kc*Y(7)*Y(48)

      kf=KFARRAY(207)
      FWDK(207)=kf*Y(3)*Y(35)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(35)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(207)=kf/Kc*Y(5)*Y(48)

      kf=KFARRAY(208)
      FWDK(208)=kf*Y(3)*Y(35)
      X=0
      X=X+g(3)*(-1)
      X=X+g(32)*(1)
      X=X+g(35)*(-1)
      X=X+g(36)*(1)
      Kc=EXP(-X/T)
      REVK(208)=kf/Kc*Y(32)*Y(36)

      kf=KFARRAY(209)
      FWDK(209)=kf*Y(2)*Y(35)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(35)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(209)=kf/Kc*Y(1)*Y(48)

      kf=KFARRAY(210)
      FWDK(210)=kf*Y(5)*Y(35)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(35)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(210)=kf/Kc*Y(6)*Y(48)

      kf=KFARRAY(211)
      FWDK(211)=kf*Y(13)*Y(35)
      X=0
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      X=X+g(35)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(211)=kf/Kc*Y(14)*Y(48)

      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      kf=KFARRAY(212)*Z
      FWDK(212)=kf*Y(2)*Y(36)
      X=0
      X=X+g(2)*(-1)
      X=X+g(36)*(-1)
      X=X+g(39)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(212)=kf/Kc*Y(39)

      kf=KFARRAY(213)
      FWDK(213)=kf*Y(3)*Y(39)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(36)*(1)
      X=X+g(39)*(-1)
      Kc=EXP(-X/T)
      REVK(213)=kf/Kc*Y(5)*Y(36)

      kf=KFARRAY(214)
      FWDK(214)=kf*Y(2)*Y(39)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(36)*(1)
      X=X+g(39)*(-1)
      Kc=EXP(-X/T)
      REVK(214)=kf/Kc*Y(1)*Y(36)

      kf=KFARRAY(215)
      FWDK(215)=kf*Y(5)*Y(39)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(36)*(1)
      X=X+g(39)*(-1)
      Kc=EXP(-X/T)
      REVK(215)=kf/Kc*Y(6)*Y(36)

      kf=KFARRAY(216)
      FWDK(216)=kf*Y(4)*Y(39)
      X=0
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      X=X+g(36)*(1)
      X=X+g(39)*(-1)
      Kc=EXP(-X/T)
      REVK(216)=kf/Kc*Y(7)*Y(36)

      kf=KFARRAY(217)
      FWDK(217)=kf*Y(3)*Y(40)
      X=0
      X=X+g(3)*(-1)
      X=X+g(15)*(1)
      X=X+g(31)*(1)
      X=X+g(40)*(-1)
      Kc=EXP(-X/T)
      REVK(217)=kf/Kc*Y(15)*Y(31)

      kf=KFARRAY(218)
      FWDK(218)=kf*Y(5)*Y(40)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(40)*(-1)
      X=X+g(47)*(1)
      Kc=EXP(-X/T)
      REVK(218)=kf/Kc*Y(2)*Y(47)

      kf=KFARRAY(219)
      FWDK(219)=kf*Y(6)*Y(40)
      X=0
      X=X+g(5)*(1)
      X=X+g(6)*(-1)
      X=X+g(40)*(-1)
      X=X+g(41)*(1)
      Kc=EXP(-X/T)
      REVK(219)=kf/Kc*Y(5)*Y(41)

      kf=KFARRAY(220)
      FWDK(220)=kf*Y(4)*Y(40)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(40)*(-1)
      X=X+g(47)*(1)
      Kc=EXP(-X/T)
      REVK(220)=kf/Kc*Y(3)*Y(47)

      kf=KFARRAY(221)
      FWDK(221)=kf*Y(1)*Y(40)
      X=0
      X=X+g(1)*(-1)
      X=X+g(2)*(1)
      X=X+g(40)*(-1)
      X=X+g(41)*(1)
      Kc=EXP(-X/T)
      REVK(221)=kf/Kc*Y(2)*Y(41)

      kf=KFARRAY(222)
      FWDK(222)=kf*Y(3)*Y(47)
      X=0
      X=X+g(3)*(-1)
      X=X+g(15)*(1)
      X=X+g(36)*(1)
      X=X+g(47)*(-1)
      Kc=EXP(-X/T)
      REVK(222)=kf/Kc*Y(15)*Y(36)

      kf=KFARRAY(223)
      FWDK(223)=kf*Y(2)*Y(47)
      X=0
      X=X+g(2)*(-1)
      X=X+g(15)*(1)
      X=X+g(32)*(1)
      X=X+g(47)*(-1)
      Kc=EXP(-X/T)
      REVK(223)=kf/Kc*Y(15)*Y(32)

      kf=KFARRAY(224)
      FWDK(224)=kf*Y(5)*Y(47)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(15)*(1)
      X=X+g(36)*(1)
      X=X+g(47)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(224)=kf/Kc*Y(2)*Y(15)*Y(36)

      kf=KFARRAY(225)
      FWDK(225)=kf*Y(31)*Y(47)
      X=0
      X=X+g(15)*(1)
      X=X+g(31)*(-1)
      X=X+g(47)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(225)=kf/Kc*Y(15)*Y(48)

      kf=KFARRAY(226)
      FWDK(226)=kf*Y(4)*Y(47)
      X=0
      X=X+g(4)*(-1)
      X=X+g(16)*(1)
      X=X+g(36)*(1)
      X=X+g(47)*(-1)
      Kc=EXP(-X/T)
      REVK(226)=kf/Kc*Y(16)*Y(36)

      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      kf=KFARRAY(227)*Z
      FWDK(227)=kf*Y(47)
      X=0
      X=X+g(15)*(1)
      X=X+g(31)*(1)
      X=X+g(47)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(227)=kf/Kc*Y(15)*Y(31)

      kf=KFARRAY(228)
      FWDK(228)=kf*Y(36)*Y(47)
      X=0
      X=X+g(15)*(1)
      X=X+g(36)*(-1)
      X=X+g(38)*(1)
      X=X+g(47)*(-1)
      Kc=EXP(-X/T)
      REVK(228)=kf/Kc*Y(15)*Y(38)

      kf=KFARRAY(229)
      FWDK(229)=kf*Y(36)*Y(47)
      X=0
      X=X+g(16)*(1)
      X=X+g(36)*(-1)
      X=X+g(47)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(229)=kf/Kc*Y(16)*Y(48)

      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      kf=KFARRAY(230)*Z
      FWDK(230)=kf*Y(41)
      X=0
      X=X+g(2)*(1)
      X=X+g(40)*(1)
      X=X+g(41)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(230)=kf/Kc*Y(2)*Y(40)

      kf=KFARRAY(231)
      FWDK(231)=kf*Y(3)*Y(41)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(41)*(-1)
      X=X+g(47)*(1)
      Kc=EXP(-X/T)
      REVK(231)=kf/Kc*Y(2)*Y(47)

      kf=KFARRAY(232)
      FWDK(232)=kf*Y(3)*Y(41)
      X=0
      X=X+g(3)*(-1)
      X=X+g(15)*(1)
      X=X+g(32)*(1)
      X=X+g(41)*(-1)
      Kc=EXP(-X/T)
      REVK(232)=kf/Kc*Y(15)*Y(32)

      kf=KFARRAY(233)
      FWDK(233)=kf*Y(3)*Y(41)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(40)*(1)
      X=X+g(41)*(-1)
      Kc=EXP(-X/T)
      REVK(233)=kf/Kc*Y(5)*Y(40)

      kf=KFARRAY(234)
      FWDK(234)=kf*Y(5)*Y(41)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(41)*(-1)
      X=X+g(45)*(1)
      Kc=EXP(-X/T)
      REVK(234)=kf/Kc*Y(2)*Y(45)

      kf=KFARRAY(235)
      FWDK(235)=kf*Y(5)*Y(41)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(41)*(-1)
      X=X+g(46)*(1)
      Kc=EXP(-X/T)
      REVK(235)=kf/Kc*Y(2)*Y(46)

      kf=KFARRAY(236)
      FWDK(236)=kf*Y(5)*Y(41)
      X=0
      X=X+g(5)*(-1)
      X=X+g(15)*(1)
      X=X+g(33)*(1)
      X=X+g(41)*(-1)
      Kc=EXP(-X/T)
      REVK(236)=kf/Kc*Y(15)*Y(33)

C     lindemann approach
      k0=(T**(-3.400))*EXP(6.020368e+01 - (9.561117e+02)/T)
      kinf=KFARRAY(237)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
      F=1
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(237)=kf*Y(2)*Y(41)
      X=0
      X=X+g(2)*(-1)
      X=X+g(41)*(-1)
      X=X+g(42)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(237)=kf/Kc*Y(42)

      kf=KFARRAY(238)
      FWDK(238)=kf*Y(31)*Y(42)
      X=0
      X=X+g(11)*(1)
      X=X+g(31)*(-1)
      X=X+g(42)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(238)=kf/Kc*Y(11)*Y(48)

      kf=KFARRAY(239)
      FWDK(239)=kf*Y(9)*Y(48)
      X=0
      X=X+g(9)*(-1)
      X=X+g(31)*(1)
      X=X+g(40)*(1)
      X=X+g(48)*(-1)
      Kc=EXP(-X/T)
      REVK(239)=kf/Kc*Y(31)*Y(40)

      kf=KFARRAY(240)
      FWDK(240)=kf*Y(10)*Y(48)
      X=0
      X=X+g(10)*(-1)
      X=X+g(31)*(1)
      X=X+g(41)*(1)
      X=X+g(48)*(-1)
      Kc=EXP(-X/T)
      REVK(240)=kf/Kc*Y(31)*Y(41)

C     lindemann approach
      k0=(T**(-3.160))*EXP(5.782699e+01 - (3.723803e+02)/T)
      kinf=KFARRAY(241)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(0)
      Pr=k0*Z/kinf
C     Troe form
      a=.6670
      T3=235.00
      T1=2117.00
      T2=4536.00
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(241)=kf*Y(10)*Y(48)
      X=0
      X=X+g(10)*(-1)
      X=X+g(43)*(1)
      X=X+g(48)*(-1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(241)=kf/Kc*Y(43)

      kf=KFARRAY(242)
      FWDK(242)=kf*Y(11)*Y(48)
      X=0
      X=X+g(11)*(-1)
      X=X+g(32)*(1)
      X=X+g(41)*(1)
      X=X+g(48)*(-1)
      Kc=EXP(-X/T)
      REVK(242)=kf/Kc*Y(32)*Y(41)

      kf=KFARRAY(243)
      FWDK(243)=kf*Y(12)*Y(48)
      X=0
      X=X+g(12)*(-1)
      X=X+g(32)*(1)
      X=X+g(41)*(1)
      X=X+g(48)*(-1)
      Kc=EXP(-X/T)
      REVK(243)=kf/Kc*Y(32)*Y(41)

      kf=KFARRAY(244)
      FWDK(244)=kf*Y(9)*Y(36)
      X=0
      X=X+g(3)*(1)
      X=X+g(9)*(-1)
      X=X+g(36)*(-1)
      X=X+g(40)*(1)
      Kc=EXP(-X/T)
      REVK(244)=kf/Kc*Y(3)*Y(40)

      kf=KFARRAY(245)
      FWDK(245)=kf*Y(9)*Y(36)
      X=0
      X=X+g(9)*(-1)
      X=X+g(15)*(1)
      X=X+g(31)*(1)
      X=X+g(36)*(-1)
      Kc=EXP(-X/T)
      REVK(245)=kf/Kc*Y(15)*Y(31)

      kf=KFARRAY(246)
      FWDK(246)=kf*Y(10)*Y(36)
      X=0
      X=X+g(3)*(1)
      X=X+g(10)*(-1)
      X=X+g(36)*(-1)
      X=X+g(41)*(1)
      Kc=EXP(-X/T)
      REVK(246)=kf/Kc*Y(3)*Y(41)

      kf=KFARRAY(247)
      FWDK(247)=kf*Y(10)*Y(36)
      X=0
      X=X+g(2)*(1)
      X=X+g(10)*(-1)
      X=X+g(36)*(-1)
      X=X+g(47)*(1)
      Kc=EXP(-X/T)
      REVK(247)=kf/Kc*Y(2)*Y(47)

      kf=KFARRAY(248)
      FWDK(248)=kf*Y(10)*Y(36)
      X=0
      X=X+g(10)*(-1)
      X=X+g(17)*(1)
      X=X+g(31)*(1)
      X=X+g(36)*(-1)
      Kc=EXP(-X/T)
      REVK(248)=kf/Kc*Y(17)*Y(31)

      kf=KFARRAY(249)
      FWDK(249)=kf*Y(11)*Y(36)
      X=0
      X=X+g(2)*(1)
      X=X+g(11)*(-1)
      X=X+g(36)*(-1)
      X=X+g(46)*(1)
      Kc=EXP(-X/T)
      REVK(249)=kf/Kc*Y(2)*Y(46)

      kf=KFARRAY(250)
      FWDK(250)=kf*Y(11)*Y(36)
      X=0
      X=X+g(5)*(1)
      X=X+g(11)*(-1)
      X=X+g(36)*(-1)
      X=X+g(41)*(1)
      Kc=EXP(-X/T)
      REVK(250)=kf/Kc*Y(5)*Y(41)

      kf=KFARRAY(251)
      FWDK(251)=kf*Y(11)*Y(36)
      X=0
      X=X+g(2)*(1)
      X=X+g(11)*(-1)
      X=X+g(36)*(-1)
      X=X+g(44)*(1)
      Kc=EXP(-X/T)
      REVK(251)=kf/Kc*Y(2)*Y(44)

      kf=KFARRAY(252)
      FWDK(252)=kf*Y(12)*Y(36)
      X=0
      X=X+g(2)*(1)
      X=X+g(12)*(-1)
      X=X+g(36)*(-1)
      X=X+g(46)*(1)
      Kc=EXP(-X/T)
      REVK(252)=kf/Kc*Y(2)*Y(46)

      kf=KFARRAY(253)
      FWDK(253)=kf*Y(12)*Y(36)
      X=0
      X=X+g(5)*(1)
      X=X+g(12)*(-1)
      X=X+g(36)*(-1)
      X=X+g(41)*(1)
      Kc=EXP(-X/T)
      REVK(253)=kf/Kc*Y(5)*Y(41)

      kf=KFARRAY(254)
      FWDK(254)=kf*Y(12)*Y(36)
      X=0
      X=X+g(2)*(1)
      X=X+g(12)*(-1)
      X=X+g(36)*(-1)
      X=X+g(44)*(1)
      Kc=EXP(-X/T)
      REVK(254)=kf/Kc*Y(2)*Y(44)

      kf=KFARRAY(255)
      FWDK(255)=kf*Y(13)*Y(36)
      X=0
      X=X+g(6)*(1)
      X=X+g(13)*(-1)
      X=X+g(36)*(-1)
      X=X+g(41)*(1)
      Kc=EXP(-X/T)
      REVK(255)=kf/Kc*Y(6)*Y(41)

      kf=KFARRAY(256)
      FWDK(256)=kf*Y(13)*Y(36)
      X=0
      X=X+g(5)*(1)
      X=X+g(13)*(-1)
      X=X+g(36)*(-1)
      X=X+g(42)*(1)
      Kc=EXP(-X/T)
      REVK(256)=kf/Kc*Y(5)*Y(42)

      kf=KFARRAY(257)
      FWDK(257)=kf*Y(3)*Y(43)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(15)*(1)
      X=X+g(43)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(257)=kf/Kc*Y(2)*Y(15)*Y(48)

      kf=KFARRAY(258)
      FWDK(258)=kf*Y(3)*Y(43)
      X=0
      X=X+g(3)*(-1)
      X=X+g(36)*(1)
      X=X+g(41)*(1)
      X=X+g(43)*(-1)
      Kc=EXP(-X/T)
      REVK(258)=kf/Kc*Y(36)*Y(41)

      kf=KFARRAY(259)
      FWDK(259)=kf*Y(4)*Y(43)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(17)*(1)
      X=X+g(43)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(259)=kf/Kc*Y(3)*Y(17)*Y(48)

      kf=KFARRAY(260)
      FWDK(260)=kf*Y(5)*Y(43)
      X=0
      X=X+g(2)*(1)
      X=X+g(5)*(-1)
      X=X+g(17)*(1)
      X=X+g(43)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(260)=kf/Kc*Y(2)*Y(17)*Y(48)

      kf=KFARRAY(261)
      FWDK(261)=kf*Y(2)*Y(43)
      X=0
      X=X+g(2)*(-1)
      X=X+g(11)*(1)
      X=X+g(43)*(-1)
      X=X+g(48)*(1)
      Kc=EXP(-X/T)
      REVK(261)=kf/Kc*Y(11)*Y(48)

      kf=KFARRAY(262)
      FWDK(262)=kf*Y(3)*Y(46)
      X=0
      X=X+g(3)*(-1)
      X=X+g(16)*(1)
      X=X+g(32)*(1)
      X=X+g(46)*(-1)
      Kc=EXP(-X/T)
      REVK(262)=kf/Kc*Y(16)*Y(32)

      kf=KFARRAY(263)
      FWDK(263)=kf*Y(3)*Y(46)
      X=0
      X=X+g(3)*(-1)
      X=X+g(15)*(1)
      X=X+g(39)*(1)
      X=X+g(46)*(-1)
      Kc=EXP(-X/T)
      REVK(263)=kf/Kc*Y(15)*Y(39)

      kf=KFARRAY(264)
      FWDK(264)=kf*Y(3)*Y(46)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(46)*(-1)
      X=X+g(47)*(1)
      Kc=EXP(-X/T)
      REVK(264)=kf/Kc*Y(5)*Y(47)

      kf=KFARRAY(265)
      FWDK(265)=kf*Y(2)*Y(46)
      X=0
      X=X+g(2)*(-1)
      X=X+g(15)*(1)
      X=X+g(33)*(1)
      X=X+g(46)*(-1)
      Kc=EXP(-X/T)
      REVK(265)=kf/Kc*Y(15)*Y(33)

      kf=KFARRAY(266)
      FWDK(266)=kf*Y(2)*Y(46)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(46)*(-1)
      X=X+g(47)*(1)
      Kc=EXP(-X/T)
      REVK(266)=kf/Kc*Y(1)*Y(47)

      kf=KFARRAY(267)
      FWDK(267)=kf*Y(5)*Y(46)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(46)*(-1)
      X=X+g(47)*(1)
      Kc=EXP(-X/T)
      REVK(267)=kf/Kc*Y(6)*Y(47)

      kf=KFARRAY(268)
      FWDK(268)=kf*Y(5)*Y(46)
      X=0
      X=X+g(5)*(-1)
      X=X+g(16)*(1)
      X=X+g(33)*(1)
      X=X+g(46)*(-1)
      Kc=EXP(-X/T)
      REVK(268)=kf/Kc*Y(16)*Y(33)

      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      kf=KFARRAY(269)*Z
      FWDK(269)=kf*Y(46)
      X=0
      X=X+g(15)*(1)
      X=X+g(32)*(1)
      X=X+g(46)*(-1)
      Kc=EXP(-X/T)/(T/PATM)
      REVK(269)=kf/Kc*Y(15)*Y(32)

      kf=KFARRAY(270)
      FWDK(270)=kf*Y(2)*Y(44)
      X=0
      X=X+g(44)*(-1)
      X=X+g(46)*(1)
      Kc=EXP(-X/T)
      REVK(270)=kf/Kc*Y(2)*Y(46)

      kf=KFARRAY(271)
      FWDK(271)=kf*Y(2)*Y(44)
      X=0
      X=X+g(2)*(-1)
      X=X+g(5)*(1)
      X=X+g(41)*(1)
      X=X+g(44)*(-1)
      Kc=EXP(-X/T)
      REVK(271)=kf/Kc*Y(5)*Y(41)

      kf=KFARRAY(272)
      FWDK(272)=kf*Y(2)*Y(44)
      X=0
      X=X+g(2)*(-1)
      X=X+g(15)*(1)
      X=X+g(33)*(1)
      X=X+g(44)*(-1)
      Kc=EXP(-X/T)
      REVK(272)=kf/Kc*Y(15)*Y(33)

      kf=KFARRAY(273)
      FWDK(273)=kf*Y(2)*Y(45)
      X=0
      X=X+g(45)*(-1)
      X=X+g(46)*(1)
      Kc=EXP(-X/T)
      REVK(273)=kf/Kc*Y(2)*Y(46)

      kf=KFARRAY(274)
      FWDK(274)=kf*Y(28)*Y(36)
      X=0
      X=X+g(15)*(1)
      X=X+g(28)*(-1)
      X=X+g(36)*(-1)
      X=X+g(44)*(1)
      Kc=EXP(-X/T)
      REVK(274)=kf/Kc*Y(15)*Y(44)

      kf=KFARRAY(275)
      FWDK(275)=kf*Y(13)*Y(31)
      X=0
      X=X+g(2)*(1)
      X=X+g(13)*(-1)
      X=X+g(31)*(-1)
      X=X+g(42)*(1)
      Kc=EXP(-X/T)
      REVK(275)=kf/Kc*Y(2)*Y(42)

      kf=KFARRAY(276)
      FWDK(276)=kf*Y(13)*Y(31)
      X=0
      X=X+g(1)*(1)
      X=X+g(13)*(-1)
      X=X+g(31)*(-1)
      X=X+g(41)*(1)
      Kc=EXP(-X/T)
      REVK(276)=kf/Kc*Y(1)*Y(41)

      kf=KFARRAY(277)
      FWDK(277)=kf*Y(2)*Y(34)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(33)*(1)
      X=X+g(34)*(-1)
      Kc=EXP(-X/T)
      REVK(277)=kf/Kc*Y(1)*Y(33)

      kf=KFARRAY(278)
      FWDK(278)=kf*Y(5)*Y(34)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(33)*(1)
      X=X+g(34)*(-1)
      Kc=EXP(-X/T)
      REVK(278)=kf/Kc*Y(6)*Y(33)

      kf=KFARRAY(279)
      FWDK(279)=kf*Y(3)*Y(34)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(33)*(1)
      X=X+g(34)*(-1)
      Kc=EXP(-X/T)
      REVK(279)=kf/Kc*Y(5)*Y(33)

      kf=KFARRAY(280)
      FWDK(280)=kf*Y(16)*Y(32)
      X=0
      X=X+g(15)*(1)
      X=X+g(16)*(-1)
      X=X+g(32)*(-1)
      X=X+g(39)*(1)
      Kc=EXP(-X/T)
      REVK(280)=kf/Kc*Y(15)*Y(39)

      kf=KFARRAY(281)
      FWDK(281)=kf*Y(37)*Y(40)
      X=0
      X=X+g(36)*(1)
      X=X+g(37)*(-1)
      X=X+g(40)*(-1)
      X=X+g(47)*(1)
      Kc=EXP(-X/T)
      REVK(281)=kf/Kc*Y(36)*Y(47)

      kf=KFARRAY(282)
      FWDK(282)=kf*Y(37)*Y(47)
      X=0
      X=X+g(16)*(1)
      X=X+g(37)*(-1)
      X=X+g(38)*(1)
      X=X+g(47)*(-1)
      Kc=EXP(-X/T)
      REVK(282)=kf/Kc*Y(16)*Y(38)

      kf=KFARRAY(283)
      FWDK(283)=kf*Y(16)*Y(31)
      X=0
      X=X+g(15)*(1)
      X=X+g(16)*(-1)
      X=X+g(31)*(-1)
      X=X+g(36)*(1)
      Kc=EXP(-X/T)
      REVK(283)=kf/Kc*Y(15)*Y(36)

      kf=KFARRAY(284)
      FWDK(284)=kf*Y(3)*Y(13)

      kf=KFARRAY(285)
      FWDK(285)=kf*Y(3)*Y(25)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(25)*(-1)
      X=X+g(52)*(1)
      Kc=EXP(-X/T)
      REVK(285)=kf/Kc*Y(2)*Y(52)

      kf=KFARRAY(286)
      FWDK(286)=kf*Y(3)*Y(26)
      X=0
      X=X+g(2)*(1)
      X=X+g(3)*(-1)
      X=X+g(26)*(-1)
      X=X+g(53)*(1)
      Kc=EXP(-X/T)
      REVK(286)=kf/Kc*Y(2)*Y(53)

      kf=KFARRAY(287)
      FWDK(287)=kf*Y(5)*Y(7)
      X=0
      X=X+g(4)*(1)
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(7)*(-1)
      Kc=EXP(-X/T)
      REVK(287)=kf/Kc*Y(4)*Y(6)

      kf=KFARRAY(288)
      FWDK(288)=kf*Y(5)*Y(13)

C     lindemann approach
      k0=(T**(-2.80))*EXP(5.913740e+01 - (2.968978e+02)/T)
      kinf=KFARRAY(289)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.578
      T3=122.0
      T1=2535.0
      T2=9365.0
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(289)=kf*Y(1)*Y(10)
      X=0
      X=X+g(1)*(-1)
      X=X+g(10)*(-1)
      X=X+g(13)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(289)=kf/Kc*Y(13)

      kf=KFARRAY(290)
      FWDK(290)=kf*Y(4)*Y(11)

      kf=KFARRAY(291)
      FWDK(291)=kf*Y(4)*Y(11)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(11)*(-1)
      X=X+g(18)*(1)
      Kc=EXP(-X/T)
      REVK(291)=kf/Kc*Y(3)*Y(18)

      kf=KFARRAY(292)
      FWDK(292)=kf*Y(11)*Y(11)

      kf=KFARRAY(293)
      FWDK(293)=kf*Y(6)*Y(12)

      kf=KFARRAY(294)
      FWDK(294)=kf*Y(4)*Y(24)
      X=0
      X=X+g(3)*(1)
      X=X+g(4)*(-1)
      X=X+g(24)*(-1)
      X=X+g(52)*(1)
      Kc=EXP(-X/T)
      REVK(294)=kf/Kc*Y(3)*Y(52)

      kf=KFARRAY(295)
      FWDK(295)=kf*Y(4)*Y(24)
      X=0
      X=X+g(4)*(-1)
      X=X+g(7)*(1)
      X=X+g(23)*(1)
      X=X+g(24)*(-1)
      Kc=EXP(-X/T)
      REVK(295)=kf/Kc*Y(7)*Y(23)

      kf=KFARRAY(296)
      FWDK(296)=kf*Y(3)*Y(53)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(52)*(1)
      X=X+g(53)*(-1)
      Kc=EXP(-X/T)
      REVK(296)=kf/Kc*Y(5)*Y(52)

      kf=KFARRAY(297)
      FWDK(297)=kf*Y(3)*Y(53)

      kf=KFARRAY(298)
      FWDK(298)=kf*Y(4)*Y(53)

      kf=KFARRAY(299)
      FWDK(299)=kf*Y(2)*Y(53)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(52)*(1)
      X=X+g(53)*(-1)
      Kc=EXP(-X/T)
      REVK(299)=kf/Kc*Y(1)*Y(52)

      kf=KFARRAY(300)
      FWDK(300)=kf*Y(2)*Y(53)

      kf=KFARRAY(301)
      FWDK(301)=kf*Y(5)*Y(53)

      kf=KFARRAY(302)
      FWDK(302)=kf*Y(7)*Y(53)

      kf=KFARRAY(303)
      FWDK(303)=kf*Y(13)*Y(53)

C     lindemann approach
      k0=(T**(-7.63))*EXP(9.672050e+01 - (1.939397e+03)/T)
      kinf=KFARRAY(304)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=0.465
      T3=201.0
      T1=1773.0
      T2=5333.0
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(304)=kf*Y(2)*Y(29)
      X=0
      X=X+g(2)*(-1)
      X=X+g(29)*(-1)
      X=X+g(52)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(304)=kf/Kc*Y(52)

      kf=KFARRAY(305)
      FWDK(305)=kf*Y(3)*Y(52)

      kf=KFARRAY(306)
      FWDK(306)=kf*Y(4)*Y(52)

      kf=KFARRAY(307)
      FWDK(307)=kf*Y(4)*Y(52)

      kf=KFARRAY(308)
      FWDK(308)=kf*Y(2)*Y(52)
      X=0
      X=X+g(2)*(-1)
      X=X+g(13)*(1)
      X=X+g(17)*(1)
      X=X+g(52)*(-1)
      Kc=EXP(-X/T)
      REVK(308)=kf/Kc*Y(13)*Y(17)

      kf=KFARRAY(309)
      FWDK(309)=kf*Y(2)*Y(52)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(29)*(1)
      X=X+g(52)*(-1)
      Kc=EXP(-X/T)
      REVK(309)=kf/Kc*Y(1)*Y(29)

      kf=KFARRAY(310)
      FWDK(310)=kf*Y(5)*Y(52)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(29)*(1)
      X=X+g(52)*(-1)
      Kc=EXP(-X/T)
      REVK(310)=kf/Kc*Y(6)*Y(29)

      kf=KFARRAY(311)
      FWDK(311)=kf*Y(5)*Y(52)
      X=0
      X=X+g(5)*(-1)
      X=X+g(17)*(1)
      X=X+g(19)*(1)
      X=X+g(52)*(-1)
      Kc=EXP(-X/T)
      REVK(311)=kf/Kc*Y(17)*Y(19)

C     lindemann approach
      k0=(T**(-16.82))*EXP(1.713882e+02 - (6.574526e+03)/T)
      kinf=KFARRAY(312)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.1527
      T3=291.0
      T1=2742.0
      T2=7748.0
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(312)=kf*Y(13)*Y(26)
      X=0
      X=X+g(13)*(-1)
      X=X+g(26)*(-1)
      X=X+g(51)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(312)=kf/Kc*Y(51)

      kf=KFARRAY(313)
      FWDK(313)=kf*Y(3)*Y(51)
      X=0
      X=X+g(3)*(-1)
      X=X+g(5)*(1)
      X=X+g(50)*(1)
      X=X+g(51)*(-1)
      Kc=EXP(-X/T)
      REVK(313)=kf/Kc*Y(5)*Y(50)

      kf=KFARRAY(314)
      FWDK(314)=kf*Y(2)*Y(51)
      X=0
      X=X+g(1)*(1)
      X=X+g(2)*(-1)
      X=X+g(50)*(1)
      X=X+g(51)*(-1)
      Kc=EXP(-X/T)
      REVK(314)=kf/Kc*Y(1)*Y(50)

      kf=KFARRAY(315)
      FWDK(315)=kf*Y(5)*Y(51)
      X=0
      X=X+g(5)*(-1)
      X=X+g(6)*(1)
      X=X+g(50)*(1)
      X=X+g(51)*(-1)
      Kc=EXP(-X/T)
      REVK(315)=kf/Kc*Y(6)*Y(50)

      kf=KFARRAY(316)
      FWDK(316)=kf*Y(8)*Y(50)
      X=0
      X=X+g(7)*(1)
      X=X+g(8)*(-1)
      X=X+g(50)*(-1)
      X=X+g(51)*(1)
      Kc=EXP(-X/T)
      REVK(316)=kf/Kc*Y(7)*Y(51)

      kf=KFARRAY(317)
      FWDK(317)=kf*Y(13)*Y(51)
      X=0
      X=X+g(13)*(-1)
      X=X+g(14)*(1)
      X=X+g(50)*(1)
      X=X+g(51)*(-1)
      Kc=EXP(-X/T)
      REVK(317)=kf/Kc*Y(14)*Y(50)

C     lindemann approach
      k0=(T**(-14.6))*EXP(1.461615e+02 - (9.143447e+03)/T)
      kinf=KFARRAY(318)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.1894
      T3=277.0
      T1=8748.0
      T2=7891.0
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(318)=kf*Y(13)*Y(25)
      X=0
      X=X+g(13)*(-1)
      X=X+g(25)*(-1)
      X=X+g(50)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(318)=kf/Kc*Y(50)

      kf=KFARRAY(319)
      FWDK(319)=kf*Y(3)*Y(50)
      X=0
      X=X+g(3)*(-1)
      X=X+g(18)*(1)
      X=X+g(26)*(1)
      X=X+g(50)*(-1)
      Kc=EXP(-X/T)
      REVK(319)=kf/Kc*Y(18)*Y(26)

C     lindemann approach
      k0=(T**(-13.545))*EXP(1.419438e+02 - (5.715032e+03)/T)
      kinf=KFARRAY(320)
      Z=sumY
      Z=Z+Y(1)*(1)
      Z=Z+Y(6)*(5)
      Z=Z+Y(14)*(1)
      Z=Z+Y(15)*(0.5)
      Z=Z+Y(16)*(1)
      Z=Z+Y(27)*(2)
      Z=Z+Y(49)*(-0.3)
      Pr=k0*Z/kinf
C     Troe form
      a=.315
      T3=369.0
      T1=3285.0
      T2=6667.0
      Fcent=(1-a)*EXP(-T/T3)+a*EXP(-T/T1)+EXP(-T2/T);
      c = -0.4 - 0.67*LOG10(Fcent)
      n = 0.75 - 1.27*LOG10(Fcent)
      F=10**(LOG10(Fcent)*(1+(  (LOG10(Pr)+c)
     *  /(n-d*(LOG10(Pr)+c))  )**2)**(-1))
      kf=kinf*(Pr/(1+Pr))*F
      FWDK(320)=kf*Y(2)*Y(50)
      X=0
      X=X+g(2)*(-1)
      X=X+g(50)*(-1)
      X=X+g(51)*(1)
      Kc=EXP(-X/T)*(T/PATM)
      REVK(320)=kf/Kc*Y(51)

      kf=KFARRAY(321)
      FWDK(321)=kf*Y(2)*Y(50)
      X=0
      X=X+g(2)*(-1)
      X=X+g(13)*(1)
      X=X+g(26)*(1)
      X=X+g(50)*(-1)
      Kc=EXP(-X/T)
      REVK(321)=kf/Kc*Y(13)*Y(26)

      kf=KFARRAY(322)
      FWDK(322)=kf*Y(5)*Y(50)
      X=0
      X=X+g(5)*(-1)
      X=X+g(19)*(1)
      X=X+g(26)*(1)
      X=X+g(50)*(-1)
      Kc=EXP(-X/T)
      REVK(322)=kf/Kc*Y(19)*Y(26)

      kf=KFARRAY(323)
      FWDK(323)=kf*Y(7)*Y(50)
      X=0
      X=X+g(4)*(1)
      X=X+g(7)*(-1)
      X=X+g(50)*(-1)
      X=X+g(51)*(1)
      Kc=EXP(-X/T)
      REVK(323)=kf/Kc*Y(4)*Y(51)

      kf=KFARRAY(324)
      FWDK(324)=kf*Y(7)*Y(50)

      kf=KFARRAY(325)
      FWDK(325)=kf*Y(13)*Y(50)
      X=0
      X=X+g(13)*(-1)
      X=X+g(26)*(2)
      X=X+g(50)*(-1)
      Kc=EXP(-X/T)
      REVK(325)=kf/Kc*Y(26)*Y(26)

      RETURN
      END
