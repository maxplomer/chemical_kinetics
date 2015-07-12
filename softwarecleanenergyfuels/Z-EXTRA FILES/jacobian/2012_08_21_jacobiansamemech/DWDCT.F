      SUBROUTINE DWDCT(T, C, A, LDA)
C
C     Analytic Jacobian for species net production rate w(T, C)
C     T: input, Temperature, K, scalar
C     C: input, species mole concentration, mole/cm^3, vector, size at least KK
C     A: output, Jacobian matrix [dg/dc dg/dT], size at least KK*(KK+1)
C     LDA: input, leading dimension of matrix A, scalar, value at least KK
C
      IMPLICIT DOUBLE PRECISION (A-H, O-Z), INTEGER (I-N)
      PARAMETER (EXPARG=690.776, SMALL=1.D-300)
      DIMENSION C(*), A(LDA,*)
C
      IF (LDA .LT. 9) THEN
      WRITE(*,*) 'LDA TOO SMALL'
      STOP
      ENDIF
      DO M=1,10
      DO N=1,9
      A(N,M) = 0.D0
      ENDDO
      ENDDO
C
      CTOT = 0.D0
      DO N=1,9
      CTOT = CTOT +C(N)
      ENDDO
C
      T2 = T*T
      T3 = T2*T
      T4 = T3*T
      TLOG = LOG(T)
      TI = 1.D0/T
      TI2 = TI/T
C
      PFAC = 1.218652692702276D-2/T
      PFAC2 = PFAC*PFAC
C
C------ thermo data for species
C H2
      IF ( T .LT. 1000 ) THEN
      G0_1 = -3.298124D0*TLOG -4.124721D-4*T +1.357169166666667D-7*T2
     *  +7.896194999999999D-12*T3 -2.067436D-14*T4 -1.012521D3/T
     *  +6.592218D0
      DG_1 = -3.298124D0/T -4.124721D-4 +2.714338333333333D-7*T
     *  +2.3688585D-11*T2 -8.269744D-14*T3 +1.012521D3/T2
      ELSE
      G0_1 = -2.991423D0*TLOG -3.500322D-4*T +9.389715D-9*T2
     *  +7.692981666666667D-13*T3 -7.91376D-17*T4 -8.35034D2/T
     *  +4.346533D0
      DG_1 = -2.991423D0/T -3.500322D-4 +1.877943D-8*T +2.3078945D-12*T2
     *  -3.165504D-16*T3 +8.35034D2/T2
      ENDIF
C O2
      IF ( T .LT. 1000 ) THEN
      G0_2 = -3.212936D0*TLOG -5.63743D-4*T +9.593583333333334D-8*T2
     *  -1.0948975D-10*T3 +4.384277D-14*T4 -1.005249D3/T -2.821802D0
      DG_2 = -3.212936D0/T -5.63743D-4 +1.918716666666667D-7*T
     *  -3.2846925D-10*T2 +1.7537108D-13*T3 +1.005249D3/T2
      ELSE
      G0_2 = -3.697578D0*TLOG -3.0675985D-4*T +2.09807D-8*T2
     *  -1.479400833333333D-12*T3 +5.682175000000001D-17*T4 -1.23393D3/T
     *  +5.084119999999999D-1
      DG_2 = -3.697578D0/T -3.0675985D-4 +4.19614D-8*T -4.4382025D-12*T2
     *  +2.27287D-16*T3 +1.23393D3/T2
      ENDIF
C O
      IF ( T .LT. 1000 ) THEN
      G0_3 = -2.946429D0*TLOG +8.19083D-4*T -4.035053333333333D-7*T2
     *  +1.3357025D-10*T3 -1.945348D-14*T4 +2.914764D4/T
     *  -1.756599999999997D-2
      DG_3 = -2.946429D0/T +8.19083D-4 -8.070106666666666D-7*T
     *  +4.0071075D-10*T2 -7.781392D-14*T3 -2.914764D4/T2
      ELSE
      G0_3 = -2.54206D0*TLOG +1.377531D-5*T +5.171338333333333D-10*T2
     *  -3.792555833333334D-13*T3 +2.184026D-17*T4 +2.92308D4/T
     *  -2.378248D0
      DG_3 = -2.54206D0/T +1.377531D-5 +1.034267666666667D-9*T
     *  -1.13776675D-12*T2 +8.736104D-17*T3 -2.92308D4/T2
      ENDIF
C OH
      IF ( T .LT. 1000 ) THEN
      G0_4 = -4.12530561D0*TLOG +1.612724695D-3*T
     *  -1.087941151666667D-6*T2 +4.832113691666666D-10*T3
     *  -1.031186895D-13*T4 +3.34630913D3/T +4.81573857D0
      DG_4 = -4.12530561D0/T +1.612724695D-3 -2.175882303333333D-6*T
     *  +1.4496341075D-9*T2 -4.12474758D-13*T3 -3.34630913D3/T2
      ELSE
      G0_4 = -2.86472886D0*TLOG -5.2825224D-4*T +4.318045966666667D-8*T2
     *  -2.54348895D-12*T3 +6.6597938D-17*T4 +3.68362875D3/T
     *  -2.83691187D0
      DG_4 = -2.86472886D0/T -5.2825224D-4 +8.636091933333334D-8*T
     *  -7.63046685D-12*T2 +2.66391752D-16*T3 -3.68362875D3/T2
      ENDIF
C H2O
      IF ( T .LT. 1000 ) THEN
      G0_5 = -3.386842D0*TLOG -1.737491D-3*T +1.059116D-6*T2
     *  -5.807150833333333D-10*T3 +1.253294D-13*T4 -3.020811D4/T
     *  +7.966090000000001D-1
      DG_5 = -3.386842D0/T -1.737491D-3 +2.118232D-6*T -1.74214525D-9*T2
     *  +5.013176000000001D-13*T3 +3.020811D4/T2
      ELSE
      G0_5 = -2.672146D0*TLOG -1.5281465D-3*T +1.455043333333334D-7*T2
     *  -1.00083D-11*T3 +3.195809D-16*T4 -2.989921D4/T -4.190671D0
      DG_5 = -2.672146D0/T -1.5281465D-3 +2.910086666666667D-7*T
     *  -3.00249D-11*T2 +1.2783236D-15*T3 +2.989921D4/T2
      ENDIF
C H
      IF ( T .LT. 1000 ) THEN
      G0_6 = -2.5D0*TLOG -0.D0*T -0.D0*T2 -0.D0*T3 -0.D0*T4
     *  +2.547163D4/T +2.9601176D0
      DG_6 = -2.5D0/T -0.D0 -0.D0*T -0.D0*T2 -0.D0*T3 -2.547163D4/T2
      ELSE
      G0_6 = -2.5D0*TLOG -0.D0*T -0.D0*T2 -0.D0*T3 -0.D0*T4
     *  +2.547163D4/T +2.9601176D0
      DG_6 = -2.5D0/T -0.D0 -0.D0*T -0.D0*T2 -0.D0*T3 -2.547163D4/T2
      ENDIF
C HO2
      IF ( T .LT. 1000 ) THEN
      G0_7 = -4.30179801D0*TLOG +2.374560255D-3*T
     *  -3.526381516666666D-6*T2 +2.02303245D-9*T3
     *  -4.646125620000001D-13*T4 +2.9480804D2/T +5.851355599999999D-1
      DG_7 = -4.30179801D0/T +2.374560255D-3 -7.052763033333333D-6*T
     *  +6.06909735D-9*T2 -1.858450248D-12*T3 -2.9480804D2/T2
      ELSE
      G0_7 = -4.0172109D0*TLOG -1.119910065D-3*T
     *  +1.056096916666667D-7*T2 -9.520530833333334D-12*T3
     *  +5.39542675D-16*T4 +1.11856713D2/T +2.321087500000001D-1
      DG_7 = -4.0172109D0/T -1.119910065D-3 +2.112193833333333D-7*T
     *  -2.85615925D-11*T2 +2.1581707D-15*T3 -1.11856713D2/T2
      ENDIF
C H2O2
      IF ( T .LT. 1000 ) THEN
      G0_8 = -3.388754D0*TLOG -3.284613D-3*T +2.475021666666666D-8*T2
     *  +3.854838333333333D-10*T3 -1.2357575D-13*T4 -1.766315D4/T
     *  -3.396609D0
      DG_8 = -3.388754D0/T -3.284613D-3 +4.950043333333333D-8*T
     *  +1.1564515D-9*T2 -4.943029999999999D-13*T3 +1.766315D4/T2
      ELSE
      G0_8 = -4.573167D0*TLOG -2.168068D-3*T +2.457815D-7*T2
     *  -1.95742D-11*T3 +7.15827D-16*T4 -1.800696D4/T +4.07203D0
      DG_8 = -4.573167D0/T -2.168068D-3 +4.91563D-7*T
     *  -5.872259999999999D-11*T2 +2.863308D-15*T3 +1.800696D4/T2
      ENDIF
C------ reaction rates
C R1
      RF = EXP(3.580487857070217D1 -4.06D-1*TLOG -8.352893466155027D3/T)
      RFLGDT = -4.06D-1/T +8.352893466155027D3/T2
      G0_SUM = -G0_2 +G0_3 +G0_4 -G0_6
      DGDT = -DG_2 +DG_3 +DG_4 -DG_6
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(2)*C(6)
      WFDT = WF*RFLGDT
      WB = RB*C(3)*C(4)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(2, 10) = A(2, 10) -DWDT
      A(3, 10) = A(3, 10) +DWDT
      A(4, 10) = A(4, 10) +DWDT
      A(6, 10) = A(6, 10) -DWDT
      WFDC = RF*C(6)
      A(2, 2) = A(2, 2) -WFDC
      A(3, 2) = A(3, 2) +WFDC
      A(4, 2) = A(4, 2) +WFDC
      A(6, 2) = A(6, 2) -WFDC
      WFDC = RF*C(2)
      A(2, 6) = A(2, 6) -WFDC
      A(3, 6) = A(3, 6) +WFDC
      A(4, 6) = A(4, 6) +WFDC
      A(6, 6) = A(6, 6) -WFDC
      WBDC = RB*C(4)
      A(2, 3) = A(2, 3) +WBDC
      A(3, 3) = A(3, 3) -WBDC
      A(4, 3) = A(4, 3) -WBDC
      A(6, 3) = A(6, 3) +WBDC
      WBDC = RB*C(3)
      A(2, 4) = A(2, 4) +WBDC
      A(3, 4) = A(3, 4) -WBDC
      A(4, 4) = A(4, 4) -WBDC
      A(6, 4) = A(6, 4) +WBDC
C
C R2
      RF = EXP(1.083565163356657D1 +2.67D0*TLOG -3.165232839455095D3/T)
      RFLGDT = 2.67D0/T +3.165232839455095D3/T2
      G0_SUM = -G0_1 -G0_3 +G0_4 +G0_6
      DGDT = -DG_1 -DG_3 +DG_4 +DG_6
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(1)*C(3)
      WFDT = WF*RFLGDT
      WB = RB*C(4)*C(6)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(1, 10) = A(1, 10) -DWDT
      A(3, 10) = A(3, 10) -DWDT
      A(4, 10) = A(4, 10) +DWDT
      A(6, 10) = A(6, 10) +DWDT
      WFDC = RF*C(3)
      A(1, 1) = A(1, 1) -WFDC
      A(3, 1) = A(3, 1) -WFDC
      A(4, 1) = A(4, 1) +WFDC
      A(6, 1) = A(6, 1) +WFDC
      WFDC = RF*C(1)
      A(1, 3) = A(1, 3) -WFDC
      A(3, 3) = A(3, 3) -WFDC
      A(4, 3) = A(4, 3) +WFDC
      A(6, 3) = A(6, 3) +WFDC
      WBDC = RB*C(6)
      A(1, 4) = A(1, 4) +WBDC
      A(3, 4) = A(3, 4) +WBDC
      A(4, 4) = A(4, 4) -WBDC
      A(6, 4) = A(6, 4) -WBDC
      WBDC = RB*C(4)
      A(1, 6) = A(1, 6) +WBDC
      A(3, 6) = A(3, 6) +WBDC
      A(4, 6) = A(4, 6) -WBDC
      A(6, 6) = A(6, 6) -WBDC
C
C R3
      RF = EXP(1.919078896564844D1 +1.51D0*TLOG -1.726033170004924D3/T)
      RFLGDT = 1.51D0/T +1.726033170004924D3/T2
      G0_SUM = -G0_1 -G0_4 +G0_5 +G0_6
      DGDT = -DG_1 -DG_4 +DG_5 +DG_6
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(1)*C(4)
      WFDT = WF*RFLGDT
      WB = RB*C(5)*C(6)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(1, 10) = A(1, 10) -DWDT
      A(4, 10) = A(4, 10) -DWDT
      A(5, 10) = A(5, 10) +DWDT
      A(6, 10) = A(6, 10) +DWDT
      WFDC = RF*C(4)
      A(1, 1) = A(1, 1) -WFDC
      A(4, 1) = A(4, 1) -WFDC
      A(5, 1) = A(5, 1) +WFDC
      A(6, 1) = A(6, 1) +WFDC
      WFDC = RF*C(1)
      A(1, 4) = A(1, 4) -WFDC
      A(4, 4) = A(4, 4) -WFDC
      A(5, 4) = A(5, 4) +WFDC
      A(6, 4) = A(6, 4) +WFDC
      WBDC = RB*C(6)
      A(1, 5) = A(1, 5) +WBDC
      A(4, 5) = A(4, 5) +WBDC
      A(5, 5) = A(5, 5) -WBDC
      A(6, 5) = A(6, 5) -WBDC
      WBDC = RB*C(5)
      A(1, 6) = A(1, 6) +WBDC
      A(4, 6) = A(4, 6) +WBDC
      A(5, 6) = A(5, 6) -WBDC
      A(6, 6) = A(6, 6) -WBDC
C
C R4
      RF = EXP(1.490407251077888D1 +2.02D0*TLOG -6.743103346374924D3/T)
      RFLGDT = 2.02D0/T +6.743103346374924D3/T2
      G0_SUM = -G0_3 +G0_4 +G0_4 -G0_5
      DGDT = -DG_3 +DG_4 +DG_4 -DG_5
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(3)*C(5)
      WFDT = WF*RFLGDT
      WB = RB*C(4)*C(4)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(3, 10) = A(3, 10) -DWDT
      A(4, 10) = A(4, 10) +DWDT +DWDT
      A(5, 10) = A(5, 10) -DWDT
      WFDC = RF*C(5)
      A(3, 3) = A(3, 3) -WFDC
      A(4, 3) = A(4, 3) +WFDC +WFDC
      A(5, 3) = A(5, 3) -WFDC
      WFDC = RF*C(3)
      A(3, 5) = A(3, 5) -WFDC
      A(4, 5) = A(4, 5) +WFDC +WFDC
      A(5, 5) = A(5, 5) -WFDC
      WBDC = (RB +RB)*C(4)
      A(3, 4) = A(3, 4) +WBDC
      A(4, 4) = A(4, 4) -WBDC -WBDC
      A(5, 4) = A(5, 4) +WBDC
C
C R5
C have third body
      RF = EXP(4.527016052855837D1 -1.4D0*TLOG -5.252575576825481D4/T)
      RFLGDT = -1.4D0/T +5.252575576825481D4/T2
      G0_SUM = -G0_1 +G0_6 +G0_6
      DGDT = -DG_1 +DG_6 +DG_6
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      EQINV = EQINV/PFAC
      DGDT = DGDT + TI
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      CM = CTOT +1.5D0*C(1) +1.1D1*C(5)
      PC = CM
      PCLGDM = 1.D0/MAX(CM, SMALL)
C
      RF = RF*PC
      WF = RF*C(1)
      WFDT = WF*RFLGDT
      WFDM = WF*PCLGDM
      RB = RB*PC
      WB = RB*C(6)*C(6)
      WBDT = WB*RBLGDT
      WBDM = WB*PCLGDM
      DWDT = WFDT -WBDT
      A(1, 10) = A(1, 10) -DWDT
      A(6, 10) = A(6, 10) +DWDT +DWDT
      WFDC = RF
      A(1, 1) = A(1, 1) -WFDC
      A(6, 1) = A(6, 1) +WFDC +WFDC
      WBDC = (RB +RB)*C(6)
      A(1, 6) = A(1, 6) +WBDC
      A(6, 6) = A(6, 6) -WBDC -WBDC
      DWDM = WFDM -WBDM
      A(1, 1) = A(1, 1) -1.5D0*DWDM
      A(1, 5) = A(1, 5) -1.1D1*DWDM
      DM_1 = -DWDM
      A(6, 1) = A(6, 1) +3.D0*DWDM
      A(6, 5) = A(6, 5) +2.2D1*DWDM
      DM_6 = DWDM +DWDM
C
C R6
C have third body
      RF = EXP(3.635766453152699D1 -5.D-1*TLOG)
      RFLGDT = -5.D-1/T
      G0_SUM = +G0_2 -G0_3 -G0_3
      DGDT = +DG_2 -DG_3 -DG_3
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      EQINV = EQINV*PFAC
      DGDT = DGDT - TI
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      CM = CTOT +1.5D0*C(1) +1.1D1*C(5)
      PC = CM
      PCLGDM = 1.D0/MAX(CM, SMALL)
C
      RF = RF*PC
      WF = RF*C(3)*C(3)
      WFDT = WF*RFLGDT
      WFDM = WF*PCLGDM
      RB = RB*PC
      WB = RB*C(2)
      WBDT = WB*RBLGDT
      WBDM = WB*PCLGDM
      DWDT = WFDT -WBDT
      A(2, 10) = A(2, 10) +DWDT
      A(3, 10) = A(3, 10) -DWDT -DWDT
      WFDC = (RF +RF)*C(3)
      A(2, 3) = A(2, 3) +WFDC
      A(3, 3) = A(3, 3) -WFDC -WFDC
      WBDC = RB
      A(2, 2) = A(2, 2) -WBDC
      A(3, 2) = A(3, 2) +WBDC +WBDC
      DWDM = WFDM -WBDM
      A(2, 1) = A(2, 1) +1.5D0*DWDM
      A(2, 5) = A(2, 5) +1.1D1*DWDM
      DM_2 = DWDM
      A(3, 1) = A(3, 1) -3.D0*DWDM
      A(3, 5) = A(3, 5) -2.2D1*DWDM
      DM_3 = -DWDM -DWDM
C
C R7
C have third body
      RF = 4.714D18/T
      RFLGDT = -1.D0/T
      G0_SUM = -G0_3 +G0_4 -G0_6
      DGDT = -DG_3 +DG_4 -DG_6
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      EQINV = EQINV*PFAC
      DGDT = DGDT - TI
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      CM = CTOT +1.5D0*C(1) +1.1D1*C(5)
      PC = CM
      PCLGDM = 1.D0/MAX(CM, SMALL)
C
      RF = RF*PC
      WF = RF*C(3)*C(6)
      WFDT = WF*RFLGDT
      WFDM = WF*PCLGDM
      RB = RB*PC
      WB = RB*C(4)
      WBDT = WB*RBLGDT
      WBDM = WB*PCLGDM
      DWDT = WFDT -WBDT
      A(3, 10) = A(3, 10) -DWDT
      A(4, 10) = A(4, 10) +DWDT
      A(6, 10) = A(6, 10) -DWDT
      WFDC = RF*C(6)
      A(3, 3) = A(3, 3) -WFDC
      A(4, 3) = A(4, 3) +WFDC
      A(6, 3) = A(6, 3) -WFDC
      WFDC = RF*C(3)
      A(3, 6) = A(3, 6) -WFDC
      A(4, 6) = A(4, 6) +WFDC
      A(6, 6) = A(6, 6) -WFDC
      WBDC = RB
      A(3, 4) = A(3, 4) +WBDC
      A(4, 4) = A(4, 4) -WBDC
      A(6, 4) = A(6, 4) +WBDC
      DWDM = WFDM -WBDM
      A(3, 1) = A(3, 1) -1.5D0*DWDM
      A(3, 5) = A(3, 5) -1.1D1*DWDM
      DM_3 = DM_3 -DWDM
      A(4, 1) = A(4, 1) +1.5D0*DWDM
      A(4, 5) = A(4, 5) +1.1D1*DWDM
      DM_4 = DWDM
      A(6, 1) = A(6, 1) -1.5D0*DWDM
      A(6, 5) = A(6, 5) -1.1D1*DWDM
      DM_6 = DM_6 -DWDM
C
C R8
C have third body
      RF = 3.8D22/T2
      RFLGDT = -2.D0/T
      G0_SUM = -G0_4 +G0_5 -G0_6
      DGDT = -DG_4 +DG_5 -DG_6
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      EQINV = EQINV*PFAC
      DGDT = DGDT - TI
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      CM = CTOT +1.5D0*C(1) +1.1D1*C(5)
      PC = CM
      PCLGDM = 1.D0/MAX(CM, SMALL)
C
      RF = RF*PC
      WF = RF*C(4)*C(6)
      WFDT = WF*RFLGDT
      WFDM = WF*PCLGDM
      RB = RB*PC
      WB = RB*C(5)
      WBDT = WB*RBLGDT
      WBDM = WB*PCLGDM
      DWDT = WFDT -WBDT
      A(4, 10) = A(4, 10) -DWDT
      A(5, 10) = A(5, 10) +DWDT
      A(6, 10) = A(6, 10) -DWDT
      WFDC = RF*C(6)
      A(4, 4) = A(4, 4) -WFDC
      A(5, 4) = A(5, 4) +WFDC
      A(6, 4) = A(6, 4) -WFDC
      WFDC = RF*C(4)
      A(4, 6) = A(4, 6) -WFDC
      A(5, 6) = A(5, 6) +WFDC
      A(6, 6) = A(6, 6) -WFDC
      WBDC = RB
      A(4, 5) = A(4, 5) +WBDC
      A(5, 5) = A(5, 5) -WBDC
      A(6, 5) = A(6, 5) +WBDC
      DWDM = WFDM -WBDM
      A(4, 1) = A(4, 1) -1.5D0*DWDM
      A(4, 5) = A(4, 5) -1.1D1*DWDM
      DM_4 = DM_4 -DWDM
      A(5, 1) = A(5, 1) +1.5D0*DWDM
      A(5, 5) = A(5, 5) +1.1D1*DWDM
      DM_5 = DWDM
      A(6, 1) = A(6, 1) -1.5D0*DWDM
      A(6, 5) = A(6, 5) -1.1D1*DWDM
      DM_6 = DM_6 -DWDM
C
C R9
C have third body
C falloff: Troe 6-parameters
      RF = EXP(2.801967910572033D1 +6.D-1*TLOG)
      RFLGDT = 6.D-1/T
      G0_SUM = -G0_2 -G0_6 +G0_7
      DGDT = -DG_2 -DG_6 +DG_7
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      EQINV = EQINV*PFAC
      DGDT = DGDT - TI
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      CM = CTOT +C(1) +1.D1*C(5) -2.2D-1*C(2)
      RL = EXP(4.790267318874081D1 -1.72D0*TLOG -2.640881071774298D2/T)
      PR = RL*CM/RF
      PRLGDT = -2.32D0/T +2.640881071774298D2/T2
      PRLGDM = 1.D0/MAX(CM, SMALL)
      TEMP1 = 1.D0 +PR
      PC = PR/TEMP1
      PCLGDT = PRLGDT/TEMP1
      PCLGDM = PRLGDM/TEMP1
      PRLG = LOG(MAX(PR, SMALL))
      TEMP1 = 2.D-1*EXP(-T/1.D-30)
      TEMP2 = 8.D-1*EXP(-T/1.D30)
      FCENT = TEMP1 + TEMP2
      FCNTDT = -9.999999999999999D29*TEMP1 -9.999999999999999D-31*TEMP2
      FTLGDT = FCNTDT/FCENT
      FTLG = LOG(MAX(FCENT, SMALL))
      TEMP1 = -9.210340371976185D-1 +PRLG -6.7D-1*FTLG
      TEMP2 = 1.855883584953201D0 -1.4D-1*PRLG -1.1762D0*FTLG
      XP = TEMP1/TEMP2
      TEMP3 = 1.D0 +1.4D-1*XP
      XPDT = (TEMP3*PRLGDT +(1.1762D0*XP -6.7D-1)*FTLGDT)/TEMP2
      XPDM = TEMP3*PRLGDM/TEMP2
      TEMP1 = 1.D0 +XP*XP
      FCLG = FTLG/TEMP1
      FC = EXP(FCLG)
      TEMP2 = FCLG*(XP +XP)
      FCLGDT = (FTLGDT -TEMP2*XPDT)/TEMP1
      FCLGDM = -TEMP2*XPDM/TEMP1
      PC = PC*FC
      PCLGDT = PCLGDT +FCLGDT
      PCLGDM = PCLGDM +FCLGDM
C
      RF = RF*PC
      WF = RF*C(2)*C(6)
      WFDT = WF*(RFLGDT +PCLGDT)
      WFDM = WF*PCLGDM
      RB = RB*PC
      WB = RB*C(7)
      WBDT = WB*(RBLGDT +PCLGDT)
      WBDM = WB*PCLGDM
      DWDT = WFDT -WBDT
      A(2, 10) = A(2, 10) -DWDT
      A(6, 10) = A(6, 10) -DWDT
      A(7, 10) = A(7, 10) +DWDT
      WFDC = RF*C(6)
      A(2, 2) = A(2, 2) -WFDC
      A(6, 2) = A(6, 2) -WFDC
      A(7, 2) = A(7, 2) +WFDC
      WFDC = RF*C(2)
      A(2, 6) = A(2, 6) -WFDC
      A(6, 6) = A(6, 6) -WFDC
      A(7, 6) = A(7, 6) +WFDC
      WBDC = RB
      A(2, 7) = A(2, 7) +WBDC
      A(6, 7) = A(6, 7) +WBDC
      A(7, 7) = A(7, 7) -WBDC
      DWDM = WFDM -WBDM
      A(2, 1) = A(2, 1) -DWDM
      A(2, 5) = A(2, 5) -1.D1*DWDM
      A(2, 2) = A(2, 2) +2.2D-1*DWDM
      DM_2 = DM_2 -DWDM
      A(6, 1) = A(6, 1) -DWDM
      A(6, 5) = A(6, 5) -1.D1*DWDM
      A(6, 2) = A(6, 2) +2.2D-1*DWDM
      DM_6 = DM_6 -DWDM
      A(7, 1) = A(7, 1) +DWDM
      A(7, 5) = A(7, 5) +1.D1*DWDM
      A(7, 2) = A(7, 2) -2.2D-1*DWDM
      DM_7 = DWDM
C
C R10
      RF = EXP(3.044042381129104D1 -4.141473174676539D2/T)
      RFLGDT = 4.141473174676539D2/T2
      G0_SUM = +G0_1 +G0_2 -G0_6 -G0_7
      DGDT = +DG_1 +DG_2 -DG_6 -DG_7
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(6)*C(7)
      WFDT = WF*RFLGDT
      WB = RB*C(1)*C(2)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(1, 10) = A(1, 10) +DWDT
      A(2, 10) = A(2, 10) +DWDT
      A(6, 10) = A(6, 10) -DWDT
      A(7, 10) = A(7, 10) -DWDT
      WFDC = RF*C(7)
      A(1, 6) = A(1, 6) +WFDC
      A(2, 6) = A(2, 6) +WFDC
      A(6, 6) = A(6, 6) -WFDC
      A(7, 6) = A(7, 6) -WFDC
      WFDC = RF*C(6)
      A(1, 7) = A(1, 7) +WFDC
      A(2, 7) = A(2, 7) +WFDC
      A(6, 7) = A(6, 7) -WFDC
      A(7, 7) = A(7, 7) -WFDC
      WBDC = RB*C(2)
      A(1, 1) = A(1, 1) -WBDC
      A(2, 1) = A(2, 1) -WBDC
      A(6, 1) = A(6, 1) +WBDC
      A(7, 1) = A(7, 1) +WBDC
      WBDC = RB*C(1)
      A(1, 2) = A(1, 2) -WBDC
      A(2, 2) = A(2, 2) -WBDC
      A(6, 2) = A(6, 2) +WBDC
      A(7, 2) = A(7, 2) +WBDC
C
C R11
      RF = EXP(3.189073886371465D1 -1.484489169537763D2/T)
      RFLGDT = 1.484489169537763D2/T2
      G0_SUM = +G0_4 +G0_4 -G0_6 -G0_7
      DGDT = +DG_4 +DG_4 -DG_6 -DG_7
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(6)*C(7)
      WFDT = WF*RFLGDT
      WB = RB*C(4)*C(4)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(4, 10) = A(4, 10) +DWDT +DWDT
      A(6, 10) = A(6, 10) -DWDT
      A(7, 10) = A(7, 10) -DWDT
      WFDC = RF*C(7)
      A(4, 6) = A(4, 6) +WFDC +WFDC
      A(6, 6) = A(6, 6) -WFDC
      A(7, 6) = A(7, 6) -WFDC
      WFDC = RF*C(6)
      A(4, 7) = A(4, 7) +WFDC +WFDC
      A(6, 7) = A(6, 7) -WFDC
      A(7, 7) = A(7, 7) -WFDC
      WBDC = (RB +RB)*C(4)
      A(4, 4) = A(4, 4) -WBDC -WBDC
      A(6, 4) = A(6, 4) +WBDC
      A(7, 4) = A(7, 4) +WBDC
C
C R12
      RF = 3.25D13
      G0_SUM = +G0_2 -G0_3 +G0_4 -G0_7
      DGDT = +DG_2 -DG_3 +DG_4 -DG_7
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = DGDT
C
      WB = RB*C(2)*C(4)
      WBDT = WB*RBLGDT
      DWDT = -WBDT
      A(2, 10) = A(2, 10) +DWDT
      A(3, 10) = A(3, 10) -DWDT
      A(4, 10) = A(4, 10) +DWDT
      A(7, 10) = A(7, 10) -DWDT
      WFDC = RF*C(7)
      A(2, 3) = A(2, 3) +WFDC
      A(3, 3) = A(3, 3) -WFDC
      A(4, 3) = A(4, 3) +WFDC
      A(7, 3) = A(7, 3) -WFDC
      WFDC = RF*C(3)
      A(2, 7) = A(2, 7) +WFDC
      A(3, 7) = A(3, 7) -WFDC
      A(4, 7) = A(4, 7) +WFDC
      A(7, 7) = A(7, 7) -WFDC
      WBDC = RB*C(4)
      A(2, 2) = A(2, 2) -WBDC
      A(3, 2) = A(3, 2) +WBDC
      A(4, 2) = A(4, 2) -WBDC
      A(7, 2) = A(7, 2) +WBDC
      WBDC = RB*C(2)
      A(2, 4) = A(2, 4) -WBDC
      A(3, 4) = A(3, 4) +WBDC
      A(4, 4) = A(4, 4) -WBDC
      A(7, 4) = A(7, 4) +WBDC
C
C R13
      RF = EXP(3.099486271104694D1 +2.500986838170401D2/T)
      RFLGDT = -2.500986838170401D2/T2
      G0_SUM = +G0_2 -G0_4 +G0_5 -G0_7
      DGDT = +DG_2 -DG_4 +DG_5 -DG_7
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(4)*C(7)
      WFDT = WF*RFLGDT
      WB = RB*C(2)*C(5)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(2, 10) = A(2, 10) +DWDT
      A(4, 10) = A(4, 10) -DWDT
      A(5, 10) = A(5, 10) +DWDT
      A(7, 10) = A(7, 10) -DWDT
      WFDC = RF*C(7)
      A(2, 4) = A(2, 4) +WFDC
      A(4, 4) = A(4, 4) -WFDC
      A(5, 4) = A(5, 4) +WFDC
      A(7, 4) = A(7, 4) -WFDC
      WFDC = RF*C(4)
      A(2, 7) = A(2, 7) +WFDC
      A(4, 7) = A(4, 7) -WFDC
      A(5, 7) = A(5, 7) +WFDC
      A(7, 7) = A(7, 7) -WFDC
      WBDC = RB*C(5)
      A(2, 2) = A(2, 2) -WBDC
      A(4, 2) = A(4, 2) +WBDC
      A(5, 2) = A(5, 2) -WBDC
      A(7, 2) = A(7, 2) +WBDC
      WBDC = RB*C(2)
      A(2, 5) = A(2, 5) -WBDC
      A(4, 5) = A(4, 5) +WBDC
      A(5, 5) = A(5, 5) -WBDC
      A(7, 5) = A(7, 5) +WBDC
C
C R14
      RF = EXP(3.367127582720596D1 -6.029542111661518D3/T)
      RFLGDT = 6.029542111661518D3/T2
      G0_SUM = +G0_2 -G0_7 -G0_7 +G0_8
      DGDT = +DG_2 -DG_7 -DG_7 +DG_8
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(7)*C(7)
      WFDT = WF*RFLGDT
      WB = RB*C(2)*C(8)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(2, 10) = A(2, 10) +DWDT
      A(7, 10) = A(7, 10) -DWDT -DWDT
      A(8, 10) = A(8, 10) +DWDT
      WFDC = (RF +RF)*C(7)
      A(2, 7) = A(2, 7) +WFDC
      A(7, 7) = A(7, 7) -WFDC -WFDC
      A(8, 7) = A(8, 7) +WFDC
      WBDC = RB*C(8)
      A(2, 2) = A(2, 2) -WBDC
      A(7, 2) = A(7, 2) +WBDC +WBDC
      A(8, 2) = A(8, 2) -WBDC
      WBDC = RB*C(2)
      A(2, 8) = A(2, 8) -WBDC
      A(7, 8) = A(7, 8) +WBDC +WBDC
      A(8, 8) = A(8, 8) -WBDC
C
C R15
      RF = EXP(2.559080028740199D1 +8.198909165857211D2/T)
      RFLGDT = -8.198909165857211D2/T2
      G0_SUM = +G0_2 -G0_7 -G0_7 +G0_8
      DGDT = +DG_2 -DG_7 -DG_7 +DG_8
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(7)*C(7)
      WFDT = WF*RFLGDT
      WB = RB*C(2)*C(8)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(2, 10) = A(2, 10) +DWDT
      A(7, 10) = A(7, 10) -DWDT -DWDT
      A(8, 10) = A(8, 10) +DWDT
      WFDC = (RF +RF)*C(7)
      A(2, 7) = A(2, 7) +WFDC
      A(7, 7) = A(7, 7) -WFDC -WFDC
      A(8, 7) = A(8, 7) +WFDC
      WBDC = RB*C(8)
      A(2, 2) = A(2, 2) -WBDC
      A(7, 2) = A(7, 2) +WBDC +WBDC
      A(8, 2) = A(8, 2) -WBDC
      WBDC = RB*C(2)
      A(2, 8) = A(2, 8) -WBDC
      A(7, 8) = A(7, 8) +WBDC +WBDC
      A(8, 8) = A(8, 8) -WBDC
C
C R16
C have third body
C falloff: Troe 6-parameters
      RF = EXP(3.331833539787744D1 -2.437078321380131D4/T)
      RFLGDT = 2.437078321380131D4/T2
      G0_SUM = +G0_4 +G0_4 -G0_8
      DGDT = +DG_4 +DG_4 -DG_8
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      EQINV = EQINV/PFAC
      DGDT = DGDT + TI
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      CM = CTOT +1.5D0*C(1) +1.1D1*C(5)
      RL = EXP(3.932793341701179D1 -2.289635837761635D4/T)
      PR = RL*CM/RF
      PRLGDT = 0.D0/T -1.474424836184964D3/T2
      PRLGDM = 1.D0/MAX(CM, SMALL)
      TEMP1 = 1.D0 +PR
      PC = PR/TEMP1
      PCLGDT = PRLGDT/TEMP1
      PCLGDM = PRLGDM/TEMP1
      PRLG = LOG(MAX(PR, SMALL))
      TEMP1 = 5.D-1*EXP(-T/1.D-30)
      TEMP2 = 5.D-1*EXP(-T/1.D30)
      FCENT = TEMP1 + TEMP2
      FCNTDT = -9.999999999999999D29*TEMP1 -9.999999999999999D-31*TEMP2
      FTLGDT = FCNTDT/FCENT
      FTLG = LOG(MAX(FCENT, SMALL))
      TEMP1 = -9.210340371976185D-1 +PRLG -6.7D-1*FTLG
      TEMP2 = 1.855883584953201D0 -1.4D-1*PRLG -1.1762D0*FTLG
      XP = TEMP1/TEMP2
      TEMP3 = 1.D0 +1.4D-1*XP
      XPDT = (TEMP3*PRLGDT +(1.1762D0*XP -6.7D-1)*FTLGDT)/TEMP2
      XPDM = TEMP3*PRLGDM/TEMP2
      TEMP1 = 1.D0 +XP*XP
      FCLG = FTLG/TEMP1
      FC = EXP(FCLG)
      TEMP2 = FCLG*(XP +XP)
      FCLGDT = (FTLGDT -TEMP2*XPDT)/TEMP1
      FCLGDM = -TEMP2*XPDM/TEMP1
      PC = PC*FC
      PCLGDT = PCLGDT +FCLGDT
      PCLGDM = PCLGDM +FCLGDM
C
      RF = RF*PC
      WF = RF*C(8)
      WFDT = WF*(RFLGDT +PCLGDT)
      WFDM = WF*PCLGDM
      RB = RB*PC
      WB = RB*C(4)*C(4)
      WBDT = WB*(RBLGDT +PCLGDT)
      WBDM = WB*PCLGDM
      DWDT = WFDT -WBDT
      A(4, 10) = A(4, 10) +DWDT +DWDT
      A(8, 10) = A(8, 10) -DWDT
      WFDC = RF
      A(4, 8) = A(4, 8) +WFDC +WFDC
      A(8, 8) = A(8, 8) -WFDC
      WBDC = (RB +RB)*C(4)
      A(4, 4) = A(4, 4) -WBDC -WBDC
      A(8, 4) = A(8, 4) +WBDC
      DWDM = WFDM -WBDM
      A(4, 1) = A(4, 1) +3.D0*DWDM
      A(4, 5) = A(4, 5) +2.2D1*DWDM
      DM_4 = DM_4 +DWDM +DWDM
      A(8, 1) = A(8, 1) -1.5D0*DWDM
      A(8, 5) = A(8, 5) -1.1D1*DWDM
      DM_8 = -DWDM
C
C R17
      RF = EXP(3.081323295642516D1 -1.997770170530481D3/T)
      RFLGDT = 1.997770170530481D3/T2
      G0_SUM = +G0_4 +G0_5 -G0_6 -G0_8
      DGDT = +DG_4 +DG_5 -DG_6 -DG_8
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(6)*C(8)
      WFDT = WF*RFLGDT
      WB = RB*C(4)*C(5)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(4, 10) = A(4, 10) +DWDT
      A(5, 10) = A(5, 10) +DWDT
      A(6, 10) = A(6, 10) -DWDT
      A(8, 10) = A(8, 10) -DWDT
      WFDC = RF*C(8)
      A(4, 6) = A(4, 6) +WFDC
      A(5, 6) = A(5, 6) +WFDC
      A(6, 6) = A(6, 6) -WFDC
      A(8, 6) = A(8, 6) -WFDC
      WFDC = RF*C(6)
      A(4, 8) = A(4, 8) +WFDC
      A(5, 8) = A(5, 8) +WFDC
      A(6, 8) = A(6, 8) -WFDC
      A(8, 8) = A(8, 8) -WFDC
      WBDC = RB*C(5)
      A(4, 4) = A(4, 4) -WBDC
      A(5, 4) = A(5, 4) -WBDC
      A(6, 4) = A(6, 4) +WBDC
      A(8, 4) = A(8, 4) +WBDC
      WBDC = RB*C(4)
      A(4, 5) = A(4, 5) -WBDC
      A(5, 5) = A(5, 5) -WBDC
      A(6, 5) = A(6, 5) +WBDC
      A(8, 5) = A(8, 5) +WBDC
C
C R18
      RF = EXP(3.15063801369851D1 -4.000572507737361D3/T)
      RFLGDT = 4.000572507737361D3/T2
      G0_SUM = +G0_1 -G0_6 +G0_7 -G0_8
      DGDT = +DG_1 -DG_6 +DG_7 -DG_8
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(6)*C(8)
      WFDT = WF*RFLGDT
      WB = RB*C(1)*C(7)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(1, 10) = A(1, 10) +DWDT
      A(6, 10) = A(6, 10) -DWDT
      A(7, 10) = A(7, 10) +DWDT
      A(8, 10) = A(8, 10) -DWDT
      WFDC = RF*C(8)
      A(1, 6) = A(1, 6) +WFDC
      A(6, 6) = A(6, 6) -WFDC
      A(7, 6) = A(7, 6) +WFDC
      A(8, 6) = A(8, 6) -WFDC
      WFDC = RF*C(6)
      A(1, 8) = A(1, 8) +WFDC
      A(6, 8) = A(6, 8) -WFDC
      A(7, 8) = A(7, 8) +WFDC
      A(8, 8) = A(8, 8) -WFDC
      WBDC = RB*C(7)
      A(1, 1) = A(1, 1) -WBDC
      A(6, 1) = A(6, 1) +WBDC
      A(7, 1) = A(7, 1) -WBDC
      A(8, 1) = A(8, 1) +WBDC
      WBDC = RB*C(1)
      A(1, 7) = A(1, 7) -WBDC
      A(6, 7) = A(6, 7) +WBDC
      A(7, 7) = A(7, 7) -WBDC
      A(8, 7) = A(8, 7) +WBDC
C
C R19
      RF = EXP(1.607205171245691D1 +2.D0*TLOG -1.997770170530481D3/T)
      RFLGDT = 2.D0/T +1.997770170530481D3/T2
      G0_SUM = -G0_3 +G0_4 +G0_7 -G0_8
      DGDT = -DG_3 +DG_4 +DG_7 -DG_8
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(3)*C(8)
      WFDT = WF*RFLGDT
      WB = RB*C(4)*C(7)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(3, 10) = A(3, 10) -DWDT
      A(4, 10) = A(4, 10) +DWDT
      A(7, 10) = A(7, 10) +DWDT
      A(8, 10) = A(8, 10) -DWDT
      WFDC = RF*C(8)
      A(3, 3) = A(3, 3) -WFDC
      A(4, 3) = A(4, 3) +WFDC
      A(7, 3) = A(7, 3) +WFDC
      A(8, 3) = A(8, 3) -WFDC
      WFDC = RF*C(3)
      A(3, 8) = A(3, 8) -WFDC
      A(4, 8) = A(4, 8) +WFDC
      A(7, 8) = A(7, 8) +WFDC
      A(8, 8) = A(8, 8) -WFDC
      WBDC = RB*C(7)
      A(3, 4) = A(3, 4) +WBDC
      A(4, 4) = A(4, 4) -WBDC
      A(7, 4) = A(7, 4) -WBDC
      A(8, 4) = A(8, 4) +WBDC
      WBDC = RB*C(4)
      A(3, 7) = A(3, 7) +WBDC
      A(4, 7) = A(4, 7) -WBDC
      A(7, 7) = A(7, 7) -WBDC
      A(8, 7) = A(8, 7) +WBDC
C
C R20
      RF = 1.D12
      G0_SUM = -G0_4 +G0_5 +G0_7 -G0_8
      DGDT = -DG_4 +DG_5 +DG_7 -DG_8
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = DGDT
C
      WB = RB*C(5)*C(7)
      WBDT = WB*RBLGDT
      DWDT = -WBDT
      A(4, 10) = A(4, 10) -DWDT
      A(5, 10) = A(5, 10) +DWDT
      A(7, 10) = A(7, 10) +DWDT
      A(8, 10) = A(8, 10) -DWDT
      WFDC = RF*C(8)
      A(4, 4) = A(4, 4) -WFDC
      A(5, 4) = A(5, 4) +WFDC
      A(7, 4) = A(7, 4) +WFDC
      A(8, 4) = A(8, 4) -WFDC
      WFDC = RF*C(4)
      A(4, 8) = A(4, 8) -WFDC
      A(5, 8) = A(5, 8) +WFDC
      A(7, 8) = A(7, 8) +WFDC
      A(8, 8) = A(8, 8) -WFDC
      WBDC = RB*C(7)
      A(4, 5) = A(4, 5) +WBDC
      A(5, 5) = A(5, 5) -WBDC
      A(7, 5) = A(7, 5) -WBDC
      A(8, 5) = A(8, 5) +WBDC
      WBDC = RB*C(5)
      A(4, 7) = A(4, 7) +WBDC
      A(5, 7) = A(5, 7) -WBDC
      A(7, 7) = A(7, 7) -WBDC
      A(8, 7) = A(8, 7) +WBDC
C
C R21
      RF = EXP(3.399404921946901D1 -4.809241692634712D3/T)
      RFLGDT = 4.809241692634712D3/T2
      G0_SUM = -G0_4 +G0_5 +G0_7 -G0_8
      DGDT = -DG_4 +DG_5 +DG_7 -DG_8
      EQINV = EXP(MIN(G0_SUM, EXPARG))
      RB = RF*EQINV
      RBLGDT = RFLGDT +DGDT
C
      WF = RF*C(4)*C(8)
      WFDT = WF*RFLGDT
      WB = RB*C(5)*C(7)
      WBDT = WB*RBLGDT
      DWDT = WFDT -WBDT
      A(4, 10) = A(4, 10) -DWDT
      A(5, 10) = A(5, 10) +DWDT
      A(7, 10) = A(7, 10) +DWDT
      A(8, 10) = A(8, 10) -DWDT
      WFDC = RF*C(8)
      A(4, 4) = A(4, 4) -WFDC
      A(5, 4) = A(5, 4) +WFDC
      A(7, 4) = A(7, 4) +WFDC
      A(8, 4) = A(8, 4) -WFDC
      WFDC = RF*C(4)
      A(4, 8) = A(4, 8) -WFDC
      A(5, 8) = A(5, 8) +WFDC
      A(7, 8) = A(7, 8) +WFDC
      A(8, 8) = A(8, 8) -WFDC
      WBDC = RB*C(7)
      A(4, 5) = A(4, 5) +WBDC
      A(5, 5) = A(5, 5) -WBDC
      A(7, 5) = A(7, 5) -WBDC
      A(8, 5) = A(8, 5) +WBDC
      WBDC = RB*C(5)
      A(4, 7) = A(4, 7) +WBDC
      A(5, 7) = A(5, 7) -WBDC
      A(7, 7) = A(7, 7) -WBDC
      A(8, 7) = A(8, 7) +WBDC
C
      DO N=1,9
      A(1, N) = A(1, N) + DM_1
      A(2, N) = A(2, N) + DM_2
      A(3, N) = A(3, N) + DM_3
      A(4, N) = A(4, N) + DM_4
      A(5, N) = A(5, N) + DM_5
      A(6, N) = A(6, N) + DM_6
      A(7, N) = A(7, N) + DM_7
      A(8, N) = A(8, N) + DM_8
      ENDDO
      END
