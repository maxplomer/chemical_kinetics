clc;clear;
cd ..
chem=ckinit;

C =[4.506851665737935e-006
    2.253425824082086e-006
    1.480701008485403e-019
    2.784516914365540e-020
    1.287650490261801e-015
    2.364418123526601e-018
    1.040676019580956e-014
    1.964644709391239e-015
    8.472881147486237e-006];


T= 8.000000012609640e+002;

J=getjacobian(T,C,chem);
J

TLOG=log(T);
T2 = T*T;
T3 = T2*T;
T4 = T3*T;
TI = 1/T;
TI2 = TI/T;
%thermo data
% H2
   if T<1000
      G0_1 = -3.298124D0*TLOG -4.124721D-4*T +1.357169166666667D-7*T2 ...
       +7.896194999999999D-12*T3 -2.067436D-14*T4 -1.012521D3/T ...
       +6.592218D0;
      DG_1 = -3.298124D0/T -4.124721D-4 +2.714338333333333D-7*T ...
       +2.3688585D-11*T2 -8.269744D-14*T3 +1.012521D3/T2;
      else
      G0_1 = -2.991423D0*TLOG -3.500322D-4*T +9.389715D-9*T2 ...
       +7.692981666666667D-13*T3 -7.91376D-17*T4 -8.35034D2/T ...
       +4.346533D0;
      DG_1 = -2.991423D0/T -3.500322D-4 +1.877943D-8*T +2.3078945D-12*T2 ...
       -3.165504D-16*T3 +8.35034D2/T2;
      end
% O2
      if T<1000
      G0_2 = -3.212936D0*TLOG -5.63743D-4*T +9.593583333333334D-8*T2 ...
       -1.0948975D-10*T3 +4.384277D-14*T4 -1.005249D3/T -2.821802D0;
      DG_2 = -3.212936D0/T -5.63743D-4 +1.918716666666667D-7*T ...
       -3.2846925D-10*T2 +1.7537108D-13*T3 +1.005249D3/T2;
      else
      G0_2 = -3.697578D0*TLOG -3.0675985D-4*T +2.09807D-8*T2 ...
       -1.479400833333333D-12*T3 +5.682175000000001D-17*T4 -1.23393D3/T ...
       +5.084119999999999D-1;
      DG_2 = -3.697578D0/T -3.0675985D-4 +4.19614D-8*T -4.4382025D-12*T2 ...
       +2.27287D-16*T3 +1.23393D3/T2;
      end
% O
      if T<1000
      G0_3 = -2.946429D0*TLOG +8.19083D-4*T -4.035053333333333D-7*T2 ...
       +1.3357025D-10*T3 -1.945348D-14*T4 +2.914764D4/T ...
       -1.756599999999997D-2;
      DG_3 = -2.946429D0/T +8.19083D-4 -8.070106666666666D-7*T ...
       +4.0071075D-10*T2 -7.781392D-14*T3 -2.914764D4/T2;
      else
      G0_3 = -2.54206D0*TLOG +1.377531D-5*T +5.171338333333333D-10*T2 ...
       -3.792555833333334D-13*T3 +2.184026D-17*T4 +2.92308D4/T ...
       -2.378248D0;
      DG_3 = -2.54206D0/T +1.377531D-5 +1.034267666666667D-9*T ...
       -1.13776675D-12*T2 +8.736104D-17*T3 -2.92308D4/T2;
      end
% OH
      if T<1000
      G0_4 = -4.12530561D0*TLOG +1.612724695D-3*T ...
       -1.087941151666667D-6*T2 +4.832113691666666D-10*T3 ...
       -1.031186895D-13*T4 +3.34630913D3/T +4.81573857D0;
      DG_4 = -4.12530561D0/T +1.612724695D-3 -2.175882303333333D-6*T ...
       +1.4496341075D-9*T2 -4.12474758D-13*T3 -3.34630913D3/T2;
      else
      G0_4 = -2.86472886D0*TLOG -5.2825224D-4*T +4.318045966666667D-8*T2 ...
       -2.54348895D-12*T3 +6.6597938D-17*T4 +3.68362875D3/T ...
       -2.83691187D0;
      DG_4 = -2.86472886D0/T -5.2825224D-4 +8.636091933333334D-8*T ...
       -7.63046685D-12*T2 +2.66391752D-16*T3 -3.68362875D3/T2;
      end
% H2O
      if T<1000
      G0_5 = -3.386842D0*TLOG -1.737491D-3*T +1.059116D-6*T2 ...
       -5.807150833333333D-10*T3 +1.253294D-13*T4 -3.020811D4/T ...
       +7.966090000000001D-1;
      DG_5 = -3.386842D0/T -1.737491D-3 +2.118232D-6*T -1.74214525D-9*T2 ...
       +5.013176000000001D-13*T3 +3.020811D4/T2;
      else
      G0_5 = -2.672146D0*TLOG -1.5281465D-3*T +1.455043333333334D-7*T2 ...
       -1.00083D-11*T3 +3.195809D-16*T4 -2.989921D4/T -4.190671D0;
      DG_5 = -2.672146D0/T -1.5281465D-3 +2.910086666666667D-7*T ...
       -3.00249D-11*T2 +1.2783236D-15*T3 +2.989921D4/T2;
      end
% H
      if T<1000
      G0_6 = -2.5D0*TLOG -0.D0*T -0.D0*T2 -0.D0*T3 -0.D0*T4 ...
       +2.547163D4/T +2.9601176D0;
      DG_6 = -2.5D0/T -0.D0 -0.D0*T -0.D0*T2 -0.D0*T3 -2.547163D4/T2;
      else
      G0_6 = -2.5D0*TLOG -0.D0*T -0.D0*T2 -0.D0*T3 -0.D0*T4 ...
       +2.547163D4/T +2.9601176D0;
      DG_6 = -2.5D0/T -0.D0 -0.D0*T -0.D0*T2 -0.D0*T3 -2.547163D4/T2;
      end
% HO2
      if T<1000
      G0_7 = -4.30179801D0*TLOG +2.374560255D-3*T ...
       -3.526381516666666D-6*T2 +2.02303245D-9*T3 ...
       -4.646125620000001D-13*T4 +2.9480804D2/T +5.851355599999999D-1;
      DG_7 = -4.30179801D0/T +2.374560255D-3 -7.052763033333333D-6*T ...
       +6.06909735D-9*T2 -1.858450248D-12*T3 -2.9480804D2/T2;
      else
      G0_7 = -4.0172109D0*TLOG -1.119910065D-3*T ...
       +1.056096916666667D-7*T2 -9.520530833333334D-12*T3 ...
       +5.39542675D-16*T4 +1.11856713D2/T +2.321087500000001D-1;
      DG_7 = -4.0172109D0/T -1.119910065D-3 +2.112193833333333D-7*T ...
       -2.85615925D-11*T2 +2.1581707D-15*T3 -1.11856713D2/T2;
      end
% H2O2
      if T<1000
      G0_8 = -3.388754D0*TLOG -3.284613D-3*T +2.475021666666666D-8*T2 ...
       +3.854838333333333D-10*T3 -1.2357575D-13*T4 -1.766315D4/T ...
       -3.396609D0;
      DG_8 = -3.388754D0/T -3.284613D-3 +4.950043333333333D-8*T ...
       +1.1564515D-9*T2 -4.943029999999999D-13*T3 +1.766315D4/T2;
      else
      G0_8 = -4.573167D0*TLOG -2.168068D-3*T +2.457815D-7*T2 ...
       -1.95742D-11*T3 +7.15827D-16*T4 -1.800696D4/T +4.07203D0;
      DG_8 = -4.573167D0/T -2.168068D-3 +4.91563D-7*T ...
       -5.872259999999999D-11*T2 +2.863308D-15*T3 +1.800696D4/T2;
      end
      
      
      RF = exp(3.580487857070217D1 -4.06D-1*TLOG -8.352893466155027D3/T);%%%%%%%%%
      RFLGDT = -4.06D-1/T +8.352893466155027D3/T2;
      G0_SUM = -G0_2 +G0_3 +G0_4 -G0_6;
      DGDT = -DG_2 +DG_3 +DG_4 -DG_6;
      EQINV = exp(G0_SUM);%%%%%%%%%%%%
      RB = RF*EQINV;
      RBLGDT = RFLGDT +DGDT;

      WF = RF*C(2)*C(6);
      WFDT = WF*RFLGDT;
      WB = RB*C(3)*C(4);
      WBDT = WB*RBLGDT;
      DWDT = WFDT -WBDT;
      
      A=zeros(10,10);%%%%%%%%%%%%
      
      A(2, 10) = A(2, 10) -DWDT;
      A(3, 10) = A(3, 10) +DWDT;
      A(4, 10) = A(4, 10) +DWDT;
      A(6, 10) = A(6, 10) -DWDT;
      WFDC = RF*C(6);
      A(2, 2) = A(2, 2) -WFDC;
      A(3, 2) = A(3, 2) +WFDC;
      A(4, 2) = A(4, 2) +WFDC;
      A(6, 2) = A(6, 2) -WFDC;
      WFDC = RF*C(2);
      A(2, 6) = A(2, 6) -WFDC;
      A(3, 6) = A(3, 6) +WFDC;
      A(4, 6) = A(4, 6) +WFDC;
      A(6, 6) = A(6, 6) -WFDC;
      WBDC = RB*C(4);
      A(2, 3) = A(2, 3) +WBDC;
      A(3, 3) = A(3, 3) -WBDC;
      A(4, 3) = A(4, 3) -WBDC;
      A(6, 3) = A(6, 3) +WBDC;
      WBDC = RB*C(3);
      A(2, 4) = A(2, 4) +WBDC;
      A(3, 4) = A(3, 4) -WBDC;
      A(4, 4) = A(4, 4) -WBDC;
      A(6, 4) = A(6, 4) +WBDC;
      
      A