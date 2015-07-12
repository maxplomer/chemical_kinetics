%clc;clear;
chem=ckinit;
C=[0.000000045068517
   0.000000022534258
                   1e-22
                   1e-22
                   1e-22
                   1e-22
                   1e-22
                   1e-22
   0.000000084728811];
C=C*1.0e+002;

T= 800;
format long e;

[fwdk, revk]=ckkfkrc(T,C,chem);
fwdk(9)
[fwdk, revk]=getkfkr(T,C);
fwdk(9)