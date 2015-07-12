%reaction 1
%H+O2=O+OH                 3.547e+15 -0.406  1.6599E+4
%-units of E are cal/mol , 1 calorie = 41 840 000 erg
%double precision should have 16 sig figs
clc;clear;
format long e;
T=800;TLOG = log(T);

RU=83145100; %erg/(mol*K)
A=3.547e+15;B=-0.406;E=1.6599E+4*41840000;
kf=A*T^B*exp(-E/(RU*T))%textbook
%6.865738291409698e+009
RF = exp(3.580487857070217D1 -4.06D-1*TLOG -8.352893466155027D3/T)%DWDCT
%6.865738029975403e+009
RF = exp(3.58048786D1 -4.06D-1*TLOG -8.35289347D3/T) %ckwyp
%6.865738198128399e+009





%results
%textbook:6.865738291409698e+009
%DWDCT:   6.865738029975403e+009
%CKWYP:   6.865738198128399e+009