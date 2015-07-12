function callgetwc
clc;clear;format long e;

global N commonT;
[N,commonT]=getthermodata;
T=1000;

C=ones(1,53);
C=C';

wdot=getwc(T, C); %molelar production rate


wdot



end