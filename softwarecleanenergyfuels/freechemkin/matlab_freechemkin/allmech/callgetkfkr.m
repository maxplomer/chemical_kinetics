function callgetkfkr
clc;clear;format long e;

global N commonT;
[N,commonT]=getthermodata;
T=1000;

C=ones(1,53);
C=C';

[fwdk, revk] = getkfkr(T, C); %molelar production rate

fwdk'
revk'

end
