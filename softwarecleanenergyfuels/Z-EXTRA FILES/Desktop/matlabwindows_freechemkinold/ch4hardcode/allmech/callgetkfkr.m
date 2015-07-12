function callgetkfkr
clc;clear;format long e;
tic

T=1000;

C=ones(1,53);
C=C';
for i=1:4225
[fwdk, revk] = getkfkr(T, C); %molelar production rate
end
fwdk';
revk';
toc
end
