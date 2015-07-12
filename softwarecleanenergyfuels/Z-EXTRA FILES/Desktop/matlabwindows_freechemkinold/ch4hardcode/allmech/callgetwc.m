function callgetwc
clc;clear;format long e;
tic

T=1000;

C=ones(1,53);
C=C';
for i=1:3600
wdot=getwc(T, C); %molelar production rate
end
toc



end