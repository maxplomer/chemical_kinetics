function callkfkr
clc;clear;format long e;
tic
% add the path of the toolbox
addpath(strcat(pwd, '/chemkin'));

% load the chemkin workspace
Q=load('GRI3.0\CH4.mat');

chem = Q.chem;

T=1000;

C=ones(1,53);
P=sum(C)*Q.RU*T;
X=C/sum(C);
for i=1:22000
[fwdk, revk] = ckkfkr(P, T, X, chem); %molelar production rate
end
fwdk;
revk;
toc
end