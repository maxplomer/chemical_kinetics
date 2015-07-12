function callkfkr
clc;clear;format long e;

% add the path of the toolbox
addpath(strcat(pwd, '/chemkin'));

% load the chemkin workspace
Q=load('GRI3.0\CH4.mat');

chem = Q.chem;

T=1000;

C=ones(1,53);
P=sum(C)*Q.RU*T;
X=C/sum(C);

[fwdk, revk] = ckkfkr(P, T, X, chem); %molelar production rate

fwdk
revk

end