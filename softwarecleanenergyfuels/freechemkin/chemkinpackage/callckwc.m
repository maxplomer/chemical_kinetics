function callckwc
clc;clear;format long e;
tic
% add the path of the toolbox
addpath(strcat(pwd, '/chemkin'));

% load the chemkin workspace
Q=load('GRI3.0\CH4.mat');

chem = Q.chem;

T=1000;

C=ones(1,53);


wdot=ckwc(T, C, chem);%molelar production rate

wdot

end