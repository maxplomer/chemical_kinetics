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
%C(48)=0;
for i=1:20000
wdot=ckwc(T, C, chem);%molelar production rate
end
toc

%Q

end