function [fwdk,revk]=getkfkr(P,T,X)
%number of reactions
N=21;
addpath (fullfile(pwd,'kfkr'));
for j=1:N
   cmd=sprintf('[fwdk(%g),revk(%g)]=getkfkr%g(P, T, X);',j,j,j);
   eval(cmd);
end




end