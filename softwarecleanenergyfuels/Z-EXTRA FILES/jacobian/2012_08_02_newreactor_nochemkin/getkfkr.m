function [fwdk,revk]=getkfkr(T,Y)
%number of reactions
N=21;
addpath (fullfile(pwd,'kfkr'));
for j=1:N
   cmd=sprintf('[fwdk(%g),revk(%g)]=getkfkr%g(T, Y);',j,j,j);
   eval(cmd);
end




end