function [DfwdkDT,DrevkDT]=getDkfkrDT(T,C)

%number of reactions
N=21;
addpath (fullfile(pwd,'jacobian\DkfkrDT'));
for i=1:N
   cmd=sprintf('[DfwdkDT(%g),DrevkDT(%g)]=getDkfkrDT%g(T, C);',i,i,i);
   eval(cmd);
end


end