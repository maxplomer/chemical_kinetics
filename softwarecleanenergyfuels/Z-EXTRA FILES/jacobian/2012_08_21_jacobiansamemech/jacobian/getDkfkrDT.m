function [DfwdkDT,DrevkDT]=getDkfkrDT(T,C)
addpath (fullfile(pwd,'jacobian\DkfkrDT'));
for i=1:21
   cmd=sprintf('[DfwdkDT(%g),DrevkDT(%g)]=getDkfkrDT%g(T, C);',i,i,i);
   eval(cmd);
end


end