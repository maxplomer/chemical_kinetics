function [DfwdkDC,DrevkDC]=getDkfkrDC(T,C,j)
addpath (fullfile(pwd,'jacobian\DkfkrDC'));
for i=1:21
   cmd=sprintf('[DfwdkDC(%g),DrevkDC(%g)]=getDkfkrDC%g(T, C,j);',i,i,i);
   eval(cmd);
end


end