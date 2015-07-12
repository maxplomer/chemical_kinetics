function [DfwdkDC,DrevkDC]=getDkfkrDC(T,C,j)

%number of reactions
N=21;
addpath (fullfile(pwd,'jacobian\DkfkrDC'));
for i=1:N
   cmd=sprintf('[DfwdkDC(%g),DrevkDC(%g)]=getDkfkrDC%g(T, C,j);',i,i,i);
   eval(cmd);
end


end