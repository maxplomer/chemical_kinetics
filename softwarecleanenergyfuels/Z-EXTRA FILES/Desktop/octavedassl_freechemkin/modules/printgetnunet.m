function printgetnunet
global nuf nur;
nu=nur-nuf;
[m,n] = size(nuf);
fid = fopen('hardcode/getnunet.m', 'wt');
fprintf(fid, 'function [nu]=getnunet\n');
%%%%%%%%%%
fprintf(fid, 'nu=[');
for i=1:m
    for j=1:n
        fprintf(fid, ' %g ',nu(i,j));
    end
    if i<m
        fprintf(fid, '\n     ');
    end 
end
fprintf(fid, '];\n');
%%%%%%%%%%

fprintf(fid, 'end');
fclose(fid);
end