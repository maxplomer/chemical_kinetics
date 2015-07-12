function printgetcomposition
global composition;

[m,n] = size(composition);
fid = fopen('hardcode/getcomposition.m', 'wt');
fprintf(fid, 'function [comp]=getcomposition\n');
%%%%%%%%%%
fprintf(fid, 'comp=[');
for i=1:m
    for j=1:n
        fprintf(fid, ' %g ',composition(i,j));
    end
    if i<m
        fprintf(fid, '\n      ');
    end 
end
fprintf(fid, '];\n');

fprintf(fid, 'end');
fclose(fid);
end