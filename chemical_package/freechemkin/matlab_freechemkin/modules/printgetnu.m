function printgetnu
global nuf nur;

[m,n] = size(nuf);
fid = fopen('hardcode/getnu.m', 'wt');
fprintf(fid, 'function [nuf,nur]=getnu\n');
%%%%%%%%%%
fprintf(fid, 'nuf=[');
for i=1:m
    for j=1:n
        fprintf(fid, ' %g ',nuf(i,j));
    end
    if i<m
        fprintf(fid, '\n     ');
    end 
end
fprintf(fid, '];\n');
%%%%%%%%%%
fprintf(fid, 'nur=[');
for i=1:m
    for j=1:n
        fprintf(fid, ' %g ',nur(i,j));
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